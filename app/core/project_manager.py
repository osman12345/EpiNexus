"""
Project Manager - Save/Load Analysis Sessions

Handles:
- Project creation and metadata
- Session state serialization
- Results export/import
- Version tracking
"""

import json
import pickle
import hashlib
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import pandas as pd
import numpy as np


@dataclass
class ProjectMetadata:
    """Project metadata structure."""
    name: str
    description: str
    created_at: str
    modified_at: str
    version: str = "1.0.0"
    author: str = ""
    genome: str = "hg38"
    marks: List[str] = None
    samples: int = 0
    status: str = "active"

    def __post_init__(self):
        if self.marks is None:
            self.marks = []


class ProjectManager:
    """Manage EpiNexus analysis projects."""

    PROJECT_VERSION = "1.0.0"

    def __init__(self, projects_dir: str = "projects"):
        self.projects_dir = Path(projects_dir)
        self.projects_dir.mkdir(parents=True, exist_ok=True)

    def create_project(
        self,
        name: str,
        description: str = "",
        author: str = "",
        genome: str = "hg38",
        marks: List[str] = None
    ) -> str:
        """Create a new project."""

        # Generate project ID
        project_id = self._generate_id(name)
        project_path = self.projects_dir / project_id
        project_path.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        (project_path / "data").mkdir(exist_ok=True)
        (project_path / "results").mkdir(exist_ok=True)
        (project_path / "exports").mkdir(exist_ok=True)

        # Create metadata
        now = datetime.now().isoformat()
        metadata = ProjectMetadata(
            name=name,
            description=description,
            created_at=now,
            modified_at=now,
            version=self.PROJECT_VERSION,
            author=author,
            genome=genome,
            marks=marks or []
        )

        # Save metadata
        self._save_metadata(project_path, metadata)

        return project_id

    def save_session(
        self,
        project_id: str,
        session_state: Dict[str, Any],
        checkpoint_name: str = None
    ) -> str:
        """Save current session state to project."""

        project_path = self.projects_dir / project_id
        if not project_path.exists():
            raise ValueError(f"Project {project_id} not found")

        # Generate checkpoint name if not provided
        if checkpoint_name is None:
            checkpoint_name = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Prepare session data
        session_data = {
            "checkpoint_name": checkpoint_name,
            "saved_at": datetime.now().isoformat(),
            "state": {}
        }

        # Serialize session state
        for key, value in session_state.items():
            session_data["state"][key] = self._serialize_value(value)

        # Save session file
        session_file = project_path / f"session_{checkpoint_name}.json"
        with open(session_file, 'w') as f:
            json.dump(session_data, f, indent=2, default=str)

        # Save DataFrames separately as parquet
        for key, value in session_state.items():
            if isinstance(value, pd.DataFrame):
                df_path = project_path / "data" / f"{key}_{checkpoint_name}.parquet"
                value.to_parquet(df_path)

        # Update metadata
        metadata = self._load_metadata(project_path)
        metadata.modified_at = datetime.now().isoformat()
        self._save_metadata(project_path, metadata)

        return checkpoint_name

    def load_session(
        self,
        project_id: str,
        checkpoint_name: str = None
    ) -> Dict[str, Any]:
        """Load a saved session state."""

        project_path = self.projects_dir / project_id
        if not project_path.exists():
            raise ValueError(f"Project {project_id} not found")

        # Find latest checkpoint if not specified
        if checkpoint_name is None:
            checkpoints = self.list_checkpoints(project_id)
            if not checkpoints:
                return {}
            checkpoint_name = checkpoints[-1]["name"]

        # Load session file
        session_file = project_path / f"session_{checkpoint_name}.json"
        if not session_file.exists():
            raise ValueError(f"Checkpoint {checkpoint_name} not found")

        with open(session_file, 'r') as f:
            session_data = json.load(f)

        # Deserialize state
        state = {}
        for key, value in session_data["state"].items():
            state[key] = self._deserialize_value(value, project_path, checkpoint_name)

        return state

    def list_projects(self) -> List[Dict[str, Any]]:
        """List all projects."""
        projects = []

        for project_path in self.projects_dir.iterdir():
            if project_path.is_dir():
                try:
                    metadata = self._load_metadata(project_path)
                    projects.append({
                        "id": project_path.name,
                        **asdict(metadata)
                    })
                except Exception:
                    continue

        # Sort by modified date
        projects.sort(key=lambda x: x["modified_at"], reverse=True)
        return projects

    def list_checkpoints(self, project_id: str) -> List[Dict[str, Any]]:
        """List all checkpoints for a project."""
        project_path = self.projects_dir / project_id
        checkpoints = []

        for f in project_path.glob("session_*.json"):
            checkpoint_name = f.stem.replace("session_", "")
            with open(f, 'r') as file:
                data = json.load(file)
                checkpoints.append({
                    "name": checkpoint_name,
                    "saved_at": data.get("saved_at", ""),
                    "file": str(f)
                })

        checkpoints.sort(key=lambda x: x["saved_at"])
        return checkpoints

    def export_project(self, project_id: str, output_path: str) -> str:
        """Export project as a zip archive."""
        import shutil

        project_path = self.projects_dir / project_id
        if not project_path.exists():
            raise ValueError(f"Project {project_id} not found")

        # Create zip archive
        output_file = shutil.make_archive(
            output_path.replace('.zip', ''),
            'zip',
            project_path
        )

        return output_file

    def import_project(self, zip_path: str) -> str:
        """Import a project from zip archive."""
        import shutil
        import zipfile

        # Extract to temp location
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Read metadata to get project name
            temp_dir = self.projects_dir / "_temp_import"
            zip_ref.extractall(temp_dir)

            metadata = self._load_metadata(temp_dir)
            project_id = self._generate_id(metadata.name)

            # Move to final location
            final_path = self.projects_dir / project_id
            if final_path.exists():
                shutil.rmtree(final_path)
            temp_dir.rename(final_path)

        return project_id

    def delete_project(self, project_id: str) -> bool:
        """Delete a project."""
        import shutil

        project_path = self.projects_dir / project_id
        if project_path.exists():
            shutil.rmtree(project_path)
            return True
        return False

    def _generate_id(self, name: str) -> str:
        """Generate unique project ID."""
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        hash_str = hashlib.md5(f"{name}{timestamp}".encode()).hexdigest()[:8]
        safe_name = "".join(c if c.isalnum() else "_" for c in name.lower())[:20]
        return f"{safe_name}_{hash_str}"

    def _save_metadata(self, project_path: Path, metadata: ProjectMetadata):
        """Save project metadata."""
        with open(project_path / "metadata.json", 'w') as f:
            json.dump(asdict(metadata), f, indent=2)

    def _load_metadata(self, project_path: Path) -> ProjectMetadata:
        """Load project metadata."""
        with open(project_path / "metadata.json", 'r') as f:
            data = json.load(f)
            return ProjectMetadata(**data)

    def _serialize_value(self, value: Any) -> Dict[str, Any]:
        """Serialize a value for JSON storage."""
        if isinstance(value, pd.DataFrame):
            return {"_type": "DataFrame", "columns": list(value.columns), "shape": list(value.shape)}
        elif isinstance(value, np.ndarray):
            return {"_type": "ndarray", "shape": list(value.shape), "data": value.tolist()}
        elif isinstance(value, (list, dict, str, int, float, bool, type(None))):
            return {"_type": "primitive", "data": value}
        else:
            return {"_type": "unknown", "repr": repr(value)}

    def _deserialize_value(self, value: Dict, project_path: Path, checkpoint: str) -> Any:
        """Deserialize a value from JSON storage."""
        val_type = value.get("_type", "primitive")

        if val_type == "DataFrame":
            # Load from parquet if exists
            # For now return empty DataFrame with correct columns
            return pd.DataFrame(columns=value.get("columns", []))
        elif val_type == "ndarray":
            return np.array(value["data"])
        elif val_type == "primitive":
            return value["data"]
        else:
            return None
