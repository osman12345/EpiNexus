"""
Workflow Manager - Track, export, and reproduce analysis workflows.

Captures analysis steps with full parameters, tool versions, and commands
for reproducibility and publication-ready methods generation.

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
import json
import hashlib
import subprocess
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional
import uuid

# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass
class WorkflowStep:
    """Represents a single analysis step in the workflow."""

    step_id: str
    step_type: str  # alignment, peak_calling, differential, methylation, annotation
    timestamp: str
    status: str  # completed, failed, running
    parameters: dict = field(default_factory=dict)
    inputs: dict = field(default_factory=dict)
    outputs: dict = field(default_factory=dict)
    tool_versions: dict = field(default_factory=dict)
    runtime_seconds: float = 0.0
    command_line: str = ""
    output_metadata: dict = field(default_factory=dict)
    error_message: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "WorkflowStep":
        """Create from dictionary."""
        return cls(**data)


@dataclass
class ValidationResult:
    """Result of workflow validation."""

    is_valid: bool
    schema_valid: bool = True
    files_found: dict = field(default_factory=dict)  # filename -> exists
    files_missing: list = field(default_factory=list)
    checksum_matches: dict = field(default_factory=dict)  # filename -> matches
    tool_versions_compatible: bool = True
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)


# =============================================================================
# WORKFLOW MANAGER
# =============================================================================


class WorkflowManager:
    """
    Manages workflow tracking, export, and import for EpiNexus.

    All methods are class methods operating on Streamlit session state.
    """

    SCHEMA_VERSION = "2.0"
    STEP_TYPES = [
        "alignment",
        "peak_calling",
        "differential",
        "annotation",
        "methylation_alignment",
        "methylation_calling",
        "dmr_analysis",
        "super_enhancer",
        "tf_analysis",
        "motif_analysis",
    ]

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------

    @classmethod
    def init_session_state(cls):
        """Initialize workflow-related session state."""
        if "workflow_steps" not in st.session_state:
            st.session_state.workflow_steps = []
        if "workflow_metadata" not in st.session_state:
            st.session_state.workflow_metadata = {}

    # -------------------------------------------------------------------------
    # Step Recording
    # -------------------------------------------------------------------------

    @classmethod
    def record_step(
        cls,
        step_type: str,
        parameters: dict,
        inputs: dict = None,
        outputs: dict = None,
        tool_versions: dict = None,
        runtime_seconds: float = 0.0,
        command_line: str = "",
        output_metadata: dict = None,
        status: str = "completed",
    ) -> WorkflowStep:
        """
        Record a completed analysis step.

        Args:
            step_type: Type of analysis (alignment, peak_calling, etc.)
            parameters: Analysis parameters used
            inputs: Input file paths
            outputs: Output file paths
            tool_versions: Version of tools used
            runtime_seconds: How long the step took
            command_line: Exact command executed (for transparency)
            output_metadata: Summary stats about outputs
            status: Step status (completed, failed)

        Returns:
            The recorded WorkflowStep
        """
        cls.init_session_state()

        step = WorkflowStep(
            step_id=f"{step_type}_{uuid.uuid4().hex[:8]}",
            step_type=step_type,
            timestamp=datetime.now().isoformat(),
            status=status,
            parameters=parameters or {},
            inputs=inputs or {},
            outputs=outputs or {},
            tool_versions=tool_versions or {},
            runtime_seconds=runtime_seconds,
            command_line=command_line,
            output_metadata=output_metadata or {},
        )

        st.session_state.workflow_steps.append(step.to_dict())
        return step

    @classmethod
    def record_failed_step(
        cls, step_type: str, parameters: dict, error_message: str, inputs: dict = None, command_line: str = ""
    ) -> WorkflowStep:
        """Record a failed analysis step."""
        cls.init_session_state()

        step = WorkflowStep(
            step_id=f"{step_type}_{uuid.uuid4().hex[:8]}",
            step_type=step_type,
            timestamp=datetime.now().isoformat(),
            status="failed",
            parameters=parameters or {},
            inputs=inputs or {},
            command_line=command_line,
            error_message=error_message,
        )

        st.session_state.workflow_steps.append(step.to_dict())
        return step

    # -------------------------------------------------------------------------
    # Workflow History
    # -------------------------------------------------------------------------

    @classmethod
    def get_workflow_history(cls) -> list:
        """Get all recorded workflow steps."""
        cls.init_session_state()
        return st.session_state.workflow_steps

    @classmethod
    def get_steps_by_type(cls, step_type: str) -> list:
        """Get all steps of a specific type."""
        return [s for s in cls.get_workflow_history() if s["step_type"] == step_type]

    @classmethod
    def get_latest_step(cls, step_type: str = None) -> Optional[dict]:
        """Get the most recent step, optionally filtered by type."""
        history = cls.get_workflow_history()
        if step_type:
            history = [s for s in history if s["step_type"] == step_type]
        return history[-1] if history else None

    @classmethod
    def clear_history(cls):
        """Clear all workflow history."""
        st.session_state.workflow_steps = []

    # -------------------------------------------------------------------------
    # Export
    # -------------------------------------------------------------------------

    @classmethod
    def export_workflow(cls, include_checksums: bool = True, include_commands: bool = True) -> dict:
        """
        Export the complete workflow as a dictionary.

        Args:
            include_checksums: Calculate SHA256 for output files
            include_commands: Include command line strings

        Returns:
            Complete workflow dictionary ready for JSON serialization
        """
        cls.init_session_state()

        # Get project info from session state
        project = {
            "assay_type": st.session_state.get("assay_type", "Unknown"),
            "target_type": st.session_state.get("target_type", "Unknown"),
            "genome": st.session_state.get("selected_genome", "Unknown"),
            "samples": st.session_state.get("samples", []),
        }

        # Process steps
        steps = []
        for step_dict in cls.get_workflow_history():
            step = step_dict.copy()
            if not include_commands:
                step["command_line"] = ""
            steps.append(step)

        # Calculate checksums if requested
        checksums = {}
        if include_checksums:
            for step in steps:
                for name, filepath in step.get("outputs", {}).items():
                    if filepath and Path(filepath).exists():
                        checksums[filepath] = cls._calculate_checksum(filepath)

        export_data = {
            "format_version": cls.SCHEMA_VERSION,
            "metadata": {
                "name": st.session_state.get("project_name", "Untitled"),
                "description": st.session_state.get("project_description", ""),
                "exported_at": datetime.now().isoformat(),
                "epinexus_version": "0.2.0",
                "genome": project["genome"],
            },
            "project": project,
            "workflow_steps": steps,
            "checksums": checksums,
        }

        return export_data

    @classmethod
    def export_to_json(cls, include_checksums: bool = True, include_commands: bool = True) -> str:
        """Export workflow as JSON string."""
        data = cls.export_workflow(include_checksums, include_commands)
        return json.dumps(data, indent=2, default=str)

    # -------------------------------------------------------------------------
    # Import & Validation
    # -------------------------------------------------------------------------

    @classmethod
    def import_workflow(cls, data: dict) -> ValidationResult:
        """
        Import and validate a workflow file.

        Args:
            data: Parsed workflow dictionary

        Returns:
            ValidationResult with status and any issues found
        """
        result = ValidationResult(is_valid=True)

        # Check schema version
        version = data.get("format_version", "1.0")
        if version not in ["1.0", "2.0"]:
            result.warnings.append(f"Unknown schema version: {version}")

        # Validate required fields
        required = ["metadata", "project", "workflow_steps"]
        for field in required:  # noqa: F402
            if field not in data:
                result.errors.append(f"Missing required field: {field}")
                result.schema_valid = False
                result.is_valid = False

        if not result.schema_valid:
            return result

        # Check input files exist
        for step in data.get("workflow_steps", []):
            for name, filepath in step.get("inputs", {}).items():
                if filepath:
                    exists = Path(filepath).exists()
                    result.files_found[filepath] = exists
                    if not exists:
                        result.files_missing.append(filepath)

        # Validate checksums if present
        checksums = data.get("checksums", {})
        for filepath, expected_hash in checksums.items():
            if Path(filepath).exists():
                actual_hash = cls._calculate_checksum(filepath)
                matches = actual_hash == expected_hash
                result.checksum_matches[filepath] = matches
                if not matches:
                    result.warnings.append(f"Checksum mismatch for {filepath}")

        # Check if any critical files are missing
        if result.files_missing:
            result.warnings.append(f"{len(result.files_missing)} input file(s) not found")

        return result

    @classmethod
    def load_workflow_from_file(cls, file_content: str) -> tuple:
        """
        Load workflow from JSON string.

        Returns:
            Tuple of (workflow_data, ValidationResult)
        """
        try:
            data = json.loads(file_content)
            validation = cls.import_workflow(data)
            return data, validation
        except json.JSONDecodeError as e:
            return None, ValidationResult(is_valid=False, schema_valid=False, errors=[f"Invalid JSON: {str(e)}"])

    # -------------------------------------------------------------------------
    # Methods Text Generation
    # -------------------------------------------------------------------------

    @classmethod
    def generate_methods_text(cls) -> str:
        """
        Generate publication-ready methods section from workflow history.

        Returns:
            Markdown-formatted methods text
        """
        history = cls.get_workflow_history()
        if not history:
            return "No analysis steps recorded."

        sections = []
        genome = st.session_state.get("selected_genome", "reference genome")

        # Group steps by type
        alignment_steps = [s for s in history if s["step_type"] == "alignment"]
        peak_steps = [s for s in history if s["step_type"] == "peak_calling"]
        diff_steps = [s for s in history if s["step_type"] == "differential"]
        meth_steps = [s for s in history if "methylation" in s["step_type"]]

        # Alignment section
        if alignment_steps:
            step = alignment_steps[-1]
            params = step["parameters"]
            versions = step.get("tool_versions", {})
            aligner = params.get("aligner", "Bowtie2")
            aligner_version = versions.get(aligner.lower(), "")

            sections.append(f"""### Read Alignment

Reads were aligned to the {genome} reference genome using {aligner}{"" if not aligner_version else f" (v{aligner_version})"} with {params.get("threads", "default")} threads. {"Paired-end" if params.get("paired", True) else "Single-end"} alignment was performed with default parameters{" and additional options: " + params.get("extra_params", "") if params.get("extra_params") else ""}.""")

        # Peak calling section
        if peak_steps:
            step = peak_steps[-1]
            params = step["parameters"]
            versions = step.get("tool_versions", {})
            caller = params.get("caller", "MACS2")
            caller_version = versions.get(caller.lower(), "")
            meta = step.get("output_metadata", {})

            sections.append(f"""### Peak Calling

Peaks were called using {caller}{"" if not caller_version else f" (v{caller_version})"} with q-value threshold of {params.get("qvalue", 0.05)}. {"Broad peak calling was enabled for histone marks." if params.get("broad", False) else "Narrow peak calling was used."} {"Summit calling was enabled." if params.get("call_summits", True) else ""}{f" A total of {meta.get('peaks_found', 'N')} peaks were identified." if meta.get("peaks_found") else ""}""")

        # Differential analysis section
        if diff_steps:
            step = diff_steps[-1]
            params = step["parameters"]
            meta = step.get("output_metadata", {})

            sections.append(f"""### Differential Analysis

Differential binding analysis was performed using {params.get("method", "DESeq2")} with {params.get("norm_method", "RLE")} normalization. Peaks with FDR < {params.get("fdr_threshold", 0.05)} and |log2FC| > {params.get("fc_threshold", 1.0)} were considered significant.{f" {meta.get('significant_peaks', 'N')} differentially bound regions were identified ({meta.get('upregulated', 'N')} gained, {meta.get('downregulated', 'N')} lost)." if meta.get("significant_peaks") else ""}""")

        # Methylation section
        if meth_steps:
            step = meth_steps[-1]
            params = step["parameters"]

            sections.append(f"""### DNA Methylation Analysis

Bisulfite sequencing reads were processed using {params.get("aligner", "Bismark")}. Methylation calling was performed with minimum coverage of {params.get("min_coverage", 5)}x. Differentially methylated regions (DMRs) were identified with minimum methylation difference of {params.get("min_diff", 20)}% and minimum {params.get("min_cpgs", 3)} CpGs per region.""")

        # Add software versions footer
        all_versions = {}
        for step in history:
            all_versions.update(step.get("tool_versions", {}))

        if all_versions:
            version_list = ", ".join([f"{k} v{v}" for k, v in all_versions.items()])
            sections.append(f"""### Software Versions

Analysis was performed using EpiNexus v0.2.0 with the following tools: {version_list}.""")

        return "\n\n".join(sections)

    # -------------------------------------------------------------------------
    # Shell Script Generation
    # -------------------------------------------------------------------------

    @classmethod
    def generate_shell_script(cls) -> str:
        """Generate a shell script to reproduce the workflow."""
        history = cls.get_workflow_history()
        if not history:
            return "#!/bin/bash\n# No analysis steps recorded.\n"

        lines = [
            "#!/bin/bash",
            "# EpiNexus Workflow Reproduction Script",
            f"# Generated: {datetime.now().isoformat()}",
            f"# Genome: {st.session_state.get('selected_genome', 'N/A')}",
            "",
            "set -e  # Exit on error",
            "",
        ]

        for i, step in enumerate(history, 1):
            if step["status"] != "completed":
                continue

            lines.append(f"# Step {i}: {step['step_type'].replace('_', ' ').title()}")
            lines.append(f"# Timestamp: {step['timestamp']}")

            if step.get("command_line"):
                lines.append(step["command_line"])
            else:
                lines.append(f"# Parameters: {json.dumps(step['parameters'], indent=2)}")
                lines.append("# (Command not recorded - run through EpiNexus UI)")

            lines.append("")

        lines.append("echo 'Workflow complete!'")
        return "\n".join(lines)

    # -------------------------------------------------------------------------
    # Utilities
    # -------------------------------------------------------------------------

    @staticmethod
    def _calculate_checksum(filepath: str, algorithm: str = "sha256") -> str:
        """Calculate file checksum."""
        hash_func = hashlib.new(algorithm)
        try:
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    hash_func.update(chunk)
            return f"{algorithm}:{hash_func.hexdigest()}"
        except (IOError, OSError):
            return ""

    @staticmethod
    def get_tool_version(tool_name: str) -> str:
        """Get version of a command-line tool."""
        version_commands = {
            "bowtie2": ["bowtie2", "--version"],
            "bwa": ["bwa"],
            "samtools": ["samtools", "--version"],
            "macs2": ["macs2", "--version"],
            "bismark": ["bismark", "--version"],
            "methyldackel": ["MethylDackel", "--version"],
        }

        cmd = version_commands.get(tool_name.lower())
        if not cmd:
            return ""

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            output = result.stdout + result.stderr
            # Extract version number (first line usually contains it)
            first_line = output.strip().split("\n")[0]
            # Try to extract version pattern
            import re

            match = re.search(r"(\d+\.\d+(?:\.\d+)?)", first_line)
            return match.group(1) if match else first_line[:50]
        except (subprocess.SubprocessError, FileNotFoundError):
            return ""
