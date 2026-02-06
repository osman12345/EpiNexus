"""
Auto-Save and State Management

Provides:
- Automatic session saving at intervals
- Crash recovery
- State versioning
- Undo/Redo functionality
"""

import json
import time
import hashlib
import threading
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, asdict
from collections import deque
import copy


@dataclass
class StateSnapshot:
    """A snapshot of application state."""
    id: str
    timestamp: str
    state: Dict[str, Any]
    description: str = ""
    checksum: str = ""

    def __post_init__(self):
        if not self.checksum:
            self.checksum = self._compute_checksum()

    def _compute_checksum(self) -> str:
        """Compute checksum of state for change detection."""
        state_str = json.dumps(self.state, sort_keys=True, default=str)
        return hashlib.md5(state_str.encode()).hexdigest()[:12]


class AutoSaveManager:
    """Manage automatic saving of session state."""

    def __init__(
        self,
        save_dir: str = "autosave",
        interval_seconds: int = 60,
        max_autosaves: int = 10
    ):
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)
        self.interval = interval_seconds
        self.max_autosaves = max_autosaves

        self._last_checksum: Optional[str] = None
        self._timer: Optional[threading.Timer] = None
        self._running = False

    def start(self, get_state_func):
        """Start auto-save timer."""
        self._running = True
        self._get_state = get_state_func
        self._schedule_save()

    def stop(self):
        """Stop auto-save timer."""
        self._running = False
        if self._timer:
            self._timer.cancel()

    def _schedule_save(self):
        """Schedule next auto-save."""
        if self._running:
            self._timer = threading.Timer(self.interval, self._do_autosave)
            self._timer.daemon = True
            self._timer.start()

    def _do_autosave(self):
        """Perform auto-save if state changed."""
        try:
            state = self._get_state()
            if state:
                snapshot = StateSnapshot(
                    id=datetime.now().strftime("%Y%m%d_%H%M%S"),
                    timestamp=datetime.now().isoformat(),
                    state=state,
                    description="Auto-save"
                )

                # Only save if state changed
                if snapshot.checksum != self._last_checksum:
                    self._save_snapshot(snapshot)
                    self._last_checksum = snapshot.checksum
                    self._cleanup_old_saves()

        except Exception as e:
            print(f"Auto-save error: {e}")

        # Schedule next save
        self._schedule_save()

    def _save_snapshot(self, snapshot: StateSnapshot):
        """Save snapshot to disk."""
        filename = f"autosave_{snapshot.id}.json"
        filepath = self.save_dir / filename

        with open(filepath, 'w') as f:
            json.dump(asdict(snapshot), f, indent=2, default=str)

    def _cleanup_old_saves(self):
        """Remove old auto-saves beyond max limit."""
        saves = sorted(self.save_dir.glob("autosave_*.json"))
        while len(saves) > self.max_autosaves:
            saves[0].unlink()
            saves = saves[1:]

    def get_autosaves(self) -> List[Dict[str, Any]]:
        """Get list of available auto-saves."""
        saves = []
        for f in sorted(self.save_dir.glob("autosave_*.json"), reverse=True):
            try:
                with open(f, 'r') as file:
                    data = json.load(file)
                    saves.append({
                        "id": data.get("id", ""),
                        "timestamp": data.get("timestamp", ""),
                        "description": data.get("description", ""),
                        "file": str(f)
                    })
            except Exception:
                continue
        return saves

    def load_autosave(self, save_id: str) -> Optional[Dict[str, Any]]:
        """Load an auto-saved state."""
        filepath = self.save_dir / f"autosave_{save_id}.json"
        if filepath.exists():
            with open(filepath, 'r') as f:
                data = json.load(f)
                return data.get("state", {})
        return None

    def force_save(self, state: Dict[str, Any], description: str = "Manual save"):
        """Force an immediate save."""
        snapshot = StateSnapshot(
            id=datetime.now().strftime("%Y%m%d_%H%M%S"),
            timestamp=datetime.now().isoformat(),
            state=state,
            description=description
        )
        self._save_snapshot(snapshot)
        self._last_checksum = snapshot.checksum


class UndoRedoManager:
    """Manage undo/redo functionality."""

    def __init__(self, max_history: int = 50):
        self.max_history = max_history
        self._undo_stack: deque = deque(maxlen=max_history)
        self._redo_stack: deque = deque(maxlen=max_history)
        self._current_state: Optional[Dict[str, Any]] = None

    def push_state(self, state: Dict[str, Any], action_name: str = ""):
        """Push a new state to the undo stack."""
        if self._current_state is not None:
            self._undo_stack.append({
                "state": copy.deepcopy(self._current_state),
                "action": action_name,
                "timestamp": datetime.now().isoformat()
            })

        self._current_state = copy.deepcopy(state)
        self._redo_stack.clear()  # Clear redo stack on new action

    def undo(self) -> Optional[Dict[str, Any]]:
        """Undo last action and return previous state."""
        if not self._undo_stack:
            return None

        # Push current state to redo stack
        if self._current_state is not None:
            self._redo_stack.append({
                "state": copy.deepcopy(self._current_state),
                "action": "Undo",
                "timestamp": datetime.now().isoformat()
            })

        # Pop from undo stack
        previous = self._undo_stack.pop()
        self._current_state = copy.deepcopy(previous["state"])

        return self._current_state

    def redo(self) -> Optional[Dict[str, Any]]:
        """Redo last undone action."""
        if not self._redo_stack:
            return None

        # Push current state to undo stack
        if self._current_state is not None:
            self._undo_stack.append({
                "state": copy.deepcopy(self._current_state),
                "action": "Redo",
                "timestamp": datetime.now().isoformat()
            })

        # Pop from redo stack
        next_state = self._redo_stack.pop()
        self._current_state = copy.deepcopy(next_state["state"])

        return self._current_state

    def can_undo(self) -> bool:
        """Check if undo is available."""
        return len(self._undo_stack) > 0

    def can_redo(self) -> bool:
        """Check if redo is available."""
        return len(self._redo_stack) > 0

    def get_undo_history(self) -> List[Dict[str, str]]:
        """Get list of undoable actions."""
        return [
            {"action": item["action"], "timestamp": item["timestamp"]}
            for item in reversed(self._undo_stack)
        ]

    def get_redo_history(self) -> List[Dict[str, str]]:
        """Get list of redoable actions."""
        return [
            {"action": item["action"], "timestamp": item["timestamp"]}
            for item in reversed(self._redo_stack)
        ]

    def clear_history(self):
        """Clear all undo/redo history."""
        self._undo_stack.clear()
        self._redo_stack.clear()


def create_streamlit_state_manager():
    """Create state manager integrated with Streamlit session state."""
    import streamlit as st

    # Initialize managers in session state
    if 'autosave_manager' not in st.session_state:
        st.session_state.autosave_manager = AutoSaveManager()

    if 'undo_manager' not in st.session_state:
        st.session_state.undo_manager = UndoRedoManager()

    return st.session_state.autosave_manager, st.session_state.undo_manager


def render_undo_redo_buttons():
    """Render undo/redo buttons in Streamlit."""
    import streamlit as st

    _, undo_manager = create_streamlit_state_manager()

    col1, col2, col3 = st.columns([1, 1, 4])

    with col1:
        undo_disabled = not undo_manager.can_undo()
        if st.button("↩️ Undo", disabled=undo_disabled, key="undo_btn"):
            state = undo_manager.undo()
            if state:
                # Restore state
                for key, value in state.items():
                    st.session_state[key] = value
                st.rerun()

    with col2:
        redo_disabled = not undo_manager.can_redo()
        if st.button("↪️ Redo", disabled=redo_disabled, key="redo_btn"):
            state = undo_manager.redo()
            if state:
                for key, value in state.items():
                    st.session_state[key] = value
                st.rerun()

    with col3:
        # Show undo stack size
        st.caption(f"History: {len(undo_manager._undo_stack)} actions")


def render_autosave_recovery():
    """Render auto-save recovery UI."""
    import streamlit as st

    autosave_manager, _ = create_streamlit_state_manager()
    saves = autosave_manager.get_autosaves()

    if saves:
        with st.expander("🔄 Recover from Auto-Save"):
            for save in saves[:5]:
                col1, col2 = st.columns([3, 1])

                with col1:
                    st.markdown(f"**{save['timestamp'][:16]}** - {save['description']}")

                with col2:
                    if st.button("Restore", key=f"restore_{save['id']}"):
                        state = autosave_manager.load_autosave(save['id'])
                        if state:
                            for key, value in state.items():
                                st.session_state[key] = value
                            st.success("State restored!")
                            st.rerun()
