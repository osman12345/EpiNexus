"""
Workflow Export & Reproduce Page

Export complete analysis workflows for reproducibility,
import workflows from collaborators, and generate methods sections.

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
import json
from pathlib import Path
from datetime import datetime
import sys

# Add paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import workflow manager
try:
    from components.workflow_manager import WorkflowManager, ValidationResult
    HAS_WORKFLOW_MANAGER = True
except ImportError:
    try:
        from frontend.components.workflow_manager import WorkflowManager, ValidationResult
        HAS_WORKFLOW_MANAGER = True
    except ImportError:
        HAS_WORKFLOW_MANAGER = False

# =============================================================================
# PAGE CONFIG
# =============================================================================

st.set_page_config(
    page_title="Workflow Export - EpiNexus",
    page_icon="📤",
    layout="wide"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_timestamp(ts: str) -> str:
    """Format ISO timestamp for display."""
    try:
        dt = datetime.fromisoformat(ts)
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, TypeError):
        return ts

def format_runtime(seconds: float) -> str:
    """Format runtime in human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def get_step_icon(step_type: str) -> str:
    """Get icon for step type."""
    icons = {
        'alignment': '🧬',
        'peak_calling': '📍',
        'differential': '📈',
        'annotation': '🏷️',
        'methylation_alignment': '🔬',
        'methylation_calling': '🧬',
        'dmr_analysis': '📊',
        'super_enhancer': '⭐',
        'tf_analysis': '🔬',
        'motif_analysis': '🔤',
    }
    return icons.get(step_type, '⚙️')

def get_status_badge(status: str) -> str:
    """Get status badge HTML."""
    colors = {
        'completed': ('green', '✅'),
        'failed': ('red', '❌'),
        'running': ('orange', '⏳'),
    }
    color, icon = colors.get(status, ('gray', '❓'))
    return f"{icon} {status.title()}"

# =============================================================================
# MAIN PAGE
# =============================================================================

def main():
    st.title("📤 Workflow Export & Reproduce")
    st.markdown("Export, import, and reproduce analysis workflows for reproducibility")

    if not HAS_WORKFLOW_MANAGER:
        st.error("WorkflowManager not available. Please check installation.")
        return

    # Initialize workflow manager
    WorkflowManager.init_session_state()

    # Tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📜 History", "📥 Export", "📤 Import", "🔄 Reproduce", "📝 Methods"
    ])

    # =========================================================================
    # TAB 1: WORKFLOW HISTORY
    # =========================================================================
    with tab1:
        render_workflow_history()

    # =========================================================================
    # TAB 2: EXPORT
    # =========================================================================
    with tab2:
        render_export_tab()

    # =========================================================================
    # TAB 3: IMPORT
    # =========================================================================
    with tab3:
        render_import_tab()

    # =========================================================================
    # TAB 4: REPRODUCE
    # =========================================================================
    with tab4:
        render_reproduce_tab()

    # =========================================================================
    # TAB 5: METHODS SECTION
    # =========================================================================
    with tab5:
        render_methods_tab()


def render_workflow_history():
    """Render workflow history timeline."""
    st.subheader("Workflow Timeline")

    history = WorkflowManager.get_workflow_history()

    if not history:
        st.info("No analysis steps recorded yet. Run some analyses to see them here.")

        st.markdown("""
        ### Getting Started

        Workflow steps are automatically recorded when you run analyses:

        1. **Alignment** - When you align FASTQ files
        2. **Peak Calling** - When you call peaks with MACS2/SEACR
        3. **Differential Analysis** - When you run differential binding
        4. **Methylation** - When you analyze DNA methylation

        Once recorded, you can export your workflow to share with collaborators
        or reproduce it later.
        """)
        return

    # Summary stats
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Steps", len(history))
    with col2:
        completed = sum(1 for s in history if s['status'] == 'completed')
        st.metric("Completed", completed)
    with col3:
        failed = sum(1 for s in history if s['status'] == 'failed')
        st.metric("Failed", failed, delta_color="inverse" if failed > 0 else "off")
    with col4:
        total_runtime = sum(s.get('runtime_seconds', 0) for s in history)
        st.metric("Total Runtime", format_runtime(total_runtime))

    st.markdown("---")

    # Timeline view
    for i, step in enumerate(history, 1):
        icon = get_step_icon(step['step_type'])
        status_badge = get_status_badge(step['status'])
        timestamp = format_timestamp(step['timestamp'])
        runtime = format_runtime(step.get('runtime_seconds', 0))

        with st.expander(
            f"{icon} **Step {i}: {step['step_type'].replace('_', ' ').title()}** — {status_badge} — {timestamp}",
            expanded=False
        ):
            col1, col2 = st.columns([2, 1])

            with col1:
                st.markdown("**Parameters:**")
                st.json(step['parameters'])

                if step.get('inputs'):
                    st.markdown("**Inputs:**")
                    for name, path in step['inputs'].items():
                        st.text(f"  {name}: {path}")

                if step.get('outputs'):
                    st.markdown("**Outputs:**")
                    for name, path in step['outputs'].items():
                        st.text(f"  {name}: {path}")

            with col2:
                st.markdown("**Metadata:**")
                st.text(f"Step ID: {step['step_id']}")
                st.text(f"Runtime: {runtime}")

                if step.get('tool_versions'):
                    st.markdown("**Tool Versions:**")
                    for tool, version in step['tool_versions'].items():
                        st.text(f"  {tool}: {version}")

                if step.get('output_metadata'):
                    st.markdown("**Results:**")
                    for key, value in step['output_metadata'].items():
                        st.text(f"  {key}: {value}")

            if step.get('command_line'):
                st.markdown("**Command:**")
                st.code(step['command_line'], language='bash')

            if step.get('error_message'):
                st.error(f"Error: {step['error_message']}")

    # Clear history button
    st.markdown("---")
    col1, col2 = st.columns([3, 1])
    with col2:
        if st.button("🗑️ Clear History", type="secondary"):
            WorkflowManager.clear_history()
            st.rerun()


def render_export_tab():
    """Render workflow export interface."""
    st.subheader("Export Workflow")

    history = WorkflowManager.get_workflow_history()

    if not history:
        st.warning("No workflow steps to export. Run some analyses first.")
        return

    st.markdown(f"**{len(history)} steps** ready for export")

    # Export options
    col1, col2 = st.columns(2)

    with col1:
        include_checksums = st.checkbox(
            "Include file checksums",
            value=True,
            help="Calculate SHA256 hashes for output files to verify data integrity"
        )

        include_commands = st.checkbox(
            "Include command lines",
            value=True,
            help="Include the exact commands executed for each step"
        )

    with col2:
        export_format = st.selectbox(
            "Export Format",
            ["JSON (.epinexus.json)", "Shell Script (.sh)"],
            help="JSON for importing into EpiNexus, Shell script for manual execution"
        )

    st.markdown("---")

    # Preview
    st.markdown("### Preview")

    if export_format.startswith("JSON"):
        with st.spinner("Generating export..."):
            export_data = WorkflowManager.export_workflow(
                include_checksums=include_checksums,
                include_commands=include_commands
            )
            export_json = json.dumps(export_data, indent=2, default=str)

        with st.expander("View JSON", expanded=False):
            st.code(export_json, language='json')

        # Download button
        project_name = st.session_state.get('project_name', 'workflow')
        filename = f"{project_name}_{datetime.now().strftime('%Y%m%d')}.epinexus.json"

        st.download_button(
            "📥 Download Workflow (JSON)",
            export_json,
            filename,
            "application/json",
            type="primary",
            use_container_width=True
        )

    else:  # Shell script
        script = WorkflowManager.generate_shell_script()

        with st.expander("View Script", expanded=False):
            st.code(script, language='bash')

        project_name = st.session_state.get('project_name', 'workflow')
        filename = f"{project_name}_{datetime.now().strftime('%Y%m%d')}.sh"

        st.download_button(
            "📥 Download Shell Script",
            script,
            filename,
            "text/x-shellscript",
            type="primary",
            use_container_width=True
        )


def render_import_tab():
    """Render workflow import interface."""
    st.subheader("Import Workflow")

    st.markdown("""
    Import a workflow file from a collaborator or previous session.
    The workflow will be validated and you can review the steps before reproducing.
    """)

    uploaded_file = st.file_uploader(
        "Select workflow file",
        type=['json', 'epinexus'],
        help="Upload a .epinexus.json workflow file"
    )

    if uploaded_file:
        try:
            content = uploaded_file.read().decode('utf-8')
            workflow_data, validation = WorkflowManager.load_workflow_from_file(content)

            if workflow_data is None:
                st.error("Failed to parse workflow file")
                for error in validation.errors:
                    st.error(error)
                return

            # Store imported workflow in session state
            st.session_state.imported_workflow = workflow_data
            st.session_state.imported_validation = validation

            # Show validation results
            st.markdown("### Validation Results")

            col1, col2, col3 = st.columns(3)
            with col1:
                if validation.schema_valid:
                    st.success("✅ Schema valid")
                else:
                    st.error("❌ Schema invalid")

            with col2:
                missing = len(validation.files_missing)
                if missing == 0:
                    st.success("✅ All files found")
                else:
                    st.warning(f"⚠️ {missing} file(s) missing")

            with col3:
                if validation.tool_versions_compatible:
                    st.success("✅ Tools compatible")
                else:
                    st.warning("⚠️ Version mismatch")

            # Show warnings
            if validation.warnings:
                st.markdown("### Warnings")
                for warning in validation.warnings:
                    st.warning(warning)

            # Show workflow summary
            st.markdown("### Workflow Summary")

            metadata = workflow_data.get('metadata', {})
            project = workflow_data.get('project', {})
            steps = workflow_data.get('workflow_steps', [])

            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Name:** {metadata.get('name', 'Unknown')}")
                st.markdown(f"**Exported:** {format_timestamp(metadata.get('exported_at', ''))}")
                st.markdown(f"**EpiNexus Version:** {metadata.get('epinexus_version', 'Unknown')}")

            with col2:
                st.markdown(f"**Assay Type:** {project.get('assay_type', 'Unknown')}")
                st.markdown(f"**Genome:** {project.get('genome', 'Unknown')}")
                st.markdown(f"**Samples:** {len(project.get('samples', []))}")

            st.markdown(f"**Total Steps:** {len(steps)}")

            # Show steps
            st.markdown("### Steps")
            for i, step in enumerate(steps, 1):
                icon = get_step_icon(step['step_type'])
                st.markdown(f"{i}. {icon} **{step['step_type'].replace('_', ' ').title()}**")

            # Missing files
            if validation.files_missing:
                st.markdown("### Missing Files")
                st.markdown("The following input files were not found:")
                for filepath in validation.files_missing:
                    st.text(f"  ❌ {filepath}")
                st.info("You can specify new paths in the Reproduce tab.")

        except Exception as e:
            st.error(f"Error reading file: {e}")


def render_reproduce_tab():
    """Render workflow reproduction interface."""
    st.subheader("Reproduce Workflow")

    if 'imported_workflow' not in st.session_state:
        st.info("Import a workflow file first in the Import tab.")
        return

    workflow = st.session_state.imported_workflow
    validation = st.session_state.get('imported_validation', ValidationResult(is_valid=True))

    steps = workflow.get('workflow_steps', [])

    if not steps:
        st.warning("No steps in imported workflow.")
        return

    st.markdown(f"**{len(steps)} steps** available for reproduction")

    # Step selection
    st.markdown("### Select Steps to Reproduce")

    selected_steps = []
    for i, step in enumerate(steps):
        icon = get_step_icon(step['step_type'])
        step_name = step['step_type'].replace('_', ' ').title()

        col1, col2, col3 = st.columns([0.5, 2, 1])
        with col1:
            selected = st.checkbox("", value=True, key=f"step_{i}")
            if selected:
                selected_steps.append(i)
        with col2:
            st.markdown(f"{icon} **{step_name}**")
            st.caption(f"Parameters: {len(step.get('parameters', {}))} | Inputs: {len(step.get('inputs', {}))}")
        with col3:
            if st.button("Edit", key=f"edit_{i}"):
                st.session_state[f"editing_step_{i}"] = True

        # Parameter editing
        if st.session_state.get(f"editing_step_{i}", False):
            with st.expander(f"Edit {step_name} Parameters", expanded=True):
                edited_params = {}
                for key, value in step.get('parameters', {}).items():
                    if isinstance(value, bool):
                        edited_params[key] = st.checkbox(key, value=value, key=f"param_{i}_{key}")
                    elif isinstance(value, int):
                        edited_params[key] = st.number_input(key, value=value, key=f"param_{i}_{key}")
                    elif isinstance(value, float):
                        edited_params[key] = st.number_input(key, value=value, key=f"param_{i}_{key}")
                    else:
                        edited_params[key] = st.text_input(key, value=str(value), key=f"param_{i}_{key}")

                if st.button("Save Changes", key=f"save_{i}"):
                    st.session_state.imported_workflow['workflow_steps'][i]['parameters'] = edited_params
                    st.session_state[f"editing_step_{i}"] = False
                    st.rerun()

    st.markdown("---")

    # File path remapping for missing files
    if validation.files_missing:
        st.markdown("### Remap Missing Files")
        st.markdown("Specify new paths for missing input files:")

        remapped_files = {}
        for filepath in validation.files_missing:
            new_path = st.text_input(
                f"New path for: {Path(filepath).name}",
                value="",
                key=f"remap_{filepath}"
            )
            if new_path:
                remapped_files[filepath] = new_path

        st.session_state.remapped_files = remapped_files

    st.markdown("---")

    # Run options
    col1, col2 = st.columns(2)

    with col1:
        if st.button("▶️ Run Selected Steps", type="primary", use_container_width=True):
            st.info("Reproduction would start here with selected steps...")
            st.warning("Note: Actual execution requires the original tools and data files.")

    with col2:
        if st.button("📄 Export as Shell Script", use_container_width=True):
            # Generate script for selected steps only
            script_lines = [
                "#!/bin/bash",
                f"# Reproduced workflow from {workflow.get('metadata', {}).get('name', 'Unknown')}",
                f"# Generated: {datetime.now().isoformat()}",
                "",
                "set -e",
                ""
            ]

            for i in selected_steps:
                step = steps[i]
                script_lines.append(f"# Step {i+1}: {step['step_type']}")
                if step.get('command_line'):
                    script_lines.append(step['command_line'])
                else:
                    script_lines.append(f"# Parameters: {json.dumps(step['parameters'])}")
                script_lines.append("")

            script = "\n".join(script_lines)

            st.download_button(
                "Download Script",
                script,
                "reproduce_workflow.sh",
                "text/x-shellscript"
            )


def render_methods_tab():
    """Render methods section generation."""
    st.subheader("Methods Section")

    st.markdown("""
    Generate a publication-ready methods section from your workflow history.
    This describes the analysis steps in a format suitable for scientific publications.
    """)

    history = WorkflowManager.get_workflow_history()

    if not history:
        st.warning("No workflow steps recorded. Run some analyses to generate methods text.")
        return

    # Generate methods
    methods_text = WorkflowManager.generate_methods_text()

    st.markdown("### Generated Methods")
    st.markdown(methods_text)

    st.markdown("---")

    # Download options
    col1, col2 = st.columns(2)

    with col1:
        st.download_button(
            "📥 Download as Markdown",
            methods_text,
            "methods_section.md",
            "text/markdown",
            use_container_width=True
        )

    with col2:
        # Plain text version
        plain_text = methods_text.replace("### ", "").replace("**", "")
        st.download_button(
            "📥 Download as Plain Text",
            plain_text,
            "methods_section.txt",
            "text/plain",
            use_container_width=True
        )

    # Tips
    st.markdown("---")
    st.markdown("""
    ### Tips for Methods Sections

    - Always review and edit the generated text before publication
    - Add specific details about your samples and experimental design
    - Include accession numbers for reference genomes and databases
    - Cite the original tools (Bowtie2, MACS2, DESeq2, etc.)
    - Mention the EpiNexus version used for analysis
    """)


if __name__ == "__main__":
    main()
