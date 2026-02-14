"""
Multi-Mark Integration Module

Integrates differential peaks from multiple histone marks to identify:
- Antagonistic chromatin changes (H3K27ac vs H3K27me3)
- Chromatin state transitions
- Coordinated regulatory changes

Supports two-mark and three-mark integration with optional RNA-seq data.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class IntegrationConfig:
    """Configuration for multi-mark integration."""

    integration_type: str = "two_mark"  # two_mark or three_mark

    # Input files (paths to DiffBind results)
    h3k27ac_file: Optional[str] = None
    h3k27me3_file: Optional[str] = None
    h3k4me1_file: Optional[str] = None

    # Optional RNA-seq
    rnaseq_file: Optional[str] = None

    # Thresholds
    peak_fdr: float = 0.1
    de_fdr: float = 0.05
    de_lfc: float = 0.3

    # Output
    output_dir: Optional[str] = None


@dataclass
class IntegrationResults:
    """Results from integration analysis."""

    integration_type: str
    total_genes: int

    # Chromatin states
    chromatin_states: Dict[str, int]

    # High-confidence targets
    high_confidence_count: int

    # Output files
    integrated_file: str
    high_confidence_file: Optional[str]
    plots_file: Optional[str]

    # Summary
    summary: Dict[str, Any]


class MarkIntegration:
    """Multi-mark integration analyzer."""

    # Chromatin state definitions
    CHROMATIN_STATES = {
        "Active_enhancer": {"ac": True, "me1": True, "me3": False},
        "Poised_enhancer": {"ac": False, "me1": True, "me3": False},
        "Active_promoter": {"ac": True, "me1": False, "me3": False},
        "Bivalent": {"ac": False, "me1": True, "me3": True},
        "Repressed": {"ac": False, "me1": False, "me3": True},
        "Complex": {"ac": True, "me1": True, "me3": True},
    }

    def __init__(self):
        """Initialize integration analyzer."""
        pass

    def load_diffbind_results(self, file_path: str, mark_name: str) -> pd.DataFrame:
        """
        Load and preprocess DiffBind results.

        Args:
            file_path: Path to DiffBind results CSV
            mark_name: Histone mark name for column prefixes

        Returns:
            Preprocessed DataFrame
        """
        df = pd.read_csv(file_path)

        # Standardize column names
        col_mapping = {"gene_name": "gene_symbol", "geneId": "ensembl_id", "Fold": "fold_change", "FDR": "fdr"}

        df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})

        # Add mark-specific prefixes
        cols_to_prefix = ["fold_change", "fdr", "significant", "direction"]
        for col in cols_to_prefix:
            if col in df.columns:
                df = df.rename(columns={col: f"{mark_name}_{col}"})

        logger.info(f"Loaded {len(df)} peaks for {mark_name}")
        return df

    def load_rnaseq_results(self, file_path: str) -> pd.DataFrame:
        """
        Load RNA-seq differential expression results.

        Args:
            file_path: Path to DE results CSV

        Returns:
            Preprocessed DataFrame
        """
        df = pd.read_csv(file_path)

        # Standardize column names (handle different DE tools)
        col_mapping = {
            "log2FoldChange": "log2FC",
            "logFC": "log2FC",
            "padj": "expression_FDR",
            "adj.P.Val": "expression_FDR",
            "FDR": "expression_FDR",
            "gene": "gene_symbol",
            "gene_name": "gene_symbol",
            "ENSEMBL": "ensembl_id",
            "ensembl_gene_id": "ensembl_id",
        }

        df = df.rename(columns={k: v for k, v in col_mapping.items() if k in df.columns})

        logger.info(f"Loaded {len(df)} genes from RNA-seq")
        return df

    def run_two_mark_integration(self, config: IntegrationConfig) -> IntegrationResults:
        """
        Run two-mark integration (H3K27ac vs H3K27me3).

        Identifies antagonistic chromatin changes:
        - Gained H3K27ac + Lost H3K27me3 = Strong activation
        - Lost H3K27ac + Gained H3K27me3 = Strong repression

        Args:
            config: Integration configuration

        Returns:
            IntegrationResults object
        """
        if not config.h3k27ac_file or not config.h3k27me3_file:
            raise ValueError("Two-mark integration requires H3K27ac and H3K27me3 files")

        # Load data
        ac_df = self.load_diffbind_results(config.h3k27ac_file, "ac")
        me3_df = self.load_diffbind_results(config.h3k27me3_file, "me3")

        # Merge on gene
        integrated = self._merge_marks(ac_df, me3_df, "ac", "me3")

        # Add RNA-seq if available
        if config.rnaseq_file:
            rnaseq_df = self.load_rnaseq_results(config.rnaseq_file)
            integrated = integrated.merge(
                rnaseq_df[["ensembl_id", "log2FC", "expression_FDR"]], on="ensembl_id", how="left"
            )

        # Classify changes
        integrated = self._classify_two_mark_changes(integrated, config)

        # Identify high-confidence targets
        high_conf = self._identify_high_confidence(integrated, config)

        # Save results
        output_dir = Path(config.output_dir) if config.output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        integrated_file = output_dir / "two_mark_integrated.csv"
        integrated.to_csv(integrated_file, index=False)

        high_conf_file = None
        if len(high_conf) > 0:
            high_conf_file = output_dir / "two_mark_high_confidence.csv"
            high_conf.to_csv(high_conf_file, index=False)

        # Generate summary
        summary = self._generate_two_mark_summary(integrated, high_conf)

        return IntegrationResults(
            integration_type="two_mark",
            total_genes=len(integrated),
            chromatin_states=summary.get("categories", {}),
            high_confidence_count=len(high_conf),
            integrated_file=str(integrated_file),
            high_confidence_file=str(high_conf_file) if high_conf_file else None,
            plots_file=None,  # TODO: Generate plots
            summary=summary,
        )

    def run_three_mark_integration(self, config: IntegrationConfig) -> IntegrationResults:
        """
        Run three-mark integration (H3K27ac + H3K27me3 + H3K4me1).

        Classifies chromatin states:
        - Active enhancer: H3K27ac + H3K4me1
        - Poised enhancer: H3K4me1 only
        - Active promoter: H3K27ac only
        - Bivalent: H3K27me3 + H3K4me1
        - Repressed: H3K27me3 only

        Args:
            config: Integration configuration

        Returns:
            IntegrationResults object
        """
        if not all([config.h3k27ac_file, config.h3k27me3_file, config.h3k4me1_file]):
            raise ValueError("Three-mark integration requires H3K27ac, H3K27me3, and H3K4me1 files")

        # Load data
        ac_df = self.load_diffbind_results(config.h3k27ac_file, "ac")
        me3_df = self.load_diffbind_results(config.h3k27me3_file, "me3")
        me1_df = self.load_diffbind_results(config.h3k4me1_file, "me1")

        # Merge all marks
        integrated = self._merge_three_marks(ac_df, me3_df, me1_df)

        # Add RNA-seq if available
        if config.rnaseq_file:
            rnaseq_df = self.load_rnaseq_results(config.rnaseq_file)
            integrated = integrated.merge(
                rnaseq_df[["ensembl_id", "log2FC", "expression_FDR"]], on="ensembl_id", how="left"
            )

        # Classify chromatin states
        integrated = self._classify_chromatin_states(integrated, config)

        # Identify chromatin transitions
        integrated = self._classify_transitions(integrated)

        # High-confidence targets
        high_conf = integrated[
            (integrated["n_marks"] >= 2)
            & (
                integrated["chromatin_transition"].isin(
                    ["Strong_activation", "Strong_repression", "Activation", "Repression"]
                )
            )
        ].copy()

        # Save results
        output_dir = Path(config.output_dir) if config.output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        integrated_file = output_dir / "three_mark_integrated.csv"
        integrated.to_csv(integrated_file, index=False)

        high_conf_file = None
        if len(high_conf) > 0:
            high_conf_file = output_dir / "three_mark_high_confidence.csv"
            high_conf.to_csv(high_conf_file, index=False)

        # Summary
        state_counts = integrated["chromatin_state"].value_counts().to_dict()
        summary = {
            "total_genes": len(integrated),
            "chromatin_states": state_counts,
            "high_confidence_targets": len(high_conf),
            "genes_with_all_three": int((integrated["n_marks"] == 3).sum()),
        }

        return IntegrationResults(
            integration_type="three_mark",
            total_genes=len(integrated),
            chromatin_states=state_counts,
            high_confidence_count=len(high_conf),
            integrated_file=str(integrated_file),
            high_confidence_file=str(high_conf_file) if high_conf_file else None,
            plots_file=None,
            summary=summary,
        )

    def _merge_marks(self, df1: pd.DataFrame, df2: pd.DataFrame, mark1: str, mark2: str) -> pd.DataFrame:
        """Merge two mark DataFrames on gene."""
        # Get gene columns
        gene_cols = ["ensembl_id", "gene_symbol"]

        # Select relevant columns
        df1_cols = [c for c in df1.columns if c.startswith(mark1) or c in gene_cols]
        df2_cols = [c for c in df2.columns if c.startswith(mark2) or c in gene_cols]

        df1_subset = df1[df1_cols].drop_duplicates(subset=["ensembl_id"])
        df2_subset = df2[df2_cols].drop_duplicates(subset=["ensembl_id"])

        merged = df1_subset.merge(df2_subset, on="ensembl_id", how="outer")

        # Coalesce gene_symbol
        if "gene_symbol_x" in merged.columns:
            merged["gene_symbol"] = merged["gene_symbol_x"].fillna(merged["gene_symbol_y"])
            merged = merged.drop(columns=["gene_symbol_x", "gene_symbol_y"])

        return merged

    def _merge_three_marks(self, ac_df: pd.DataFrame, me3_df: pd.DataFrame, me1_df: pd.DataFrame) -> pd.DataFrame:
        """Merge three mark DataFrames."""
        # First merge ac and me3
        merged = self._merge_marks(ac_df, me3_df, "ac", "me3")

        # Then add me1
        gene_cols = ["ensembl_id", "gene_symbol"]
        me1_cols = [c for c in me1_df.columns if c.startswith("me1") or c in gene_cols]
        me1_subset = me1_df[me1_cols].drop_duplicates(subset=["ensembl_id"])

        merged = merged.merge(me1_subset, on="ensembl_id", how="outer")

        # Coalesce gene_symbol
        if "gene_symbol_x" in merged.columns:
            merged["gene_symbol"] = merged["gene_symbol_x"].fillna(merged["gene_symbol_y"])
            merged = merged.drop(columns=["gene_symbol_x", "gene_symbol_y"])

        return merged

    def _classify_two_mark_changes(self, df: pd.DataFrame, config: IntegrationConfig) -> pd.DataFrame:
        """Classify two-mark antagonistic changes."""
        df = df.copy()

        # Determine significance
        df["ac_sig"] = df.get("ac_fdr", 1) < config.peak_fdr
        df["me3_sig"] = df.get("me3_fdr", 1) < config.peak_fdr

        # Determine direction
        df["ac_gained"] = (df.get("ac_fold_change", 0) > 0) & df["ac_sig"]
        df["ac_lost"] = (df.get("ac_fold_change", 0) < 0) & df["ac_sig"]
        df["me3_gained"] = (df.get("me3_fold_change", 0) > 0) & df["me3_sig"]
        df["me3_lost"] = (df.get("me3_fold_change", 0) < 0) & df["me3_sig"]

        # Classify integration category
        def classify(row):
            if row["ac_gained"] and row["me3_lost"]:
                return "Strong_activation"
            elif row["ac_lost"] and row["me3_gained"]:
                return "Strong_repression"
            elif row["ac_gained"] and not row["me3_sig"]:
                return "Activation_ac_only"
            elif row["me3_lost"] and not row["ac_sig"]:
                return "Activation_me3_only"
            elif row["ac_lost"] and not row["me3_sig"]:
                return "Repression_ac_only"
            elif row["me3_gained"] and not row["ac_sig"]:
                return "Repression_me3_only"
            elif row["ac_gained"] and row["me3_gained"]:
                return "Discordant_both_gained"
            elif row["ac_lost"] and row["me3_lost"]:
                return "Discordant_both_lost"
            else:
                return "No_significant_change"

        df["integration_category"] = df.apply(classify, axis=1)

        return df

    def _classify_chromatin_states(self, df: pd.DataFrame, config: IntegrationConfig) -> pd.DataFrame:
        """Classify chromatin states based on mark presence."""
        df = df.copy()

        # Mark presence
        df["has_ac"] = ~df.get("ac_fdr", pd.Series([np.nan])).isna()
        df["has_me3"] = ~df.get("me3_fdr", pd.Series([np.nan])).isna()
        df["has_me1"] = ~df.get("me1_fdr", pd.Series([np.nan])).isna()

        df["n_marks"] = df["has_ac"].astype(int) + df["has_me3"].astype(int) + df["has_me1"].astype(int)

        # Classify state
        def get_state(row):
            ac, me3, me1 = row["has_ac"], row["has_me3"], row["has_me1"]

            if ac and me1 and me3:
                return "Complex"
            elif ac and me1 and not me3:
                return "Active_enhancer"
            elif not ac and me1 and me3:
                return "Bivalent"
            elif not ac and me1 and not me3:
                return "Poised_enhancer"
            elif not ac and not me1 and me3:
                return "Repressed"
            elif ac and not me1 and not me3:
                return "Active_promoter"
            else:
                return "Other"

        df["chromatin_state"] = df.apply(get_state, axis=1)

        return df

    def _classify_transitions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Classify chromatin state transitions."""
        df = df.copy()

        # Get directional changes
        df["ac_gained"] = df.get("ac_fold_change", 0) > 0
        df["ac_lost"] = df.get("ac_fold_change", 0) < 0
        df["me3_gained"] = df.get("me3_fold_change", 0) > 0
        df["me3_lost"] = df.get("me3_fold_change", 0) < 0
        df["me1_gained"] = df.get("me1_fold_change", 0) > 0
        df["me1_lost"] = df.get("me1_fold_change", 0) < 0

        def get_transition(row):
            # Strong activation: lose me3, gain ac and me1
            if row.get("me3_lost", False) and row.get("ac_gained", False) and row.get("me1_gained", False):
                return "Strong_activation"
            # Strong repression: gain me3, lose ac and me1
            elif row.get("me3_gained", False) and row.get("ac_lost", False) and row.get("me1_lost", False):
                return "Strong_repression"
            # Activation
            elif row.get("me3_lost", False) or (row.get("ac_gained", False) and row.get("me1_gained", False)):
                if not (row.get("me3_gained", False) or (row.get("ac_lost", False) and row.get("me1_lost", False))):
                    return "Activation"
            # Repression
            elif row.get("me3_gained", False) or (row.get("ac_lost", False) and row.get("me1_lost", False)):
                if not (row.get("me3_lost", False) or (row.get("ac_gained", False) and row.get("me1_gained", False))):
                    return "Repression"

            return "Mixed"

        df["chromatin_transition"] = df.apply(get_transition, axis=1)

        return df

    def _identify_high_confidence(self, df: pd.DataFrame, config: IntegrationConfig) -> pd.DataFrame:
        """Identify high-confidence coordinated targets."""
        # Filter for strong changes
        high_conf = df[df["integration_category"].isin(["Strong_activation", "Strong_repression"])].copy()

        # If RNA-seq available, also require expression change
        if "expression_FDR" in df.columns:
            high_conf = high_conf[
                (high_conf["expression_FDR"] < config.de_fdr) & (abs(high_conf["log2FC"]) > config.de_lfc)
            ]

        return high_conf.sort_values(by=["ac_fdr", "me3_fdr"], ascending=True)

    def _generate_two_mark_summary(self, integrated: pd.DataFrame, high_conf: pd.DataFrame) -> Dict[str, Any]:
        """Generate summary statistics for two-mark integration."""
        categories = integrated["integration_category"].value_counts().to_dict()

        return {
            "total_genes": len(integrated),
            "genes_with_both_marks": int(integrated["ac_sig"].fillna(False) | integrated["me3_sig"].fillna(False)).sum()
            if "ac_sig" in integrated.columns
            else 0,
            "categories": categories,
            "high_confidence_targets": len(high_conf),
            "strong_activation": categories.get("Strong_activation", 0),
            "strong_repression": categories.get("Strong_repression", 0),
        }
