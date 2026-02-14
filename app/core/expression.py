"""
Gene Expression Integration Module

Integrates RNA-seq differential expression data with:
- Histone mark changes
- TF binding data
- Chromatin state analysis

Supports multiple DE tools:
- DESeq2 output
- edgeR output
- limma output
- Custom formats
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ExpressionData:
    """Container for gene expression data."""

    # Raw data
    counts: Optional[pd.DataFrame] = None  # genes x samples
    normalized: Optional[pd.DataFrame] = None  # TPM/FPKM/CPM

    # Differential expression
    de_results: Optional[pd.DataFrame] = None

    # Metadata
    sample_info: Optional[pd.DataFrame] = None
    gene_info: Optional[pd.DataFrame] = None

    # Analysis info
    comparison_name: str = ""
    de_tool: str = "unknown"


@dataclass
class IntegratedGene:
    """A gene with integrated epigenetic and expression data."""

    gene_id: str
    gene_symbol: str

    # Expression
    log2fc: float
    expression_fdr: float
    expression_direction: str  # up, down, unchanged
    base_mean: Optional[float] = None

    # Histone marks
    histone_marks: Dict[str, Dict] = field(default_factory=dict)
    # e.g., {"H3K27ac": {"direction": "gained", "fdr": 0.01, "log2fc": 1.5}}

    # TF binding
    tf_binding: List[str] = field(default_factory=list)  # TF names
    tf_motifs: List[str] = field(default_factory=list)  # Enriched motifs

    # Integration scores
    chromatin_state: Optional[str] = None
    concordance_score: float = 0.0  # How well marks agree with expression
    regulatory_category: str = "unknown"


@dataclass
class IntegrationResult:
    """Results from expression + epigenetic integration."""

    total_genes: int
    de_genes: int
    integrated_genes: int

    # Categories
    concordant_activation: int = 0  # Up expression + active marks
    concordant_repression: int = 0  # Down expression + repressive marks
    discordant: int = 0

    # Data
    genes: List[IntegratedGene] = field(default_factory=list)
    summary_df: Optional[pd.DataFrame] = None


class ExpressionLoader:
    """
    Load and parse gene expression data from various formats.

    Supports:
    - DESeq2 results
    - edgeR results
    - limma-voom results
    - Generic CSV with log2FC and FDR columns
    """

    # Column name mappings for different tools
    COLUMN_MAPPINGS = {
        "deseq2": {
            "gene_id": ["gene_id", "gene", "ensembl_gene_id", "ENSEMBL"],
            "gene_symbol": ["gene_name", "symbol", "gene_symbol", "SYMBOL"],
            "log2fc": ["log2FoldChange", "log2FC", "logFC"],
            "pvalue": ["pvalue", "PValue", "P.Value"],
            "fdr": ["padj", "FDR", "adj.P.Val", "q_value"],
            "basemean": ["baseMean", "AveExpr", "logCPM"],
        },
        "edger": {
            "gene_id": ["gene_id", "genes", "GeneID"],
            "gene_symbol": ["gene_name", "Symbol"],
            "log2fc": ["logFC", "log2FoldChange"],
            "pvalue": ["PValue", "pvalue"],
            "fdr": ["FDR", "padj"],
            "basemean": ["logCPM", "AveExpr"],
        },
        "limma": {
            "gene_id": ["gene_id", "ID", "ProbeID"],
            "gene_symbol": ["gene_name", "Symbol", "Gene.symbol"],
            "log2fc": ["logFC", "log2FoldChange"],
            "pvalue": ["P.Value", "pvalue"],
            "fdr": ["adj.P.Val", "FDR"],
            "basemean": ["AveExpr"],
        },
    }

    def __init__(self):
        """Initialize expression loader."""
        pass

    def load(
        self,
        filepath: str,
        tool: str = "auto",
        gene_id_col: Optional[str] = None,
        log2fc_col: Optional[str] = None,
        fdr_col: Optional[str] = None,
    ) -> ExpressionData:
        """
        Load differential expression results.

        Args:
            filepath: Path to DE results file
            tool: DE tool used (deseq2, edger, limma, or auto)
            gene_id_col: Override gene ID column name
            log2fc_col: Override log2FC column name
            fdr_col: Override FDR column name

        Returns:
            ExpressionData object
        """
        # Load file
        df = self._read_file(filepath)

        # Auto-detect tool if needed
        if tool == "auto":
            tool = self._detect_tool(df)

        # Standardize columns
        df = self._standardize_columns(df, tool, gene_id_col, log2fc_col, fdr_col)

        # Create ExpressionData
        expr_data = ExpressionData(de_results=df, de_tool=tool, comparison_name=Path(filepath).stem)

        logger.info(f"Loaded {len(df)} genes from {filepath} (detected: {tool})")
        return expr_data

    def _read_file(self, filepath: str) -> pd.DataFrame:
        """Read expression file (CSV, TSV, or Excel)."""
        filepath = str(filepath)

        if filepath.endswith(".xlsx") or filepath.endswith(".xls"):
            return pd.read_excel(filepath)
        elif filepath.endswith(".tsv") or filepath.endswith(".txt"):
            return pd.read_csv(filepath, sep="\t")
        else:
            # Try to detect separator
            with open(filepath) as f:
                first_line = f.readline()
            sep = "\t" if "\t" in first_line else ","
            return pd.read_csv(filepath, sep=sep)

    def _detect_tool(self, df: pd.DataFrame) -> str:
        """Auto-detect which DE tool generated the results."""
        cols = set(df.columns.str.lower())

        if "log2foldchange" in cols and "padj" in cols:
            return "deseq2"
        elif "logfc" in cols and "fdr" in cols:
            return "edger"
        elif "logfc" in cols and "adj.p.val" in cols.union(set(df.columns)):
            return "limma"
        else:
            return "generic"

    def _standardize_columns(
        self,
        df: pd.DataFrame,
        tool: str,
        gene_id_col: Optional[str] = None,
        log2fc_col: Optional[str] = None,
        fdr_col: Optional[str] = None,
    ) -> pd.DataFrame:
        """Standardize column names to common format."""
        df = df.copy()

        # Get column mappings
        mappings = self.COLUMN_MAPPINGS.get(tool, self.COLUMN_MAPPINGS["deseq2"])

        # Helper function to find matching column
        def find_col(target_names: List[str], override: Optional[str] = None) -> Optional[str]:
            if override and override in df.columns:
                return override
            for name in target_names:
                if name in df.columns:
                    return name
                # Case-insensitive search
                for col in df.columns:
                    if col.lower() == name.lower():
                        return col
            return None

        # Map columns
        rename_map = {}

        # Gene ID
        col = find_col(mappings["gene_id"], gene_id_col)
        if col:
            rename_map[col] = "gene_id"
        elif df.index.name:
            df["gene_id"] = df.index

        # Gene symbol
        col = find_col(mappings["gene_symbol"])
        if col:
            rename_map[col] = "gene_symbol"

        # Log2FC
        col = find_col(mappings["log2fc"], log2fc_col)
        if col:
            rename_map[col] = "log2FC"

        # P-value
        col = find_col(mappings["pvalue"])
        if col:
            rename_map[col] = "pvalue"

        # FDR
        col = find_col(mappings["fdr"], fdr_col)
        if col:
            rename_map[col] = "FDR"

        # Base mean
        col = find_col(mappings["basemean"])
        if col:
            rename_map[col] = "baseMean"

        # Apply renaming
        df = df.rename(columns=rename_map)

        # Add derived columns
        if "log2FC" in df.columns and "FDR" in df.columns:
            df["direction"] = np.where(df["log2FC"] > 0, "up", "down")
            df["significant"] = df["FDR"] < 0.05

        return df

    def load_count_matrix(
        self, filepath: str, sample_info: Optional[str] = None
    ) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        """
        Load raw count matrix.

        Args:
            filepath: Path to count matrix (genes x samples)
            sample_info: Optional path to sample metadata

        Returns:
            Tuple of (counts_df, sample_info_df)
        """
        counts = self._read_file(filepath)

        # Set gene ID as index if present
        if "gene_id" in counts.columns:
            counts = counts.set_index("gene_id")
        elif counts.columns[0].lower() in ["gene", "geneid", "gene_id"]:
            counts = counts.set_index(counts.columns[0])

        sample_df = None
        if sample_info:
            sample_df = self._read_file(sample_info)

        return counts, sample_df


class ExpressionAnalyzer:
    """
    Analyze gene expression data and integrate with epigenetics.
    """

    def __init__(self, fdr_threshold: float = 0.05, log2fc_threshold: float = 0.5):
        """
        Initialize analyzer.

        Args:
            fdr_threshold: FDR cutoff for significance
            log2fc_threshold: Minimum absolute log2FC
        """
        self.fdr_threshold = fdr_threshold
        self.log2fc_threshold = log2fc_threshold

    def get_de_genes(self, expr_data: ExpressionData, direction: str = "both") -> pd.DataFrame:
        """
        Get differentially expressed genes.

        Args:
            expr_data: Expression data
            direction: 'up', 'down', or 'both'

        Returns:
            DataFrame of DE genes
        """
        df = expr_data.de_results.copy()

        # Filter by significance
        mask = (df["FDR"] < self.fdr_threshold) & (np.abs(df["log2FC"]) > self.log2fc_threshold)

        if direction == "up":
            mask &= df["log2FC"] > 0
        elif direction == "down":
            mask &= df["log2FC"] < 0

        return df[mask].sort_values("FDR")

    def integrate_with_peaks(
        self, expr_data: ExpressionData, peak_data: pd.DataFrame, gene_col: str = "gene_symbol"
    ) -> pd.DataFrame:
        """
        Integrate expression with peak/histone data.

        Args:
            expr_data: Expression data
            peak_data: Peak annotation with gene assignments
            gene_col: Column containing gene names

        Returns:
            Merged DataFrame
        """
        expr_df = expr_data.de_results

        # Get unique genes from peaks
        if gene_col in peak_data.columns:
            peak_genes = (
                peak_data.groupby(gene_col)
                .agg(
                    {
                        "peak_id": "count",
                        "log2FC": "mean" if "log2FC" in peak_data.columns else "first",
                        "FDR": "min" if "FDR" in peak_data.columns else "first",
                    }
                )
                .reset_index()
            )

            peak_genes.columns = [gene_col, "n_peaks", "peak_log2FC", "peak_FDR"]
        else:
            peak_genes = peak_data

        # Merge
        merged = expr_df.merge(
            peak_genes, left_on="gene_symbol", right_on=gene_col, how="outer", suffixes=("_expr", "_peak")
        )

        return merged


class ExpressionHistoneIntegrator:
    """
    Integrate gene expression with histone modification data.

    Identifies:
    - Concordant changes (expression matches chromatin)
    - Discordant changes (potential regulatory complexity)
    - Key regulatory targets
    """

    # Expected relationships between marks and expression
    MARK_EXPRESSION_LOGIC = {
        "H3K27ac": {"gained": "up", "lost": "down"},  # Active mark
        "H3K4me3": {"gained": "up", "lost": "down"},  # Active promoter
        "H3K4me1": {"gained": "up", "lost": "down"},  # Enhancer
        "H3K27me3": {"gained": "down", "lost": "up"},  # Repressive
        "H3K9me3": {"gained": "down", "lost": "up"},  # Heterochromatin
    }

    def __init__(self, expr_fdr: float = 0.05, expr_lfc: float = 0.5, peak_fdr: float = 0.1):
        """
        Initialize integrator.

        Args:
            expr_fdr: Expression FDR threshold
            expr_lfc: Expression log2FC threshold
            peak_fdr: Peak FDR threshold
        """
        self.expr_fdr = expr_fdr
        self.expr_lfc = expr_lfc
        self.peak_fdr = peak_fdr

    def integrate(
        self, expression: ExpressionData, histone_data: Dict[str, pd.DataFrame], tf_data: Optional[pd.DataFrame] = None
    ) -> IntegrationResult:
        """
        Integrate expression with histone marks and TF binding.

        Args:
            expression: Gene expression data
            histone_data: Dict of mark name -> annotated peaks DataFrame
            tf_data: Optional TF binding data

        Returns:
            IntegrationResult
        """
        expr_df = expression.de_results
        integrated_genes = []

        # Get all genes
        all_genes = set(expr_df["gene_symbol"].dropna().unique())

        for mark_name, peak_df in histone_data.items():
            if "gene_symbol" in peak_df.columns:
                all_genes.update(peak_df["gene_symbol"].dropna().unique())

        logger.info(f"Integrating {len(all_genes)} genes across {len(histone_data)} marks")

        # Process each gene
        for gene in all_genes:
            integrated = self._integrate_single_gene(gene, expr_df, histone_data, tf_data)
            if integrated:
                integrated_genes.append(integrated)

        # Calculate summary statistics
        result = self._summarize_integration(integrated_genes, expr_df)

        return result

    def _integrate_single_gene(
        self, gene: str, expr_df: pd.DataFrame, histone_data: Dict[str, pd.DataFrame], tf_data: Optional[pd.DataFrame]
    ) -> Optional[IntegratedGene]:
        """Integrate data for a single gene."""
        # Get expression data
        gene_expr = expr_df[expr_df["gene_symbol"] == gene]

        if len(gene_expr) == 0:
            log2fc = 0.0
            expr_fdr = 1.0
            expr_dir = "unchanged"
            base_mean = None
        else:
            row = gene_expr.iloc[0]
            log2fc = row.get("log2FC", 0)
            expr_fdr = row.get("FDR", 1)
            base_mean = row.get("baseMean")

            if expr_fdr < self.expr_fdr and abs(log2fc) > self.expr_lfc:
                expr_dir = "up" if log2fc > 0 else "down"
            else:
                expr_dir = "unchanged"

        # Get histone mark data
        histone_marks = {}
        for mark_name, peak_df in histone_data.items():
            gene_peaks = peak_df[peak_df["gene_symbol"] == gene]

            if len(gene_peaks) > 0:
                # Get most significant peak
                if "FDR" in gene_peaks.columns:
                    best_peak = gene_peaks.loc[gene_peaks["FDR"].idxmin()]
                else:
                    best_peak = gene_peaks.iloc[0]

                mark_fdr = best_peak.get("FDR", best_peak.get("fdr", 1))
                mark_lfc = best_peak.get("log2FC", best_peak.get("Fold", 0))

                if mark_fdr < self.peak_fdr:
                    direction = "gained" if mark_lfc > 0 else "lost"
                else:
                    direction = "unchanged"

                histone_marks[mark_name] = {
                    "direction": direction,
                    "fdr": mark_fdr,
                    "log2fc": mark_lfc,
                    "n_peaks": len(gene_peaks),
                }

        # Get TF binding
        tf_binding = []
        if tf_data is not None and "gene_symbol" in tf_data.columns:
            gene_tf = tf_data[tf_data["gene_symbol"] == gene]
            if "tf_name" in gene_tf.columns:
                tf_binding = gene_tf["tf_name"].unique().tolist()

        # Calculate concordance
        concordance, category = self._calculate_concordance(expr_dir, histone_marks)

        # Determine chromatin state
        chromatin_state = self._determine_chromatin_state(histone_marks)

        # Get gene_id
        gene_id = ""
        if len(gene_expr) > 0 and "gene_id" in gene_expr.columns:
            gene_id = gene_expr.iloc[0]["gene_id"]

        return IntegratedGene(
            gene_id=gene_id,
            gene_symbol=gene,
            log2fc=log2fc,
            expression_fdr=expr_fdr,
            expression_direction=expr_dir,
            base_mean=base_mean,
            histone_marks=histone_marks,
            tf_binding=tf_binding,
            chromatin_state=chromatin_state,
            concordance_score=concordance,
            regulatory_category=category,
        )

    def _calculate_concordance(self, expr_direction: str, histone_marks: Dict[str, Dict]) -> Tuple[float, str]:
        """Calculate concordance between expression and marks."""
        if not histone_marks or expr_direction == "unchanged":
            return 0.0, "no_change"

        concordant_count = 0
        total_marks = 0

        for mark, data in histone_marks.items():
            if data["direction"] == "unchanged":
                continue

            total_marks += 1
            expected = self.MARK_EXPRESSION_LOGIC.get(mark, {})
            expected_expr = expected.get(data["direction"])

            if expected_expr == expr_direction:
                concordant_count += 1

        if total_marks == 0:
            return 0.0, "marks_unchanged"

        concordance = concordant_count / total_marks

        if concordance >= 0.7:
            if expr_direction == "up":
                category = "concordant_activation"
            else:
                category = "concordant_repression"
        elif concordance <= 0.3:
            category = "discordant"
        else:
            category = "mixed"

        return concordance, category

    def _determine_chromatin_state(self, histone_marks: Dict[str, Dict]) -> str:
        """Determine chromatin state from histone marks."""
        has_ac = histone_marks.get("H3K27ac", {}).get("direction") == "gained"
        has_me3 = histone_marks.get("H3K27me3", {}).get("direction") == "gained"
        has_me1 = histone_marks.get("H3K4me1", {}).get("direction") == "gained"
        has_k4me3 = histone_marks.get("H3K4me3", {}).get("direction") == "gained"

        lost_me3 = histone_marks.get("H3K27me3", {}).get("direction") == "lost"
        lost_ac = histone_marks.get("H3K27ac", {}).get("direction") == "lost"

        if has_ac and has_me1 and lost_me3:
            return "strong_activation"
        elif lost_ac and has_me3:
            return "strong_repression"
        elif has_ac and has_me1:
            return "active_enhancer"
        elif has_ac or has_k4me3:
            return "active"
        elif has_me3:
            return "repressed"
        elif has_me1:
            return "poised"
        else:
            return "unknown"

    def _summarize_integration(self, genes: List[IntegratedGene], expr_df: pd.DataFrame) -> IntegrationResult:
        """Summarize integration results."""
        # Count categories
        categories = [g.regulatory_category for g in genes]

        concordant_act = categories.count("concordant_activation")
        concordant_rep = categories.count("concordant_repression")
        discordant = categories.count("discordant")

        # Count DE genes
        de_genes = len(expr_df[(expr_df["FDR"] < self.expr_fdr) & (np.abs(expr_df["log2FC"]) > self.expr_lfc)])

        # Create summary DataFrame
        summary_data = []
        for g in genes:
            row = {
                "gene_id": g.gene_id,
                "gene_symbol": g.gene_symbol,
                "log2FC": g.log2fc,
                "expression_FDR": g.expression_fdr,
                "expression_direction": g.expression_direction,
                "chromatin_state": g.chromatin_state,
                "concordance_score": g.concordance_score,
                "regulatory_category": g.regulatory_category,
                "n_tf_bound": len(g.tf_binding),
                "tf_binding": ";".join(g.tf_binding) if g.tf_binding else "",
            }

            # Add histone mark columns
            for mark in ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3"]:
                if mark in g.histone_marks:
                    row[f"{mark}_direction"] = g.histone_marks[mark]["direction"]
                    row[f"{mark}_log2FC"] = g.histone_marks[mark]["log2fc"]
                else:
                    row[f"{mark}_direction"] = "no_data"
                    row[f"{mark}_log2FC"] = np.nan

            summary_data.append(row)

        summary_df = pd.DataFrame(summary_data)

        return IntegrationResult(
            total_genes=len(genes),
            de_genes=de_genes,
            integrated_genes=len([g for g in genes if g.regulatory_category != "no_change"]),
            concordant_activation=concordant_act,
            concordant_repression=concordant_rep,
            discordant=discordant,
            genes=genes,
            summary_df=summary_df,
        )

    def get_key_targets(
        self, result: IntegrationResult, min_concordance: float = 0.7, require_tf: bool = False
    ) -> pd.DataFrame:
        """
        Get high-confidence regulatory targets.

        Args:
            result: Integration result
            min_concordance: Minimum concordance score
            require_tf: Require TF binding evidence

        Returns:
            DataFrame of key targets
        """
        targets = []

        for gene in result.genes:
            if gene.concordance_score < min_concordance:
                continue
            if require_tf and not gene.tf_binding:
                continue
            if gene.expression_direction == "unchanged":
                continue

            targets.append(gene)

        # Convert to DataFrame
        if not targets:
            return pd.DataFrame()

        data = [
            {
                "gene_symbol": g.gene_symbol,
                "log2FC": g.log2fc,
                "expression_FDR": g.expression_fdr,
                "direction": g.expression_direction,
                "chromatin_state": g.chromatin_state,
                "concordance": g.concordance_score,
                "tf_binding": ";".join(g.tf_binding),
                "category": g.regulatory_category,
            }
            for g in targets
        ]

        return pd.DataFrame(data).sort_values("expression_FDR")


# Convenience functions
def load_expression(filepath: str, **kwargs) -> ExpressionData:
    """Load expression data from file."""
    loader = ExpressionLoader()
    return loader.load(filepath, **kwargs)


def integrate_expression_with_epigenetics(
    expression_file: str,
    histone_files: Dict[str, str],
    tf_file: Optional[str] = None,
    output_file: Optional[str] = None,
) -> IntegrationResult:
    """
    Convenience function for full integration analysis.

    Args:
        expression_file: Path to DE results
        histone_files: Dict of mark name -> peak file path
        tf_file: Optional TF binding file
        output_file: Optional output file path

    Returns:
        IntegrationResult
    """
    # Load expression
    loader = ExpressionLoader()
    expr_data = loader.load(expression_file)

    # Load histone data
    histone_data = {}
    for mark, filepath in histone_files.items():
        histone_data[mark] = pd.read_csv(filepath)

    # Load TF data
    tf_data = None
    if tf_file:
        tf_data = pd.read_csv(tf_file)

    # Integrate
    integrator = ExpressionHistoneIntegrator()
    result = integrator.integrate(expr_data, histone_data, tf_data)

    # Save if requested
    if output_file and result.summary_df is not None:
        result.summary_df.to_csv(output_file, index=False)
        logger.info(f"Saved integration results to {output_file}")

    return result
