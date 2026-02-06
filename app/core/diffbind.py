"""
DiffBind Analysis Runner

Wrapper for running DiffBind differential peak analysis via R.
Generates sample sheets, executes R scripts, and parses results.
"""

import os
import subprocess
import tempfile
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class DiffBindConfig:
    """Configuration for DiffBind analysis."""
    comparison_name: str
    group1: str
    group2: str
    histone_mark: str
    genome: str = "mm10"

    # Analysis parameters
    fdr_threshold: float = 0.1
    lfc_threshold: float = 0.5
    min_overlap: int = 1
    summit_size: int = 250

    # Normalization
    normalize_method: str = "RLE"  # RLE, TMM, or lib
    use_background: bool = True

    # Gene annotation
    tss_upstream: int = 3000
    tss_downstream: int = 3000

    # Output options
    output_dir: Optional[str] = None
    generate_plots: bool = True


@dataclass
class DiffBindResults:
    """Results from DiffBind analysis."""
    comparison_name: str
    total_peaks: int
    significant_peaks: int
    gained_peaks: int
    lost_peaks: int

    # File paths
    all_peaks_file: str
    significant_peaks_file: str
    genes_file: str
    plots_file: Optional[str]

    # Summary statistics
    summary: Dict[str, Any]


class DiffBindRunner:
    """Runner for DiffBind differential peak analysis."""

    def __init__(self, r_path: str = "Rscript", r_libs: Optional[str] = None):
        """
        Initialize DiffBind runner.

        Args:
            r_path: Path to Rscript executable
            r_libs: Optional path to R library directory
        """
        self.r_path = r_path
        self.r_libs = r_libs
        self.script_dir = Path(__file__).parent.parent.parent / "scripts" / "r_scripts"

    def create_sample_sheet(
        self,
        samples: List[Dict],
        output_path: str
    ) -> str:
        """
        Create DiffBind-compatible sample sheet.

        Args:
            samples: List of sample dictionaries with keys:
                     SampleID, Condition, Factor, Replicate, bamReads, Peaks
            output_path: Path to save CSV

        Returns:
            Path to created sample sheet
        """
        required_cols = ["SampleID", "Condition", "Factor", "Replicate", "bamReads", "Peaks"]

        # Validate samples
        for sample in samples:
            missing = [c for c in required_cols if c not in sample]
            if missing:
                raise ValueError(f"Sample missing columns: {missing}")

        df = pd.DataFrame(samples)

        # Add optional columns with defaults
        if "Treatment" not in df.columns:
            df["Treatment"] = df["Condition"]
        if "ControlID" not in df.columns:
            df["ControlID"] = ""
        if "PeakCaller" not in df.columns:
            df["PeakCaller"] = "bed"

        df.to_csv(output_path, index=False)
        logger.info(f"Created sample sheet: {output_path}")

        return output_path

    def generate_r_script(self, config: DiffBindConfig) -> str:
        """
        Generate R script for DiffBind analysis.

        Args:
            config: DiffBind configuration

        Returns:
            R script as string
        """
        # Determine TxDb and OrgDb based on genome
        txdb_map = {
            "mm10": "TxDb.Mmusculus.UCSC.mm10.knownGene",
            "mm39": "TxDb.Mmusculus.UCSC.mm39.knownGene",
            "hg38": "TxDb.Hsapiens.UCSC.hg38.knownGene",
            "hg19": "TxDb.Hsapiens.UCSC.hg19.knownGene",
            "rn6": "TxDb.Rnorvegicus.UCSC.rn6.refGene",
            "rn7": "TxDb.Rnorvegicus.UCSC.rn7.refGene",
        }

        orgdb_map = {
            "mm10": "org.Mm.eg.db", "mm39": "org.Mm.eg.db",
            "hg38": "org.Hs.eg.db", "hg19": "org.Hs.eg.db",
            "rn6": "org.Rn.eg.db", "rn7": "org.Rn.eg.db",
        }

        txdb = txdb_map.get(config.genome, txdb_map["mm10"])
        orgdb = orgdb_map.get(config.genome, orgdb_map["mm10"])

        # Normalization method mapping
        norm_map = {
            "RLE": "DBA_NORM_RLE",
            "TMM": "DBA_NORM_TMM",
            "lib": "DBA_NORM_LIB"
        }
        norm_method = norm_map.get(config.normalize_method, "DBA_NORM_RLE")

        script = f'''#!/usr/bin/env Rscript
# DiffBind Analysis Script - Auto-generated
# Comparison: {config.comparison_name}

suppressPackageStartupMessages({{
    library(DiffBind)
    library(csaw)
    library(ChIPseeker)
    library({txdb})
    library({orgdb})
    library(GenomicRanges)
    library(dplyr)
    library(ggplot2)
}})

# Configuration
config <- list(
    samplesheet = Sys.getenv("SAMPLESHEET"),
    output_dir = Sys.getenv("OUTPUT_DIR"),
    comparison_name = "{config.comparison_name}",
    group1 = "{config.group1}",
    group2 = "{config.group2}",
    fdr_threshold = {config.fdr_threshold},
    lfc_threshold = {config.lfc_threshold},
    min_overlap = {config.min_overlap},
    summit_size = {config.summit_size},
    tss_region = c(-{config.tss_upstream}, {config.tss_downstream})
)

cat("\\n=== DiffBind Analysis ===\\n")
cat(sprintf("Comparison: %s\\n", config$comparison_name))
cat(sprintf("Groups: %s vs %s\\n\\n", config$group1, config$group2))

# Create output directory
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load samples
cat("Loading sample sheet...\\n")
dba_obj <- dba(sampleSheet = config$samplesheet)
print(dba_obj)

# Count reads
cat("\\nCounting reads in peaks...\\n")
dba_count <- dba.count(
    dba_obj,
    minOverlap = config$min_overlap,
    summits = config$summit_size,
    bParallel = TRUE
)

# Normalize
cat("\\nNormalizing...\\n")
dba_norm <- dba.normalize(
    dba_count,
    normalize = {norm_method},
    library = DBA_LIBSIZE_FULL,
    background = {"TRUE" if config.use_background else "FALSE"}
)

# Set up contrast
cat("\\nSetting up contrast...\\n")
dba_contrast <- dba.contrast(
    dba_norm,
    contrast = c("Condition", config$group1, config$group2),
    minMembers = 2
)

# Run differential analysis
cat("\\nRunning differential analysis (DESeq2)...\\n")
dba_analyzed <- dba.analyze(dba_contrast, method = DBA_DESEQ2, bParallel = TRUE)

# Get results
cat("\\nExtracting results...\\n")
db_report <- dba.report(dba_analyzed, th = 1, bCounts = TRUE)
results_df <- as.data.frame(db_report)
results_df$peak_id <- sprintf("peak_%05d", seq_len(nrow(results_df)))

# Annotate peaks
cat("\\nAnnotating peaks to genes...\\n")
peaks_gr <- GRanges(
    seqnames = results_df$seqnames,
    ranges = IRanges(start = results_df$start, end = results_df$end),
    peak_id = results_df$peak_id
)

txdb <- {txdb}
peak_anno <- annotatePeak(peaks_gr, tssRegion = config$tss_region,
                          TxDb = txdb, annoDb = "{orgdb}")
anno_df <- as.data.frame(peak_anno)

# Merge results with annotation
results_annotated <- results_df %>%
    left_join(
        anno_df %>% dplyr::select(
            peak_id, annotation, geneId = ENSEMBL,
            gene_name = SYMBOL, distanceToTSS
        ),
        by = "peak_id"
    ) %>%
    mutate(
        significant = FDR < config$fdr_threshold,
        direction = ifelse(Fold > 0, "Gained", "Lost")
    )

# Save all peaks
all_peaks_file <- file.path(config$output_dir, "all_peaks_annotated.csv")
write.csv(results_annotated, all_peaks_file, row.names = FALSE)
cat(sprintf("Saved all peaks: %s\\n", all_peaks_file))

# Save significant peaks
sig_peaks <- results_annotated %>% filter(significant)
sig_peaks_file <- file.path(config$output_dir, "significant_peaks.csv")
write.csv(sig_peaks, sig_peaks_file, row.names = FALSE)
cat(sprintf("Saved significant peaks: %s\\n", sig_peaks_file))

# Create gene summary
genes_summary <- results_annotated %>%
    filter(!is.na(gene_name)) %>%
    group_by(gene_name, geneId) %>%
    summarise(
        n_peaks = n(),
        n_sig_peaks = sum(significant),
        direction = paste(unique(direction[significant]), collapse = ","),
        min_FDR = min(FDR),
        max_abs_Fold = max(abs(Fold)),
        .groups = "drop"
    ) %>%
    arrange(min_FDR)

genes_file <- file.path(config$output_dir, "genes_summary.csv")
write.csv(genes_summary, genes_file, row.names = FALSE)
cat(sprintf("Saved gene summary: %s\\n", genes_file))

# Generate plots
cat("\\nGenerating plots...\\n")
plots_file <- file.path(config$output_dir, "analysis_plots.pdf")
pdf(plots_file, width = 10, height = 8)

# MA plot
tryCatch({{
    dba.plotMA(dba_analyzed, th = config$fdr_threshold)
    title(main = sprintf("MA Plot: %s", config$comparison_name))
}}, error = function(e) cat("MA plot error:", e$message, "\\n"))

# Volcano plot
tryCatch({{
    dba.plotVolcano(dba_analyzed, th = config$fdr_threshold)
    title(main = sprintf("Volcano Plot: %s", config$comparison_name))
}}, error = function(e) cat("Volcano plot error:", e$message, "\\n"))

# PCA
tryCatch({{
    dba.plotPCA(dba_norm, label = DBA_CONDITION)
    title(main = "PCA: All Samples")
}}, error = function(e) cat("PCA plot error:", e$message, "\\n"))

# Correlation heatmap
tryCatch({{
    dba.plotHeatmap(dba_norm, correlations = TRUE)
}}, error = function(e) cat("Heatmap error:", e$message, "\\n"))

dev.off()
cat(sprintf("Saved plots: %s\\n", plots_file))

# Summary statistics
summary_stats <- list(
    comparison = config$comparison_name,
    total_peaks = nrow(results_annotated),
    significant_peaks = sum(results_annotated$significant),
    gained = sum(results_annotated$significant & results_annotated$direction == "Gained"),
    lost = sum(results_annotated$significant & results_annotated$direction == "Lost"),
    genes_with_peaks = n_distinct(results_annotated$gene_name[!is.na(results_annotated$gene_name)]),
    genes_with_sig_peaks = n_distinct(sig_peaks$gene_name[!is.na(sig_peaks$gene_name)])
)

# Save summary as JSON
summary_file <- file.path(config$output_dir, "summary.json")
writeLines(jsonlite::toJSON(summary_stats, auto_unbox = TRUE, pretty = TRUE), summary_file)
cat(sprintf("Saved summary: %s\\n", summary_file))

cat("\\n=== Analysis Complete ===\\n")
cat(sprintf("Total peaks: %d\\n", summary_stats$total_peaks))
cat(sprintf("Significant peaks (FDR < %.2f): %d\\n", config$fdr_threshold, summary_stats$significant_peaks))
cat(sprintf("  Gained: %d\\n", summary_stats$gained))
cat(sprintf("  Lost: %d\\n", summary_stats$lost))
'''

        return script

    def run(
        self,
        samples: List[Dict],
        config: DiffBindConfig,
        output_dir: Optional[str] = None
    ) -> DiffBindResults:
        """
        Run DiffBind analysis.

        Args:
            samples: List of sample dictionaries
            config: Analysis configuration
            output_dir: Output directory (overrides config.output_dir)

        Returns:
            DiffBindResults object
        """
        # Set output directory
        if output_dir:
            config.output_dir = output_dir
        elif not config.output_dir:
            config.output_dir = tempfile.mkdtemp(prefix="diffbind_")

        output_path = Path(config.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Create sample sheet
        samplesheet_path = output_path / "samplesheet.csv"
        self.create_sample_sheet(samples, str(samplesheet_path))

        # Generate R script
        script_content = self.generate_r_script(config)
        script_path = output_path / "diffbind_analysis.R"
        script_path.write_text(script_content)

        # Set environment variables
        env = os.environ.copy()
        env["SAMPLESHEET"] = str(samplesheet_path)
        env["OUTPUT_DIR"] = str(output_path)
        if self.r_libs:
            env["R_LIBS"] = self.r_libs

        # Run R script
        logger.info(f"Running DiffBind analysis: {config.comparison_name}")

        try:
            result = subprocess.run(
                [self.r_path, str(script_path)],
                env=env,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            if result.returncode != 0:
                logger.error(f"DiffBind failed: {result.stderr}")
                raise RuntimeError(f"DiffBind analysis failed: {result.stderr}")

            logger.info("DiffBind analysis completed successfully")

        except subprocess.TimeoutExpired:
            raise RuntimeError("DiffBind analysis timed out")

        # Parse results
        summary_file = output_path / "summary.json"
        if summary_file.exists():
            with open(summary_file) as f:
                summary = json.load(f)
        else:
            summary = {}

        return DiffBindResults(
            comparison_name=config.comparison_name,
            total_peaks=summary.get("total_peaks", 0),
            significant_peaks=summary.get("significant_peaks", 0),
            gained_peaks=summary.get("gained", 0),
            lost_peaks=summary.get("lost", 0),
            all_peaks_file=str(output_path / "all_peaks_annotated.csv"),
            significant_peaks_file=str(output_path / "significant_peaks.csv"),
            genes_file=str(output_path / "genes_summary.csv"),
            plots_file=str(output_path / "analysis_plots.pdf") if config.generate_plots else None,
            summary=summary
        )

    def run_from_samplesheet(
        self,
        samplesheet_path: str,
        config: DiffBindConfig
    ) -> DiffBindResults:
        """
        Run DiffBind from an existing sample sheet.

        Args:
            samplesheet_path: Path to DiffBind sample sheet CSV
            config: Analysis configuration

        Returns:
            DiffBindResults object
        """
        samples = pd.read_csv(samplesheet_path).to_dict("records")
        return self.run(samples, config)
