"""
Sample Project Generator for New User Onboarding

Generates a realistic but small demo dataset that allows new users to
explore every EpiNexus feature without needing their own data.  The
generated project includes:

- Sample sheet with treatment / control groups
- Simulated ChIP-seq peak files (narrowPeak format)
- Differential analysis results (CSV)
- A gene expression table for integration
- A mini methylation dataset for the methylation module

All data is synthetic but modelled on realistic distributions
(peak widths, fold-change distributions, bimodal methylation, etc.)
so that the UI plots look representative.

Usage::

    from app.core.sample_project import create_sample_project

    project_dir = create_sample_project("/path/to/output")
    # project_dir now contains a ready-to-load EpiNexus project

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def create_sample_project(
    output_dir: str,
    genome: str = "hg38",
    seed: int = 42,
) -> Path:
    """Generate a complete sample project ready for EpiNexus import.

    Args:
        output_dir: Root directory for the generated project.
        genome: Reference genome label (used in filenames/metadata).
        seed: Random seed for reproducibility.

    Returns:
        ``Path`` to the created project directory.
    """
    rng = np.random.default_rng(seed)
    root = Path(output_dir) / "EpiNexus_Sample_Project"
    root.mkdir(parents=True, exist_ok=True)

    logger.info("Generating sample project in %s", root)

    # 1. Sample sheet
    _write_sample_sheet(root, genome)

    # 2. Peak files
    peaks_dir = root / "peaks"
    peaks_dir.mkdir(exist_ok=True)
    for sample in ["Treatment_Rep1", "Treatment_Rep2", "Control_Rep1", "Control_Rep2"]:
        _write_narrowpeak(peaks_dir / f"{sample}_H3K27ac.narrowPeak", rng, sample)

    # 3. Differential results
    _write_differential_results(root / "differential_H3K27ac.csv", rng)

    # 4. Expression table
    _write_expression_table(root / "expression_deseq2.csv", rng)

    # 5. Methylation data
    _write_methylation_data(root / "methylation_control.cov", rng, bias=1.0)
    _write_methylation_data(root / "methylation_treatment.cov", rng, bias=0.7)

    # 6. README
    _write_readme(root, genome)

    logger.info("Sample project created: %s", root)
    return root


# ---------------------------------------------------------------------------
# Internal generators
# ---------------------------------------------------------------------------

_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
_CHR_WEIGHTS = np.array(
    [8, 7, 6, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3],
    dtype=float,
)
_CHR_WEIGHTS /= _CHR_WEIGHTS.sum()


def _write_sample_sheet(root: Path, genome: str) -> None:
    """Write a sample sheet CSV."""
    rows = [
        {"sample_name": "Treatment_Rep1", "condition": "Treatment", "replicate": 1,
         "histone_mark": "H3K27ac", "genome": genome,
         "peak_file": "peaks/Treatment_Rep1_H3K27ac.narrowPeak"},
        {"sample_name": "Treatment_Rep2", "condition": "Treatment", "replicate": 2,
         "histone_mark": "H3K27ac", "genome": genome,
         "peak_file": "peaks/Treatment_Rep2_H3K27ac.narrowPeak"},
        {"sample_name": "Control_Rep1", "condition": "Control", "replicate": 1,
         "histone_mark": "H3K27ac", "genome": genome,
         "peak_file": "peaks/Control_Rep1_H3K27ac.narrowPeak"},
        {"sample_name": "Control_Rep2", "condition": "Control", "replicate": 2,
         "histone_mark": "H3K27ac", "genome": genome,
         "peak_file": "peaks/Control_Rep2_H3K27ac.narrowPeak"},
    ]
    pd.DataFrame(rows).to_csv(root / "sample_sheet.csv", index=False)


def _write_narrowpeak(path: Path, rng: np.random.Generator, sample: str) -> None:
    """Generate a synthetic narrowPeak file (~15 000 peaks)."""
    n = 15_000
    chroms = rng.choice(_CHROMS, size=n, p=_CHR_WEIGHTS)
    starts = rng.integers(1_000_000, 200_000_000, size=n)
    widths = np.clip(rng.normal(500, 150, n).astype(int), 100, 2000)
    ends = starts + widths
    scores = rng.integers(100, 1000, size=n)
    signals = np.round(rng.exponential(5, n) + 1, 2)
    pvals = np.round(-np.log10(rng.uniform(1e-10, 0.05, n)), 2)
    qvals = np.round(pvals * rng.uniform(0.5, 1.0, n), 2)
    summits = widths // 2

    df = pd.DataFrame({
        "chrom": chroms, "start": starts, "end": ends,
        "name": [f"{sample}_peak_{i}" for i in range(n)],
        "score": scores, "strand": ".",
        "signalValue": signals, "pValue": pvals,
        "qValue": qvals, "peak": summits,
    })
    df.sort_values(["chrom", "start"], inplace=True)
    df.to_csv(path, sep="\t", header=False, index=False)


def _write_differential_results(path: Path, rng: np.random.Generator) -> None:
    """Generate simulated differential peak analysis results."""
    n = 25_000
    log2fc = rng.normal(0, 1.2, n)
    pval = 10 ** (-np.abs(log2fc) * rng.uniform(0.5, 2.5, n))
    pval = np.clip(pval, 1e-300, 1.0)
    fdr = np.minimum(pval * n / (np.argsort(np.argsort(pval)) + 1), 1.0)

    chroms = rng.choice(_CHROMS, size=n, p=_CHR_WEIGHTS)
    starts = rng.integers(1_000_000, 200_000_000, size=n)
    ends = starts + rng.integers(200, 1500, size=n)

    direction = np.where(
        fdr < 0.05,
        np.where(log2fc > 0, "Up", "Down"),
        "Not Significant",
    )

    pd.DataFrame({
        "chr": chroms, "start": starts, "end": ends,
        "log2FoldChange": np.round(log2fc, 4),
        "pvalue": pval, "FDR": fdr,
        "direction": direction,
    }).to_csv(path, index=False)


def _write_expression_table(path: Path, rng: np.random.Generator) -> None:
    """Generate simulated DESeq2 expression output."""
    n = 8_000
    gene_ids = [f"ENSG{i:011d}" for i in range(n)]
    symbols = [f"Gene{i}" for i in range(n)]

    log2fc = rng.normal(0, 1.5, n)
    basemean = 10 ** rng.normal(2, 1, n)
    pval = 10 ** (-np.abs(log2fc) * rng.uniform(0.3, 2, n))
    pval = np.clip(pval, 1e-300, 1.0)
    fdr = np.minimum(pval * n / (np.argsort(np.argsort(pval)) + 1), 1.0)

    pd.DataFrame({
        "gene_id": gene_ids, "gene_name": symbols,
        "baseMean": np.round(basemean, 2),
        "log2FoldChange": np.round(log2fc, 4),
        "pvalue": pval, "padj": fdr,
    }).to_csv(path, index=False)


def _write_methylation_data(
    path: Path, rng: np.random.Generator, bias: float = 1.0
) -> None:
    """Generate a Bismark-style .cov methylation file.

    Args:
        path: Output file path.
        rng: Numpy random generator.
        bias: Multiplicative bias applied to methylation % (<1 → hypomethylated).
    """
    n = 20_000
    chroms = rng.choice(_CHROMS, size=n, p=_CHR_WEIGHTS)
    starts = rng.integers(1_000_000, 200_000_000, size=n)
    ends = starts + 1

    # Bimodal methylation (typical CpG distribution)
    high_frac = 0.7
    meth = np.where(
        rng.random(n) < high_frac,
        np.clip(rng.beta(8, 2, n) * 100 * bias, 0, 100),
        np.clip(rng.beta(2, 8, n) * 100 * bias, 0, 100),
    )

    coverage = np.clip(rng.exponential(15, n).astype(int), 1, 200)
    count_meth = np.round(coverage * meth / 100).astype(int)
    count_unmeth = coverage - count_meth

    pd.DataFrame({
        "chr": chroms, "start": starts, "end": ends,
        "meth_pct": np.round(meth, 2),
        "count_meth": count_meth,
        "count_unmeth": count_unmeth,
    }).sort_values(["chr", "start"]).to_csv(path, sep="\t", header=False, index=False)


def _write_readme(root: Path, genome: str) -> None:
    """Write a project README."""
    (root / "README.txt").write_text(
        f"""EpiNexus Sample Project
======================

Genome assembly : {genome}
Generated by    : EpiNexus sample_project module

Contents
--------
sample_sheet.csv                  - Sample sheet (4 samples, 2 conditions)
peaks/                            - Simulated narrowPeak files (H3K27ac)
differential_H3K27ac.csv          - Differential peak analysis results
expression_deseq2.csv             - Simulated DESeq2 output (8 000 genes)
methylation_control.cov           - Bismark-style CpG methylation (control)
methylation_treatment.cov         - Bismark-style CpG methylation (treatment)

Quick Start
-----------
1. Open EpiNexus and go to Data & Project.
2. Upload ``sample_sheet.csv`` via the sample sheet uploader.
3. Explore differential analysis, expression integration,
   methylation, and more using the included data files.

Note: All data is synthetic and intended for demonstration only.
"""
    )
