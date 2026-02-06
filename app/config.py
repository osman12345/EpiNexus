"""
Configuration settings for the Histone Mark Analyzer.

Supports multiple reference genomes, configurable tool paths,
and flexible deployment options.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional
from pydantic_settings import BaseSettings
from pydantic import Field


class GenomeConfig:
    """Reference genome configuration."""

    SUPPORTED_GENOMES = {
        "mm10": {
            "name": "Mouse (mm10/GRCm38)",
            "species": "Mus musculus",
            "txdb": "TxDb.Mmusculus.UCSC.mm10.knownGene",
            "orgdb": "org.Mm.eg.db",
            "chrom_sizes": "mm10.chrom.sizes"
        },
        "mm39": {
            "name": "Mouse (mm39/GRCm39)",
            "species": "Mus musculus",
            "txdb": "TxDb.Mmusculus.UCSC.mm39.knownGene",
            "orgdb": "org.Mm.eg.db",
            "chrom_sizes": "mm39.chrom.sizes"
        },
        "hg38": {
            "name": "Human (hg38/GRCh38)",
            "species": "Homo sapiens",
            "txdb": "TxDb.Hsapiens.UCSC.hg38.knownGene",
            "orgdb": "org.Hs.eg.db",
            "chrom_sizes": "hg38.chrom.sizes"
        },
        "hg19": {
            "name": "Human (hg19/GRCh37)",
            "species": "Homo sapiens",
            "txdb": "TxDb.Hsapiens.UCSC.hg19.knownGene",
            "orgdb": "org.Hs.eg.db",
            "chrom_sizes": "hg19.chrom.sizes"
        },
        "rn6": {
            "name": "Rat (rn6)",
            "species": "Rattus norvegicus",
            "txdb": "TxDb.Rnorvegicus.UCSC.rn6.refGene",
            "orgdb": "org.Rn.eg.db",
            "chrom_sizes": "rn6.chrom.sizes"
        },
        "rn7": {
            "name": "Rat (rn7/mRatBN7.2)",
            "species": "Rattus norvegicus",
            "txdb": "TxDb.Rnorvegicus.UCSC.rn7.refGene",
            "orgdb": "org.Rn.eg.db",
            "chrom_sizes": "rn7.chrom.sizes"
        }
    }

    HISTONE_MARKS = {
        "H3K27ac": {
            "type": "active",
            "description": "Active enhancers and promoters",
            "peak_type": "sharp",
            "color": "#4DAF4A"
        },
        "H3K27me3": {
            "type": "repressive",
            "description": "Polycomb repression",
            "peak_type": "broad",
            "color": "#E41A1C"
        },
        "H3K4me1": {
            "type": "priming",
            "description": "Enhancer priming/poised",
            "peak_type": "broad",
            "color": "#FF7F00"
        },
        "H3K4me3": {
            "type": "active",
            "description": "Active promoters",
            "peak_type": "sharp",
            "color": "#377EB8"
        },
        "H3K9me3": {
            "type": "repressive",
            "description": "Heterochromatin",
            "peak_type": "broad",
            "color": "#984EA3"
        },
        "H3K36me3": {
            "type": "active",
            "description": "Transcribed gene bodies",
            "peak_type": "broad",
            "color": "#A65628"
        }
    }


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # Application
    app_name: str = "Histone Mark Analyzer"
    app_version: str = "0.1.0"
    debug: bool = False

    # Paths
    base_dir: Path = Field(default_factory=lambda: Path(__file__).parent.parent)
    data_dir: Path = Field(default_factory=lambda: Path(__file__).parent.parent / "data")
    results_dir: Path = Field(default_factory=lambda: Path(__file__).parent.parent / "results")
    uploads_dir: Path = Field(default_factory=lambda: Path(__file__).parent.parent / "data" / "uploads")
    references_dir: Path = Field(default_factory=lambda: Path(__file__).parent.parent / "data" / "references")

    # Database
    database_url: str = "sqlite:///./histone_analyzer.db"

    # Redis (for background jobs)
    redis_url: str = "redis://localhost:6379/0"

    # External Tools
    bowtie2_path: str = "bowtie2"
    samtools_path: str = "samtools"
    bedtools_path: str = "bedtools"
    picard_path: str = "picard"
    seacr_path: str = "SEACR_1.3.sh"
    macs2_path: str = "macs2"
    deeptools_path: str = "deepTools"

    # R configuration
    r_path: str = "Rscript"
    r_libs_path: Optional[str] = None

    # Default analysis parameters
    default_genome: str = "mm10"
    default_fdr_threshold: float = 0.1
    default_lfc_threshold: float = 0.5
    default_tss_region: tuple = (-3000, 3000)

    # Peak calling defaults
    seacr_threshold: float = 0.01
    seacr_mode: str = "stringent"  # or "relaxed"
    macs2_qvalue: float = 0.05

    # Server
    api_host: str = "0.0.0.0"
    api_port: int = 8000
    streamlit_port: int = 8501

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    def ensure_directories(self):
        """Create necessary directories if they don't exist."""
        for dir_path in [self.data_dir, self.results_dir, self.uploads_dir, self.references_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def get_genome_config(self, genome: str) -> Dict:
        """Get configuration for a specific genome."""
        if genome not in GenomeConfig.SUPPORTED_GENOMES:
            raise ValueError(f"Unsupported genome: {genome}. Supported: {list(GenomeConfig.SUPPORTED_GENOMES.keys())}")
        return GenomeConfig.SUPPORTED_GENOMES[genome]

    def get_histone_config(self, mark: str) -> Dict:
        """Get configuration for a specific histone mark."""
        if mark not in GenomeConfig.HISTONE_MARKS:
            raise ValueError(f"Unknown histone mark: {mark}. Known marks: {list(GenomeConfig.HISTONE_MARKS.keys())}")
        return GenomeConfig.HISTONE_MARKS[mark]

    def get_reference_path(self, genome: str) -> Path:
        """Get path to reference genome bowtie2 index."""
        return self.references_dir / genome / "genome"

    def get_chrom_sizes_path(self, genome: str) -> Path:
        """Get path to chromosome sizes file."""
        config = self.get_genome_config(genome)
        return self.references_dir / genome / config["chrom_sizes"]


# Global settings instance
settings = Settings()
