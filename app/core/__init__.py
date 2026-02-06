"""
Core analysis modules for the Histone Mark Analyzer.

Pure Python implementation - no R required!

Includes:
- Differential peak analysis (PyDESeq2)
- Peak annotation
- Multi-mark integration
- Transcription factor analysis
- Quality control
"""

# Differential analysis
from .differential import DifferentialAnalyzer, DifferentialConfig, DifferentialResults

# Peak annotation
from .annotation import PeakAnnotator, GeneAnnotation, AnnotationDatabase

# Multi-mark integration
from .integration import MarkIntegration, IntegrationConfig, IntegrationResults

# Quality control
from .qc import QCAnalyzer, QCMetrics

# Transcription factor analysis
from .transcription_factors import (
    JASPARDatabase,
    PositionWeightMatrix,
    MotifScanner,
    MotifEnrichmentAnalyzer,
    MotifEnrichmentResult,
    TFTargetPredictor,
    TFTargetPrediction,
    TFHistoneIntegrator,
    run_tf_motif_analysis
)

# Legacy R wrapper (kept for compatibility)
from .diffbind import DiffBindRunner

__all__ = [
    # Differential analysis
    "DifferentialAnalyzer",
    "DifferentialConfig",
    "DifferentialResults",

    # Annotation
    "PeakAnnotator",
    "GeneAnnotation",
    "AnnotationDatabase",

    # Integration
    "MarkIntegration",
    "IntegrationConfig",
    "IntegrationResults",

    # Quality control
    "QCAnalyzer",
    "QCMetrics",

    # TF Analysis
    "JASPARDatabase",
    "PositionWeightMatrix",
    "MotifScanner",
    "MotifEnrichmentAnalyzer",
    "MotifEnrichmentResult",
    "TFTargetPredictor",
    "TFTargetPrediction",
    "TFHistoneIntegrator",
    "run_tf_motif_analysis",

    # Legacy
    "DiffBindRunner",
]
