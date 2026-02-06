"""
FASTQ Preprocessing Pipeline

Complete pipeline from raw FASTQ to analysis-ready files:
1. Quality control (FastQC)
2. Trimming (Trimmomatic/Cutadapt)
3. Alignment (Bowtie2/BWA)
4. Post-alignment QC
5. Peak calling (MACS2/SICER)
6. BigWig generation
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum
import json


class AssayType(Enum):
    """Supported assay types."""
    CHIP_SEQ_HISTONE = "chip_histone"
    CHIP_SEQ_TF = "chip_tf"
    ATAC_SEQ = "atac"
    CUT_AND_RUN = "cut_run"
    CUT_AND_TAG = "cut_tag"


class PeakType(Enum):
    """Peak calling modes."""
    NARROW = "narrow"  # TFs, H3K4me3, H3K27ac
    BROAD = "broad"    # H3K27me3, H3K36me3, H3K9me3


@dataclass
class PipelineConfig:
    """Configuration for the preprocessing pipeline."""
    # Input
    assay_type: AssayType
    genome: str = "hg38"
    paired_end: bool = True

    # Trimming
    trim_adapter: str = "auto"
    min_length: int = 20
    quality_threshold: int = 20

    # Alignment
    aligner: str = "bowtie2"  # or 'bwa'
    max_fragment_size: int = 1000

    # Filtering
    min_mapq: int = 30
    remove_duplicates: bool = True
    remove_blacklist: bool = True

    # Peak calling
    peak_caller: str = "macs2"  # or 'sicer', 'seacr'
    peak_type: PeakType = PeakType.NARROW
    qvalue_threshold: float = 0.01

    # Resources
    threads: int = 4
    memory_gb: int = 16


@dataclass
class SampleInfo:
    """Information about a single sample."""
    sample_id: str
    fastq_r1: str
    fastq_r2: Optional[str] = None
    control_id: Optional[str] = None  # Input/IgG control
    condition: str = ""
    replicate: int = 1
    mark_or_tf: str = ""  # e.g., "H3K27ac" or "CTCF"


@dataclass
class PipelineOutput:
    """Output files from pipeline."""
    sample_id: str
    bam_file: str = ""
    bam_filtered: str = ""
    bigwig_file: str = ""
    peaks_file: str = ""
    qc_report: str = ""
    stats: Dict = field(default_factory=dict)


class FASTQPipeline:
    """
    End-to-end ChIP-seq/ATAC-seq preprocessing pipeline.

    Workflow:
    FASTQ → QC → Trim → Align → Filter → Deduplicate → Peak Call → BigWig
    """

    # Genome index paths (configure for your system)
    GENOME_INDICES = {
        "hg38": {
            "bowtie2": "/reference/bowtie2/hg38",
            "bwa": "/reference/bwa/hg38.fa",
            "chrom_sizes": "/reference/hg38.chrom.sizes",
            "blacklist": "/reference/hg38-blacklist.v2.bed"
        },
        "hg19": {
            "bowtie2": "/reference/bowtie2/hg19",
            "bwa": "/reference/bwa/hg19.fa",
            "chrom_sizes": "/reference/hg19.chrom.sizes",
            "blacklist": "/reference/hg19-blacklist.v2.bed"
        },
        "mm10": {
            "bowtie2": "/reference/bowtie2/mm10",
            "bwa": "/reference/bwa/mm10.fa",
            "chrom_sizes": "/reference/mm10.chrom.sizes",
            "blacklist": "/reference/mm10-blacklist.v2.bed"
        }
    }

    # Peak type recommendations
    MARK_PEAK_TYPES = {
        # Narrow peaks
        "H3K4me3": PeakType.NARROW,
        "H3K27ac": PeakType.NARROW,
        "H3K9ac": PeakType.NARROW,
        "H3K4me2": PeakType.NARROW,
        # Broad peaks
        "H3K27me3": PeakType.BROAD,
        "H3K36me3": PeakType.BROAD,
        "H3K9me3": PeakType.BROAD,
        "H3K79me2": PeakType.BROAD,
        "H4K20me1": PeakType.BROAD,
        # TFs are always narrow
        "TF": PeakType.NARROW
    }

    def __init__(self, config: PipelineConfig, output_dir: str):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.dirs = {
            'fastqc': self.output_dir / 'fastqc',
            'trimmed': self.output_dir / 'trimmed',
            'aligned': self.output_dir / 'aligned',
            'filtered': self.output_dir / 'filtered',
            'peaks': self.output_dir / 'peaks',
            'bigwig': self.output_dir / 'bigwig',
            'qc': self.output_dir / 'qc'
        }
        for d in self.dirs.values():
            d.mkdir(exist_ok=True)

    def run_sample(self, sample: SampleInfo, control: SampleInfo = None) -> PipelineOutput:
        """Run complete pipeline for a single sample."""

        output = PipelineOutput(sample_id=sample.sample_id)

        # Step 1: FastQC
        self._run_fastqc(sample)

        # Step 2: Trim adapters
        trimmed_r1, trimmed_r2 = self._trim_reads(sample)

        # Step 3: Align to genome
        bam_file = self._align_reads(sample, trimmed_r1, trimmed_r2)
        output.bam_file = bam_file

        # Step 4: Filter and deduplicate
        filtered_bam = self._filter_bam(sample, bam_file)
        output.bam_filtered = filtered_bam

        # Step 5: Generate BigWig
        bigwig_file = self._generate_bigwig(sample, filtered_bam)
        output.bigwig_file = bigwig_file

        # Step 6: Call peaks
        control_bam = None
        if control:
            control_bam = str(self.dirs['filtered'] / f"{control.sample_id}.filtered.bam")

        peaks_file = self._call_peaks(sample, filtered_bam, control_bam)
        output.peaks_file = peaks_file

        # Step 7: Generate QC report
        stats = self._collect_stats(sample, output)
        output.stats = stats

        return output

    def _run_fastqc(self, sample: SampleInfo) -> str:
        """Run FastQC on raw reads."""
        cmd = [
            "fastqc",
            "-o", str(self.dirs['fastqc']),
            "-t", str(self.config.threads),
            sample.fastq_r1
        ]
        if sample.fastq_r2:
            cmd.append(sample.fastq_r2)

        return self._execute(cmd, f"FastQC for {sample.sample_id}")

    def _trim_reads(self, sample: SampleInfo) -> Tuple[str, Optional[str]]:
        """Trim adapters and low-quality bases."""

        trimmed_r1 = str(self.dirs['trimmed'] / f"{sample.sample_id}_R1.trimmed.fastq.gz")
        trimmed_r2 = None

        if self.config.paired_end and sample.fastq_r2:
            trimmed_r2 = str(self.dirs['trimmed'] / f"{sample.sample_id}_R2.trimmed.fastq.gz")

            cmd = [
                "trim_galore",
                "--paired",
                "--quality", str(self.config.quality_threshold),
                "--length", str(self.config.min_length),
                "--cores", str(self.config.threads),
                "-o", str(self.dirs['trimmed']),
                sample.fastq_r1,
                sample.fastq_r2
            ]
        else:
            cmd = [
                "trim_galore",
                "--quality", str(self.config.quality_threshold),
                "--length", str(self.config.min_length),
                "--cores", str(self.config.threads),
                "-o", str(self.dirs['trimmed']),
                sample.fastq_r1
            ]

        self._execute(cmd, f"Trimming {sample.sample_id}")
        return trimmed_r1, trimmed_r2

    def _align_reads(self, sample: SampleInfo, r1: str, r2: str = None) -> str:
        """Align reads to reference genome."""

        output_bam = str(self.dirs['aligned'] / f"{sample.sample_id}.bam")
        genome_idx = self.GENOME_INDICES[self.config.genome][self.config.aligner]

        if self.config.aligner == "bowtie2":
            if r2:
                cmd = [
                    "bowtie2",
                    "-x", genome_idx,
                    "-1", r1,
                    "-2", r2,
                    "-p", str(self.config.threads),
                    "--maxins", str(self.config.max_fragment_size),
                    "--no-mixed",
                    "--no-discordant"
                ]
            else:
                cmd = [
                    "bowtie2",
                    "-x", genome_idx,
                    "-U", r1,
                    "-p", str(self.config.threads)
                ]

            # Pipe to samtools for BAM conversion and sorting
            full_cmd = f"{' '.join(cmd)} | samtools view -bS - | samtools sort -@ {self.config.threads} -o {output_bam}"
            self._execute_shell(full_cmd, f"Aligning {sample.sample_id}")

        elif self.config.aligner == "bwa":
            if r2:
                cmd = f"bwa mem -t {self.config.threads} {genome_idx} {r1} {r2}"
            else:
                cmd = f"bwa mem -t {self.config.threads} {genome_idx} {r1}"

            full_cmd = f"{cmd} | samtools view -bS - | samtools sort -@ {self.config.threads} -o {output_bam}"
            self._execute_shell(full_cmd, f"Aligning {sample.sample_id}")

        # Index BAM
        self._execute(["samtools", "index", output_bam], "Indexing BAM")

        return output_bam

    def _filter_bam(self, sample: SampleInfo, bam_file: str) -> str:
        """Filter BAM: remove duplicates, low MAPQ, blacklist regions."""

        filtered_bam = str(self.dirs['filtered'] / f"{sample.sample_id}.filtered.bam")

        # Remove duplicates with Picard
        dedup_bam = bam_file.replace(".bam", ".dedup.bam")
        if self.config.remove_duplicates:
            cmd = [
                "picard", "MarkDuplicates",
                f"I={bam_file}",
                f"O={dedup_bam}",
                f"M={bam_file.replace('.bam', '.dup_metrics.txt')}",
                "REMOVE_DUPLICATES=true"
            ]
            self._execute(cmd, f"Removing duplicates for {sample.sample_id}")
        else:
            dedup_bam = bam_file

        # Filter by MAPQ
        mapq_bam = dedup_bam.replace(".bam", ".mapq.bam")
        cmd = [
            "samtools", "view",
            "-b", "-q", str(self.config.min_mapq),
            "-o", mapq_bam,
            dedup_bam
        ]
        self._execute(cmd, f"Filtering MAPQ for {sample.sample_id}")

        # Remove blacklist regions
        if self.config.remove_blacklist:
            blacklist = self.GENOME_INDICES[self.config.genome].get('blacklist')
            if blacklist and os.path.exists(blacklist):
                cmd = [
                    "bedtools", "intersect",
                    "-v", "-abam", mapq_bam,
                    "-b", blacklist,
                    "-wa"
                ]
                with open(filtered_bam, 'wb') as f:
                    subprocess.run(cmd, stdout=f, check=True)
            else:
                filtered_bam = mapq_bam
        else:
            filtered_bam = mapq_bam

        # Index final BAM
        self._execute(["samtools", "index", filtered_bam], "Indexing filtered BAM")

        return filtered_bam

    def _generate_bigwig(self, sample: SampleInfo, bam_file: str) -> str:
        """Generate BigWig signal track."""

        bigwig_file = str(self.dirs['bigwig'] / f"{sample.sample_id}.bw")

        cmd = [
            "bamCoverage",
            "-b", bam_file,
            "-o", bigwig_file,
            "--normalizeUsing", "RPKM",
            "--binSize", "10",
            "-p", str(self.config.threads)
        ]

        self._execute(cmd, f"Generating BigWig for {sample.sample_id}")
        return bigwig_file

    def _call_peaks(self, sample: SampleInfo, bam_file: str, control_bam: str = None) -> str:
        """Call peaks using MACS2 or SICER."""

        # Determine peak type
        peak_type = self.MARK_PEAK_TYPES.get(
            sample.mark_or_tf,
            PeakType.NARROW if self.config.assay_type == AssayType.CHIP_SEQ_TF else self.config.peak_type
        )

        output_prefix = str(self.dirs['peaks'] / sample.sample_id)

        if self.config.peak_caller == "macs2":
            cmd = [
                "macs2", "callpeak",
                "-t", bam_file,
                "-n", sample.sample_id,
                "--outdir", str(self.dirs['peaks']),
                "-g", "hs" if "hg" in self.config.genome else "mm",
                "-q", str(self.config.qvalue_threshold)
            ]

            if control_bam:
                cmd.extend(["-c", control_bam])

            if peak_type == PeakType.BROAD:
                cmd.extend(["--broad", "--broad-cutoff", str(self.config.qvalue_threshold)])

            if self.config.paired_end:
                cmd.extend(["-f", "BAMPE"])

            self._execute(cmd, f"Calling peaks for {sample.sample_id}")

            # Return peak file path
            if peak_type == PeakType.BROAD:
                return f"{output_prefix}_peaks.broadPeak"
            else:
                return f"{output_prefix}_peaks.narrowPeak"

        return ""

    def _collect_stats(self, sample: SampleInfo, output: PipelineOutput) -> Dict:
        """Collect QC statistics."""

        stats = {
            'sample_id': sample.sample_id,
            'total_reads': 0,
            'mapped_reads': 0,
            'mapping_rate': 0.0,
            'duplicates': 0,
            'duplicate_rate': 0.0,
            'final_reads': 0,
            'peaks_called': 0,
            'frip': 0.0
        }

        # Get stats from samtools flagstat
        if output.bam_filtered and os.path.exists(output.bam_filtered):
            result = subprocess.run(
                ["samtools", "flagstat", output.bam_filtered],
                capture_output=True, text=True
            )
            # Parse flagstat output
            lines = result.stdout.split('\n')
            for line in lines:
                if 'total' in line and 'QC-passed' in line:
                    stats['final_reads'] = int(line.split()[0])

        # Count peaks
        if output.peaks_file and os.path.exists(output.peaks_file):
            with open(output.peaks_file) as f:
                stats['peaks_called'] = sum(1 for _ in f)

        return stats

    def _execute(self, cmd: List[str], description: str):
        """Execute a command."""
        print(f"Running: {description}")
        print(f"Command: {' '.join(cmd)}")
        # In production, uncomment:
        # subprocess.run(cmd, check=True)
        return True

    def _execute_shell(self, cmd: str, description: str):
        """Execute a shell command."""
        print(f"Running: {description}")
        print(f"Command: {cmd}")
        # In production, uncomment:
        # subprocess.run(cmd, shell=True, check=True)
        return True


def get_recommended_config(assay_type: str, mark_or_tf: str) -> PipelineConfig:
    """Get recommended pipeline configuration for a given assay/mark."""

    config = PipelineConfig(assay_type=AssayType(assay_type))

    # Set peak type based on mark
    if mark_or_tf in FASTQPipeline.MARK_PEAK_TYPES:
        config.peak_type = FASTQPipeline.MARK_PEAK_TYPES[mark_or_tf]

    # ATAC-seq specific settings
    if assay_type == "atac":
        config.remove_duplicates = True
        config.peak_caller = "macs2"
        config.peak_type = PeakType.NARROW

    # CUT&RUN/CUT&Tag settings
    if assay_type in ["cut_run", "cut_tag"]:
        config.peak_caller = "seacr"
        config.remove_duplicates = False  # Lower duplication expected

    return config
