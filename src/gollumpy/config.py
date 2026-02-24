"""Configuration and chromosome region definitions."""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Literal

import tomllib

ACROCENTRIC_CHROMS = ("chr13", "chr14", "chr15", "chr21", "chr22")


@dataclass(frozen=True)
class ChromRegion:
    """SAAC (Short Arm Acrocentric Chromosome) region boundary."""

    chrom: str
    saac_start: int
    saac_end: int


@dataclass
class FilterParams:
    """Read and alignment filtering thresholds."""

    min_mapq_r1: int = 40
    samtools_mapq: int = 40
    min_alignment_score: int = 80
    max_divergence: float = 0.06
    centromere_buffer: int = 5_000_000


@dataclass
class ClusterParams:
    """HDBSCAN clustering parameters."""

    min_cluster_size: int = 3
    min_samples: int = 5
    cluster_selection_epsilon: float = 100.0
    allow_single_cluster: bool = True
    min_supporting_reads: int = 5


@dataclass
class GollumConfig:
    """Main configuration for a gollum run."""

    target_chrom: str
    mode: Literal["t2t", "grch38"]
    fasta_t2t: Path
    input_file: Path
    output_dir: Path
    sample_name: str
    fasta_grch38: Path | None = None
    blacklist_bed: Path | None = None
    filter_params: FilterParams = field(default_factory=FilterParams)
    cluster_params: ClusterParams = field(default_factory=ClusterParams)

    @property
    def reference_fasta(self) -> Path:
        """Reference FASTA matching the input CRAM/BAM alignment.

        In T2T mode, returns fasta_t2t. In GRCh38 mode, returns fasta_grch38.
        """
        if self.mode == "grch38":
            if self.fasta_grch38 is None:
                msg = "GRCh38 reference is required in grch38 mode"
                raise ValueError(msg)
            return self.fasta_grch38
        return self.fasta_t2t

    def __post_init__(self) -> None:
        if self.target_chrom not in ACROCENTRIC_CHROMS:
            msg = f"Invalid chromosome '{self.target_chrom}'. Must be one of {ACROCENTRIC_CHROMS}"
            raise ValueError(msg)
        if not self.input_file.exists():
            msg = f"Input file not found: {self.input_file}"
            raise FileNotFoundError(msg)
        if not self.fasta_t2t.exists():
            msg = f"T2T reference not found: {self.fasta_t2t}"
            raise FileNotFoundError(msg)
        if self.mode == "grch38":
            if self.fasta_grch38 is None:
                msg = "GRCh38 reference (--ref) is required in grch38 mode"
                raise ValueError(msg)
            if not self.fasta_grch38.exists():
                msg = f"GRCh38 reference not found: {self.fasta_grch38}"
                raise FileNotFoundError(msg)
        if self.blacklist_bed is not None and not self.blacklist_bed.exists():
            msg = f"Blacklist BED not found: {self.blacklist_bed}"
            raise FileNotFoundError(msg)
        self.output_dir.mkdir(parents=True, exist_ok=True)


def load_acro_regions() -> dict[str, ChromRegion]:
    """Load acrocentric chromosome SAAC regions from bundled TOML resource."""
    resource_file = resources.files("gollumpy") / "resources" / "acro_regions.toml"
    data = tomllib.loads(resource_file.read_text(encoding="utf-8"))
    return {
        chrom: ChromRegion(chrom=chrom, saac_start=vals["saac_start"], saac_end=vals["saac_end"])
        for chrom, vals in data.items()
    }


def load_blacklist(blacklist_bed: Path) -> list[tuple[str, int, int]]:
    """Load blacklist BED file as list of (chrom, start, end) tuples."""
    regions: list[tuple[str, int, int]] = []
    with open(blacklist_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def is_blacklisted(chrom: str, pos: int, blacklist: list[tuple[str, int, int]]) -> bool:
    """Check if a position falls within any blacklisted region."""
    return any(bl_chrom == chrom and bl_start <= pos < bl_end for bl_chrom, bl_start, bl_end in blacklist)
