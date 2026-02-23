"""Command-line interface for gollum."""

from __future__ import annotations

import functools
from pathlib import Path

import click

from gollumpy import __version__
from gollumpy.config import ACROCENTRIC_CHROMS, ClusterParams, FilterParams, GollumConfig
from gollumpy.pipeline import run_pipeline


def shared_options(func):  # noqa: ANN001, ANN201
    """Decorator that adds shared CLI options for both t2t and grch38 subcommands."""

    @click.argument("cram", type=click.Path(exists=True, path_type=Path))
    @click.option(
        "-f", "--fasta", required=True, type=click.Path(exists=True, path_type=Path), help="T2T-CHM13 reference FASTA"
    )
    @click.option("-o", "--output-dir", required=True, type=click.Path(path_type=Path), help="Output directory")
    @click.option(
        "--chr", "chrom", default="chr22", type=click.Choice(list(ACROCENTRIC_CHROMS)), help="Target chromosome"
    )
    @click.option("-s", "--sample-name", default=None, help="Sample name (default: from filename)")
    @click.option(
        "--blacklist-bed", default=None, type=click.Path(exists=True, path_type=Path), help="Blacklisted regions BED"
    )
    @click.option("--min-alignment-score", default=80, type=int, help="Minimum alignment matching length (mlen)")
    @click.option("--max-divergence", default=0.06, type=float, help="Maximum sequence divergence (NM/blen)")
    @click.option("--min-cluster-size", default=3, type=int, help="HDBSCAN min_cluster_size")
    @click.option("--min-samples", default=5, type=int, help="HDBSCAN min_samples")
    @click.option("--cluster-epsilon", default=100.0, type=float, help="HDBSCAN cluster_selection_epsilon (bp)")
    @click.option("--allow-single/--no-allow-single", default=True, help="Allow single-cluster detection")
    @click.option("--min-supporting-reads", default=5, type=int, help="Minimum reads per cluster to report")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        return func(*args, **kwargs)

    return wrapper


def _build_config(
    mode: str,
    cram: Path,
    fasta: Path,
    output_dir: Path,
    chrom: str,
    sample_name: str | None,
    blacklist_bed: Path | None,
    min_alignment_score: int,
    max_divergence: float,
    min_cluster_size: int,
    min_samples: int,
    cluster_epsilon: float,
    allow_single: bool,
    min_supporting_reads: int,
    fasta_grch38: Path | None = None,
) -> GollumConfig:
    """Build GollumConfig from CLI arguments."""
    if sample_name is None:
        sample_name = cram.stem

    return GollumConfig(
        target_chrom=chrom,
        mode=mode,  # type: ignore[arg-type]
        fasta_t2t=fasta,
        input_file=cram,
        output_dir=output_dir,
        sample_name=sample_name,
        fasta_grch38=fasta_grch38,
        blacklist_bed=blacklist_bed,
        filter_params=FilterParams(
            min_alignment_score=min_alignment_score,
            max_divergence=max_divergence,
        ),
        cluster_params=ClusterParams(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_epsilon=cluster_epsilon,
            allow_single_cluster=allow_single,
            min_supporting_reads=min_supporting_reads,
        ),
    )


@click.group()
@click.version_option(version=__version__)
def main() -> None:
    """Gollum — Detection of ring chromosomes on acrocentric chromosomes."""


@main.command()
@shared_options
def t2t(
    cram: Path,
    fasta: Path,
    output_dir: Path,
    chrom: str,
    sample_name: str | None,
    blacklist_bed: Path | None,
    min_alignment_score: int,
    max_divergence: float,
    min_cluster_size: int,
    min_samples: int,
    cluster_epsilon: float,
    allow_single: bool,
    min_supporting_reads: int,
) -> None:
    """Detect ring chromosomes from T2T-CHM13 aligned BAM/CRAM."""
    config = _build_config(
        mode="t2t",
        cram=cram,
        fasta=fasta,
        output_dir=output_dir,
        chrom=chrom,
        sample_name=sample_name,
        blacklist_bed=blacklist_bed,
        min_alignment_score=min_alignment_score,
        max_divergence=max_divergence,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_epsilon=cluster_epsilon,
        allow_single=allow_single,
        min_supporting_reads=min_supporting_reads,
    )
    run_pipeline(config)


@main.command()
@shared_options
@click.option(
    "--ref", "reference", required=True, type=click.Path(exists=True, path_type=Path),
    help="GRCh38 reference FASTA (for CRAM decoding)",
)
def grch38(
    cram: Path,
    fasta: Path,
    output_dir: Path,
    chrom: str,
    sample_name: str | None,
    blacklist_bed: Path | None,
    min_alignment_score: int,
    max_divergence: float,
    min_cluster_size: int,
    min_samples: int,
    cluster_epsilon: float,
    allow_single: bool,
    min_supporting_reads: int,
    reference: Path,
) -> None:
    """Detect ring chromosomes from GRCh38-aligned BAM/CRAM.

    Requires both a GRCh38 reference (--ref, for CRAM decoding) and a
    T2T-CHM13 reference (-f, for R2 mate realignment).
    """
    config = _build_config(
        mode="grch38",
        cram=cram,
        fasta=fasta,
        output_dir=output_dir,
        chrom=chrom,
        sample_name=sample_name,
        blacklist_bed=blacklist_bed,
        min_alignment_score=min_alignment_score,
        max_divergence=max_divergence,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_epsilon=cluster_epsilon,
        allow_single=allow_single,
        min_supporting_reads=min_supporting_reads,
        fasta_grch38=reference,
    )
    run_pipeline(config)
