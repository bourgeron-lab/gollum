"""R2 mate realignment to T2T reference using mappy (minimap2 Python bindings)."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import mappy
import pandas as pd

from gollumpy.config import ACROCENTRIC_CHROMS, load_acro_regions

if TYPE_CHECKING:
    from gollumpy.config import GollumConfig

logger = logging.getLogger(__name__)


def align_mates(reads_df: pd.DataFrame, config: GollumConfig) -> pd.DataFrame:
    """Align mate sequences to T2T reference using minimap2 via mappy.

    For each mate sequence, finds all alignments to the T2T reference.
    Filters alignments by quality and checks if they fall within SAAC regions.
    Returns the input DataFrame filtered to reads with at least one SAAC hit,
    plus alignment metadata columns.
    """
    if reads_df.empty or "mate_sequence" not in reads_df.columns:
        return pd.DataFrame()

    acro_regions = load_acro_regions()

    # Build minimap2 index
    logger.info("Building minimap2 index for %s", config.fasta_t2t)
    aligner = mappy.Aligner(str(config.fasta_t2t), preset="sr", best_n=5)
    if not aligner:
        msg = f"Failed to build minimap2 index for {config.fasta_t2t}"
        raise RuntimeError(msg)

    # Collect all alignments
    alignment_records: list[dict] = []

    for _, row in reads_df.iterrows():
        read_id = row["read_id"]
        seq = row["mate_sequence"]

        if not seq or not isinstance(seq, str):
            continue

        for hit in aligner.map(seq):
            divergence = hit.NM / hit.blen if hit.blen > 0 else 1.0
            alignment_records.append({
                "read_id": read_id,
                "align_chrom": hit.ctg,
                "align_start": hit.r_st,
                "align_end": hit.r_en,
                "mapq": hit.mapq,
                "NM": hit.NM,
                "blen": hit.blen,
                "mlen": hit.mlen,
                "divergence": divergence,
                "is_primary": hit.is_primary,
            })

    if not alignment_records:
        logger.info("No alignments found for mate sequences")
        return pd.DataFrame()

    df_align = pd.DataFrame(alignment_records)
    logger.info("Total alignments: %d", len(df_align))

    # Filter on alignment quality
    # Do NOT filter on mapQ — multi-mapping is expected in SAAC regions
    df_align = df_align[df_align["mlen"] >= config.filter_params.min_alignment_score]
    df_align = df_align[df_align["divergence"] <= config.filter_params.max_divergence]
    logger.info("Alignments after quality filter: %d", len(df_align))

    if df_align.empty:
        logger.info("No alignments passed quality filters")
        return pd.DataFrame()

    # Check which alignments fall within SAAC regions
    def is_in_saac(row: pd.Series) -> bool:
        chrom = row["align_chrom"]
        if chrom not in acro_regions:
            return False
        region = acro_regions[chrom]
        # Alignment overlaps SAAC if it starts before centromere boundary
        return row["align_start"] < region.saac_end

    df_align["is_saac"] = df_align.apply(is_in_saac, axis=1)
    saac_hits = df_align[df_align["is_saac"]]
    logger.info("Alignments in SAAC regions: %d", len(saac_hits))

    if saac_hits.empty:
        logger.info("No alignments mapped to SAAC regions")
        return pd.DataFrame()

    # Get read IDs with at least one SAAC hit
    reads_with_saac = set(saac_hits["read_id"])

    # Compute per-read best alignment stats for SAAC hits
    saac_stats = (
        saac_hits.groupby("read_id")
        .agg(
            best_mlen=("mlen", "max"),
            best_divergence=("divergence", "min"),
            n_saac_hits=("read_id", "count"),
        )
        .reset_index()
    )

    # Filter original reads to those with SAAC hits and merge stats
    result = reads_df[reads_df["read_id"].isin(reads_with_saac)].merge(saac_stats, on="read_id")
    logger.info("Reads with SAAC-confirmed mates: %d", len(result))
    return result


def build_aligner(config: GollumConfig) -> mappy.Aligner:
    """Build a minimap2 aligner for the T2T reference."""
    aligner = mappy.Aligner(str(config.fasta_t2t), preset="sr", best_n=5)
    if not aligner:
        msg = f"Failed to build minimap2 index for {config.fasta_t2t}"
        raise RuntimeError(msg)
    return aligner


def compute_acro_specificity(
    reads_df: pd.DataFrame,
    config: GollumConfig,
    *,
    aligner: mappy.Aligner | None = None,
) -> float | None:
    """Compute acrocentric specificity ratio for aligned reads.

    Ratio = sum(mlen for acrocentric hits) / sum(mlen for non-acrocentric hits).
    Higher values indicate reads preferentially align to acrocentric chromosomes.
    """
    if reads_df.empty or "mate_sequence" not in reads_df.columns:
        return None

    if aligner is None:
        aligner = mappy.Aligner(str(config.fasta_t2t), preset="sr", best_n=5)
        if not aligner:
            return None

    acro_mlen_sum = 0
    non_acro_mlen_sum = 0

    for _, row in reads_df.iterrows():
        seq = row["mate_sequence"]
        if not seq or not isinstance(seq, str):
            continue

        for hit in aligner.map(seq):
            if hit.ctg in ACROCENTRIC_CHROMS:
                acro_mlen_sum += hit.mlen
            else:
                non_acro_mlen_sum += hit.mlen

    if non_acro_mlen_sum == 0:
        return float("inf") if acro_mlen_sum > 0 else None

    return acro_mlen_sum / non_acro_mlen_sum
