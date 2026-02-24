"""Pipeline orchestrator: extract → align → cluster → report."""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

from gollumpy.align import align_mates, build_aligner, compute_acro_specificity, compute_mate_concordance
from gollumpy.cluster import cluster_breakpoints, pre_cluster_reads
from gollumpy.extract import (
    extract_discordant_reads_grch38,
    extract_discordant_reads_t2t,
    extract_mate_sequences,
)
from gollumpy.report import generate_report

if TYPE_CHECKING:
    from gollumpy.config import GollumConfig
    from gollumpy.models import Breakpoint

logger = logging.getLogger(__name__)


def compute_ring_score(
    supporting_reads: int,
    confidence: float,
    acro_specificity: float | None,
    mate_concordance: float,
    span: int = 0,
) -> float:
    """Compute composite ring score for a breakpoint (0–10).

    Higher values indicate a more likely real ring chromosome breakpoint.

    Components (each normalized to ~0–1):
      - mate_concordance (weight 0.25): fraction of mates converging on one SAAC chrom
      - read_signal (weight 0.25): min(reads / 15, 1.0) — linear, caps at 15 reads
      - confidence (weight 0.10): HDBSCAN mean probability (already 0–1)
      - specificity (weight 0.15): min(acro_specificity, 20) / 20, capped at 1.0
      - span_tightness (weight 0.25): penalizes wide clusters via log10 scale
        (~1.0 at 1bp, ~0.57 at 1kb, ~0.14 at 1Mb, 0.0 at ≥10Mb)

    v3 recalibration: read_signal uses linear scaling (was log2) to better
    separate real clusters (10–20 reads) from noise clusters (3–5 reads).
    Confidence weight reduced from 0.20 to 0.10 since HDBSCAN confidence
    favours small tight clusters regardless of biological significance.
    """
    read_signal = min(supporting_reads / 15.0, 1.0)

    if acro_specificity is None:
        spec_value = 0.0
    elif acro_specificity == float("inf"):
        spec_value = 1.0
    else:
        spec_value = min(acro_specificity, 20.0) / 20.0

    span_tightness = max(1.0 - math.log10(max(span, 1)) / 7.0, 0.0)

    score = (
        0.25 * mate_concordance
        + 0.25 * read_signal
        + 0.10 * confidence
        + 0.15 * spec_value
        + 0.25 * span_tightness
    ) * 10.0

    return round(score, 3)


def run_pipeline(config: GollumConfig) -> list[Breakpoint]:
    """Run the full gollum ring detection pipeline.

    Steps:
    1. Extract discordant reads (T2T or GRCh38 mode)
    2. Extract mate sequences
    3. Align mates to T2T reference with minimap2
    4. Cluster breakpoint positions with HDBSCAN
    5. Filter, enrich, score, and report
    """
    # Configure logging
    log_path = config.output_dir / f"{config.sample_name}.log"
    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] %(name)s %(levelname)s - %(message)s", datefmt="%H:%M:%S")
    file_handler.setFormatter(formatter)
    logging.getLogger("gollumpy").addHandler(file_handler)
    logging.getLogger("gollumpy").setLevel(logging.DEBUG)

    logger.info("Gollum v2 pipeline starting")
    logger.info("Mode: %s, Target: %s, Sample: %s", config.mode, config.target_chrom, config.sample_name)

    # Step 1: Extract discordant reads
    extract_fn = extract_discordant_reads_t2t if config.mode == "t2t" else extract_discordant_reads_grch38
    reads = extract_fn(config)

    if reads.empty:
        logger.info("No discordant reads found — no ring detected")
        generate_report([], reads, config)
        return []

    # Step 1b: Pre-cluster by R1 position to remove scattered noise
    # (v1's two-pass architecture: positional pre-filter before alignment)
    before_precluster = len(reads)
    reads = pre_cluster_reads(reads)
    logger.info("Pre-clustering: %d → %d reads", before_precluster, len(reads))

    if reads.empty:
        logger.info("No reads survived pre-clustering — no ring detected")
        generate_report([], reads, config)
        return []

    # Step 2: Extract mate sequences
    reads = extract_mate_sequences(config, reads)

    if reads.empty:
        logger.info("No mate sequences extracted — no ring detected")
        generate_report([], reads, config)
        return []

    # Step 3: Align mates with minimap2
    aligned = align_mates(reads, config)

    if aligned.empty:
        logger.info("No mates aligned to SAAC regions — no ring detected")
        generate_report([], aligned, config)
        return []

    # Step 3b: Target-chrom PHR filter — keep only reads with ANY SAAC hit
    # on the target chromosome (v1's isPHR semantics: any hit, not just best)
    before_phr = len(aligned)
    aligned = aligned[aligned["has_target_saac_hit"]]
    logger.info("After target-chrom filter: %d/%d reads", len(aligned), before_phr)

    if aligned.empty:
        logger.info("No mates aligned to target SAAC — no ring detected")
        generate_report([], aligned, config)
        return []

    # Step 4: Cluster breakpoints
    breakpoints, labeled_reads = cluster_breakpoints(aligned, config.cluster_params)

    if not breakpoints:
        logger.info("No ring detected after clustering")
        generate_report([], aligned, config)
    else:
        # Step 4a: Filter clusters by min_supporting_reads
        clustered_reads = labeled_reads[labeled_reads["cluster"] != -1]
        cluster_ids = sorted(clustered_reads["cluster"].unique())

        min_reads = config.cluster_params.min_supporting_reads
        keep_indices = [i for i, bp in enumerate(breakpoints) if bp.supporting_reads >= min_reads]

        if len(keep_indices) < len(breakpoints):
            logger.info(
                "Filtered %d/%d clusters with < %d supporting reads",
                len(breakpoints) - len(keep_indices),
                len(breakpoints),
                min_reads,
            )

        breakpoints = [breakpoints[i] for i in keep_indices]
        keep_cluster_ids = [cluster_ids[i] for i in keep_indices]

        if not breakpoints:
            logger.info("No ring detected after min_supporting_reads filter")
            generate_report([], aligned, config)
        else:
            # Step 4b: Enrich surviving clusters with per-cluster metrics
            aligner = build_aligner(config)

            for bp, cid in zip(breakpoints, keep_cluster_ids, strict=True):
                group = clustered_reads[clustered_reads["cluster"] == cid]
                bp.acro_specificity = compute_acro_specificity(group, config, aligner=aligner)
                concordance, dominant_chrom = compute_mate_concordance(group)
                bp.mate_concordance = concordance
                bp.dominant_saac_chrom = dominant_chrom

            # Step 4c: Compute ring_score
            for bp in breakpoints:
                bp.ring_score = compute_ring_score(
                    bp.supporting_reads,
                    bp.confidence,
                    bp.acro_specificity,
                    bp.mate_concordance or 0.0,
                    span=bp.pos_max - bp.pos_min,
                )

            # Sort by ring_score descending (best candidate first)
            breakpoints.sort(key=lambda bp: (bp.ring_score or 0.0), reverse=True)
            logger.info("Detected %d breakpoint(s)", len(breakpoints))

            # Step 5: Report
            generate_report(breakpoints, aligned, config)

    # Clean up file handler
    logging.getLogger("gollumpy").removeHandler(file_handler)
    file_handler.close()

    return breakpoints
