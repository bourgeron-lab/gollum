"""Pipeline orchestrator: extract → align → cluster → report."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from gollumpy.align import align_mates, build_aligner, compute_acro_specificity
from gollumpy.cluster import cluster_breakpoints
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


def run_pipeline(config: GollumConfig) -> list[Breakpoint]:
    """Run the full gollum ring detection pipeline.

    Steps:
    1. Extract discordant reads (T2T or GRCh38 mode)
    2. Extract mate sequences
    3. Align mates to T2T reference with minimap2
    4. Cluster breakpoint positions with HDBSCAN
    5. Generate report
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

    # Step 4: Cluster breakpoints
    breakpoints, labeled_reads = cluster_breakpoints(aligned, config.cluster_params)

    if not breakpoints:
        logger.info("No ring detected after clustering")
        generate_report([], aligned, config)
    else:
        # Step 4b: Compute per-cluster acrocentric specificity
        aligner = build_aligner(config)
        clustered_reads = labeled_reads[labeled_reads["cluster"] != -1]

        for bp, (_cid, group) in zip(breakpoints, clustered_reads.groupby("cluster"), strict=True):
            bp.acro_specificity = compute_acro_specificity(group, config, aligner=aligner)

        breakpoints.sort(key=lambda bp: bp.position)
        logger.info("Detected %d breakpoint(s)", len(breakpoints))

        # Step 5: Report
        generate_report(breakpoints, aligned, config)

    # Clean up file handler
    logging.getLogger("gollumpy").removeHandler(file_handler)
    file_handler.close()

    return breakpoints
