"""Output generation for gollum results."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from gollumpy.config import GollumConfig
    from gollumpy.models import Breakpoint

logger = logging.getLogger(__name__)


def _fmt(value: float | None, decimals: int = 3) -> str:
    """Format a nullable float for output."""
    if value is None:
        return "NA"
    return f"{value:.{decimals}f}"


def generate_report(
    breakpoints: list[Breakpoint],
    reads_df: pd.DataFrame,
    config: GollumConfig,
) -> None:
    """Write gollum results to TSV files and log summary."""
    output_dir = config.output_dir
    sample = config.sample_name

    # Write summary TSV
    summary_path = output_dir / f"{sample}.summary.tsv"
    if breakpoints:
        summary_rows = []
        for i, bp in enumerate(breakpoints):
            summary_rows.append({
                "sample": sample,
                "cluster": f"cluster_{i}",
                "chrom": bp.chrom,
                "breakpoint": f"{bp.chrom}:{bp.position}",
                "span": f"{bp.pos_min}-{bp.pos_max}",
                "supporting_reads": bp.supporting_reads,
                "confidence": _fmt(bp.confidence),
                "acro_specificity": _fmt(bp.acro_specificity),
                "mate_concordance": _fmt(bp.mate_concordance),
                "dominant_saac_chrom": bp.dominant_saac_chrom or "NA",
                "ring_score": _fmt(bp.ring_score),
            })
        df_summary = pd.DataFrame(summary_rows)
        df_summary.to_csv(summary_path, sep="\t", index=False)
        logger.info("Summary written to %s", summary_path)
    else:
        summary_path.write_text("no ring detected.\n")
        logger.info("No ring detected")

    # Write supporting reads TSV
    if not reads_df.empty:
        reads_path = output_dir / f"{sample}.supporting_reads.tsv"
        reads_df.to_csv(reads_path, sep="\t", index=False)
        logger.info("Supporting reads written to %s", reads_path)

    # Print summary to stdout
    if breakpoints:
        header = (
            "sample\tcluster\tbreakpoint\tsupporting_reads\tconfidence"
            "\tacro_specificity\tmate_concordance\tring_score"
        )
        print(header)
        for i, bp in enumerate(breakpoints):
            print(
                f"{sample}\tcluster_{i}\t{bp.chrom}:{bp.position}\t{bp.supporting_reads}"
                f"\t{_fmt(bp.confidence)}\t{_fmt(bp.acro_specificity)}"
                f"\t{_fmt(bp.mate_concordance)}\t{_fmt(bp.ring_score)}"
            )
    else:
        print("no ring detected.")
