#!/usr/bin/env python3
"""Extract ring breakpoint reads from GRCh38 CRAMs using T2T cluster results.

For each breakpoint cluster identified by the T2T pipeline, this script finds
the same reads (by read_id) in GRCh38 CRAMs and reports where R1 and its mate
map in GRCh38 coordinates. This helps define which discordant reads to look
for in GRCh38 mode.

Usage:
    uv run python scripts/extract_grch38_reads.py \
        --t2t-outputs examples/outputs-t2t/ \
        --cram-dir /path/to/grch38/crams/ \
        --ref /path/to/grch38.fa \
        --output grch38_breakpoint_reads.tsv
"""

from __future__ import annotations

import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd
import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def _resolve_path(path: Path) -> Path:
    """Resolve a path, handling broken symlinks on mounted volumes.

    When a network volume is mounted (e.g. /Volumes/share/), symlinks inside
    it may point to absolute server paths (e.g. /pasteur/helix/.../share/...).
    These are broken locally. This function detects broken symlinks, reads the
    target, and remaps by finding a common directory name between the target
    and the parent mount path.
    """
    if path.exists():
        return path
    if not path.is_symlink():
        return path

    target = os.readlink(path)
    # Walk up the parent directories to find a mount point name that
    # also appears in the symlink target path
    for parent in [path.parent, *path.parent.parents]:
        mount_name = parent.name
        if not mount_name:
            continue
        marker = f"/{mount_name}/"
        if marker in target:
            suffix = target[target.index(marker) + len(marker) :]
            remapped = parent / suffix
            if remapped.exists():
                logger.info("Resolved broken symlink: %s -> %s", path.name, remapped)
                return remapped
    return path


def _parse_clusters(summary_path: Path) -> list[dict]:
    """Parse summary.tsv → list of cluster dicts with read_positions."""
    df = pd.read_csv(summary_path, sep="\t")
    clusters = []
    for _, row in df.iterrows():
        positions = [int(p) for p in str(row["read_positions"]).split(",")]
        clusters.append({
            "sample": row["sample"],
            "cluster": row["cluster"],
            "chrom": row["chrom"],
            "breakpoint": row["breakpoint"],
            "positions": positions,
        })
    return clusters


def _match_reads_to_clusters(
    clusters: list[dict],
    supporting_reads_path: Path,
) -> list[dict]:
    """Match read_ids from supporting_reads.tsv to cluster positions.

    Returns a flat list of dicts, one per read, with cluster info attached.
    """
    reads_df = pd.read_csv(supporting_reads_path, sep="\t")

    # Build pos → [read_id, ...] multimap (preserving order for duplicate positions)
    pos_to_reads: dict[int, list[dict]] = defaultdict(list)
    for _, row in reads_df.iterrows():
        pos_to_reads[int(row["pos"])].append({
            "read_id": row["read_id"],
            "t2t_pos": int(row["pos"]),
            "t2t_mate_chrom": row["mate_chrom"],
            "t2t_mate_pos": int(row["mate_pos"]),
        })

    matched = []
    for cluster in clusters:
        # Track consumed indices per position to handle duplicate positions
        consumed: dict[int, int] = defaultdict(int)
        for pos in cluster["positions"]:
            candidates = pos_to_reads.get(pos, [])
            idx = consumed[pos]
            if idx < len(candidates):
                read_info = candidates[idx]
                consumed[pos] += 1
                matched.append({
                    "sample": cluster["sample"],
                    "cluster": cluster["cluster"],
                    "breakpoint": cluster["breakpoint"],
                    **read_info,
                })
            else:
                logger.warning(
                    "No read found for %s %s pos=%d (index=%d)",
                    cluster["sample"], cluster["cluster"], pos, idx,
                )
                matched.append({
                    "sample": cluster["sample"],
                    "cluster": cluster["cluster"],
                    "breakpoint": cluster["breakpoint"],
                    "read_id": "NOT_FOUND",
                    "t2t_pos": pos,
                    "t2t_mate_chrom": "NA",
                    "t2t_mate_pos": 0,
                })

    return matched


def _extract_grch38_info(
    read_ids: set[str],
    cram_path: Path,
    ref_path: Path,
    chrom: str,
    fetch_start: int,
    fetch_end: int,
) -> dict[str, dict]:
    """Fetch reads from GRCh38 CRAM and return mapping info keyed by read_id."""
    grch38_info: dict[str, dict] = {}

    with pysam.AlignmentFile(str(cram_path), reference_filename=str(ref_path)) as bam:
        # Verify chromosome exists in the CRAM
        refs = bam.references
        target_chrom = chrom
        if chrom not in refs:
            # Try without 'chr' prefix or with it
            alt = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"
            if alt in refs:
                target_chrom = alt
            else:
                logger.warning("Chromosome %s not found in %s", chrom, cram_path)
                return grch38_info

        start = max(0, fetch_start)
        logger.info("Fetching %s:%d-%d from %s", target_chrom, start, fetch_end, cram_path.name)

        for read in bam.fetch(target_chrom, start, fetch_end):
            if read.query_name in read_ids and read.query_name not in grch38_info:
                mate_chrom = read.next_reference_name if read.next_reference_name else "*"
                grch38_info[read.query_name] = {
                    "grch38_chrom": read.reference_name,
                    "grch38_pos": read.reference_start,
                    "grch38_mapq": read.mapping_quality,
                    "grch38_cigar": read.cigarstring,
                    "grch38_mate_chrom": mate_chrom,
                    "grch38_mate_pos": read.next_reference_start,
                    "grch38_mate_unmapped": read.mate_is_unmapped,
                    "grch38_is_proper_pair": read.is_proper_pair,
                    "grch38_flags": read.flag,
                }

    logger.info("Found %d/%d reads in GRCh38 CRAM", len(grch38_info), len(read_ids))
    return grch38_info


NA_GRCH38 = {
    "grch38_chrom": "NA",
    "grch38_pos": "NA",
    "grch38_mapq": "NA",
    "grch38_cigar": "NA",
    "grch38_mate_chrom": "NA",
    "grch38_mate_pos": "NA",
    "grch38_mate_unmapped": "NA",
    "grch38_is_proper_pair": "NA",
    "grch38_flags": "NA",
}

OUTPUT_COLUMNS = [
    "sample", "cluster", "breakpoint", "read_id",
    "t2t_pos", "t2t_mate_chrom", "t2t_mate_pos",
    "grch38_chrom", "grch38_pos", "grch38_mapq", "grch38_cigar",
    "grch38_mate_chrom", "grch38_mate_pos", "grch38_mate_unmapped",
    "grch38_is_proper_pair", "grch38_flags",
]


@click.command()
@click.option(
    "--t2t-outputs", required=True, type=click.Path(path_type=Path),
    help="Directory with T2T output files (*.summary.tsv + *.supporting_reads.tsv)",
)
@click.option(
    "--cram-dir", required=True, type=click.Path(path_type=Path),
    help="Directory with GRCh38 CRAM files (<sample>.cram)",
)
@click.option(
    "--ref", required=True, type=click.Path(path_type=Path),
    help="GRCh38 reference FASTA",
)
@click.option(
    "--output", "-o", default=None, type=click.Path(path_type=Path),
    help="Output TSV path (default: stdout)",
)
@click.option(
    "--buffer", default=5_000_000, type=int,
    help="Position buffer (bp) for GRCh38 fetch region around cluster",
)
@click.option(
    "--chrom", default="chr22", help="Target chromosome",
)
def main(
    t2t_outputs: Path,
    cram_dir: Path,
    ref: Path,
    output: Path | None,
    buffer: int,
    chrom: str,
) -> None:
    """Extract ring breakpoint reads from GRCh38 CRAMs.

    For each breakpoint cluster identified by the T2T pipeline, finds the
    same reads (by read_id) in GRCh38 CRAMs and reports R1 + mate mapping.
    """
    # Resolve paths (handles broken symlinks on mounted volumes)
    t2t_outputs = _resolve_path(t2t_outputs)
    cram_dir = _resolve_path(cram_dir)
    ref = _resolve_path(ref)

    for label, path in [("--t2t-outputs", t2t_outputs), ("--cram-dir", cram_dir), ("--ref", ref)]:
        if not path.exists():
            logger.error("%s path does not exist: %s", label, path)
            sys.exit(1)

    # Discover samples from summary files
    summary_files = sorted(t2t_outputs.glob("*.summary.tsv"))
    if not summary_files:
        logger.error("No *.summary.tsv files found in %s", t2t_outputs)
        sys.exit(1)

    all_rows: list[dict] = []

    for summary_path in summary_files:
        sample = summary_path.stem.replace(".summary", "")
        supporting_path = t2t_outputs / f"{sample}.supporting_reads.tsv"

        if not supporting_path.exists():
            logger.warning("No supporting_reads.tsv for %s, skipping", sample)
            continue

        # Find GRCh38 CRAM (resolve broken symlinks on mounted volumes)
        cram_path = _resolve_path(cram_dir / f"{sample}.cram")
        if not cram_path.exists():
            logger.warning("No GRCh38 CRAM for %s at %s, skipping", sample, cram_path)
            continue

        logger.info("Processing %s", sample)

        # Step 1-3: Parse clusters and match read_ids
        clusters = _parse_clusters(summary_path)
        if not clusters:
            logger.info("No clusters for %s", sample)
            continue

        matched_reads = _match_reads_to_clusters(clusters, supporting_path)

        # Step 4-5: Determine fetch region and extract from GRCh38
        read_ids = {r["read_id"] for r in matched_reads if r["read_id"] != "NOT_FOUND"}
        all_positions = [r["t2t_pos"] for r in matched_reads]
        fetch_start = min(all_positions) - buffer
        fetch_end = max(all_positions) + buffer

        grch38_info = _extract_grch38_info(
            read_ids, cram_path, ref, chrom, fetch_start, fetch_end,
        )

        # Step 6: Merge T2T + GRCh38 info
        for read in matched_reads:
            rid = read["read_id"]
            g38 = grch38_info.get(rid, NA_GRCH38)
            all_rows.append({**read, **g38})

    # Output
    if not all_rows:
        logger.warning("No reads found across any sample")
        sys.exit(0)

    result = pd.DataFrame(all_rows, columns=OUTPUT_COLUMNS)

    if output:
        result.to_csv(output, sep="\t", index=False)
        logger.info("Written %d rows to %s", len(result), output)
    else:
        result.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
