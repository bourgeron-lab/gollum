"""Single-pass HDBSCAN clustering of breakpoint positions."""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import hdbscan
import numpy as np
import pandas as pd

from gollumpy.models import Breakpoint

if TYPE_CHECKING:
    from gollumpy.config import ClusterParams

logger = logging.getLogger(__name__)


def pre_cluster_reads(
    reads_df: pd.DataFrame,
    *,
    min_cluster_size: int = 2,
    min_samples: int = 2,
    cluster_selection_epsilon: float = 500.0,
) -> pd.DataFrame:
    """Pre-cluster discordant reads by R1 position to remove scattered noise.

    This replicates v1's first-pass clustering: a quick positional grouping
    that discards reads not near any positional hotspot (HDBSCAN noise label).
    Uses relaxed parameters to avoid discarding real signal.

    Returns a filtered DataFrame containing only reads assigned to a cluster.
    """
    if reads_df.empty or "pos" not in reads_df.columns:
        return reads_df

    positions = np.array(reads_df["pos"], dtype=float).reshape(-1, 1)

    if len(positions) < min_cluster_size:
        return reads_df

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        allow_single_cluster=True,
    )
    labels = clusterer.fit_predict(positions)

    n_clustered = int(np.sum(labels != -1))
    n_noise = int(np.sum(labels == -1))
    n_clusters = len(set(labels) - {-1})
    logger.info(
        "Pre-clustering: %d clusters, %d reads kept, %d noise discarded",
        n_clusters,
        n_clustered,
        n_noise,
    )

    if n_clustered == 0:
        logger.info("Pre-clustering removed all reads — returning original set")
        return reads_df

    return reads_df[labels != -1].reset_index(drop=True)


def _trim_cluster_outliers(reads_df: pd.DataFrame) -> pd.DataFrame:
    """Remove positional outliers from each cluster using IQR method.

    For each cluster, positions beyond Q1 − 1.5×IQR or Q3 + 1.5×IQR
    are relabeled as noise (cluster = −1).  Only applied to clusters
    with ≥ 5 reads to avoid trimming small but valid clusters.
    """
    reads_df = reads_df.copy()
    cluster_ids = reads_df.loc[reads_df["cluster"] != -1, "cluster"].unique()
    for cid in cluster_ids:
        mask = reads_df["cluster"] == cid
        if mask.sum() < 5:
            continue
        positions = reads_df.loc[mask, "pos"]
        q1, q3 = positions.quantile([0.25, 0.75])
        iqr = q3 - q1
        if iqr == 0:
            continue
        lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        outlier_mask = mask & (
            (reads_df["pos"] < lower) | (reads_df["pos"] > upper)
        )
        reads_df.loc[outlier_mask, "cluster"] = -1
    return reads_df


def _refine_wide_clusters(
    reads_df: pd.DataFrame,
    max_span: int,
    *,
    gap_threshold: int = 500,
    min_subcluster: int = 3,
) -> pd.DataFrame:
    """Extract tight sub-clusters from clusters exceeding *max_span* bp.

    HDBSCAN sometimes absorbs tight signal reads into wide noise clusters.
    Instead of discarding the entire cluster, we attempt gap-based splitting:

    1. Sort positions within the wide cluster.
    2. Find gaps > *gap_threshold* bp between consecutive reads.
    3. Split into sub-groups at those gaps.
    4. Keep sub-groups with ≥ *min_subcluster* reads AND span ≤ *max_span*.
    5. Assign new cluster IDs to surviving sub-groups; discard the rest as noise.

    Clusters already within *max_span* are left untouched.
    """
    reads_df = reads_df.copy()
    next_cluster_id = int(reads_df["cluster"].max()) + 1 if not reads_df.empty else 0

    for cid in reads_df.loc[reads_df["cluster"] != -1, "cluster"].unique():
        mask = reads_df["cluster"] == cid
        positions = reads_df.loc[mask, "pos"]
        span = int(positions.max() - positions.min())
        if span <= max_span:
            continue

        # Sort by position and find large gaps
        sorted_idx = positions.sort_values().index
        sorted_pos = positions.loc[sorted_idx].values
        gaps = np.diff(sorted_pos)
        split_points = np.where(gaps > gap_threshold)[0] + 1

        # Split into sub-groups
        sub_indices = np.split(sorted_idx.values, split_points)

        # Reclassify all reads in this cluster as noise first
        reads_df.loc[mask, "cluster"] = -1

        # Re-assign tight sub-groups that meet criteria
        for sub_idx in sub_indices:
            if len(sub_idx) < min_subcluster:
                continue
            sub_pos = reads_df.loc[sub_idx, "pos"]
            sub_span = int(sub_pos.max() - sub_pos.min())
            if sub_span > max_span:
                continue
            reads_df.loc[sub_idx, "cluster"] = next_cluster_id
            next_cluster_id += 1

    return reads_df


def _parse_clip_position(
    cigarstring: str,
    ref_start: int,
    ref_end: int,
    min_clip: int = 20,
) -> int | None:
    """Extract breakpoint position from CIGAR soft-clipping.

    Right-side clip (nM…kS): breakpoint at *reference_end*.
    Left-side clip  (kS…nM): breakpoint at *reference_start*.
    Returns ``None`` if no significant soft-clip (< *min_clip* bp).
    """
    if not cigarstring or not isinstance(cigarstring, str):
        return None
    # Right-side soft clip
    match = re.search(r"(\d+)S$", cigarstring)
    if match and int(match.group(1)) >= min_clip:
        return ref_end
    # Left-side soft clip
    match = re.match(r"^(\d+)S", cigarstring)
    if match and int(match.group(1)) >= min_clip:
        return ref_start
    return None


def _refine_breakpoint_position(group: pd.DataFrame) -> int:
    """Determine precise breakpoint from soft-clip consensus, else median.

    1. Parse CIGARs for significant soft-clips (≥ 20 bp).
    2. If ≥ 3 clips agree (IQR < 50 bp), use their median.
    3. Otherwise, fall back to median of all R1 positions.
    """
    clip_positions: list[int] = []

    if "cigarstring" in group.columns and "reference_end" in group.columns:
        for _, row in group.iterrows():
            clip_pos = _parse_clip_position(
                row.get("cigarstring", ""),
                int(row["pos"]),
                int(row.get("reference_end", row["pos"])),
            )
            if clip_pos is not None:
                clip_positions.append(clip_pos)

    if len(clip_positions) >= 3:
        clips = pd.Series(clip_positions)
        q1, q3 = clips.quantile([0.25, 0.75])
        if q3 - q1 < 50:  # clips agree within 50 bp
            return int(clips.median())

    # Fallback: use reference_end (read 3' end) when available — closer to
    # the ring junction than reference_start for right-aligned reads.
    if "reference_end" in group.columns:
        return int(group["reference_end"].median())

    return int(group["pos"].median())


def cluster_breakpoints(
    reads_df: pd.DataFrame,
    cluster_params: ClusterParams,
) -> tuple[list[Breakpoint], pd.DataFrame]:
    """Cluster R1 positions of filtered reads using HDBSCAN.

    Single-pass clustering on the R1 positions of reads whose R2 mate
    passed alignment filters. Returns a tuple of (breakpoints, labeled_df)
    where labeled_df has ``cluster`` and ``probability`` columns added.
    Breakpoints have ``acro_specificity=None``; the caller is responsible
    for enriching each breakpoint with its per-cluster specificity.

    After HDBSCAN, positional outliers are trimmed from each cluster
    using the IQR method. Breakpoint positions are refined using CIGAR
    soft-clip consensus when available, falling back to median.
    """
    if reads_df.empty or "pos" not in reads_df.columns:
        return [], pd.DataFrame()

    positions = np.array(reads_df["pos"], dtype=float).reshape(-1, 1)

    if len(positions) < cluster_params.min_cluster_size:
        logger.info(
            "Too few reads (%d) for clustering (min_cluster_size=%d)",
            len(positions),
            cluster_params.min_cluster_size,
        )
        return [], pd.DataFrame()

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=cluster_params.min_cluster_size,
        min_samples=cluster_params.min_samples,
        cluster_selection_epsilon=cluster_params.cluster_selection_epsilon,
        allow_single_cluster=cluster_params.allow_single_cluster,
        cluster_selection_method="leaf",
    )
    labels = clusterer.fit_predict(positions)
    probabilities = clusterer.probabilities_

    reads_df = reads_df.copy()
    reads_df["cluster"] = labels
    reads_df["probability"] = probabilities

    # Trim positional outliers from each cluster (IQR method)
    before_trim = int((labels != -1).sum())
    reads_df = _trim_cluster_outliers(reads_df)
    clustered = reads_df[reads_df["cluster"] != -1]
    n_trimmed = before_trim - len(clustered)
    if n_trimmed > 0:
        logger.info("Trimmed %d outlier reads from clusters", n_trimmed)

    # Refine wide clusters: extract tight sub-clusters via gap-based splitting
    before_filter = clustered["cluster"].nunique()
    reads_df = _refine_wide_clusters(
        reads_df,
        cluster_params.max_cluster_span,
        min_subcluster=cluster_params.min_cluster_size,
    )
    clustered = reads_df[reads_df["cluster"] != -1]
    after_filter = clustered["cluster"].nunique() if not clustered.empty else 0
    n_wide = before_filter - after_filter  # net clusters lost (may be negative if sub-clusters were extracted)
    if n_wide > 0:
        logger.info(
            "Refined %d wide cluster(s) exceeding %d bp span",
            n_wide,
            cluster_params.max_cluster_span,
        )
    elif after_filter > before_filter:
        logger.info(
            "Extracted %d sub-cluster(s) from wide clusters",
            after_filter - before_filter,
        )

    if clustered.empty:
        logger.info("No clusters found by HDBSCAN")
        return [], pd.DataFrame()

    n_clusters = clustered["cluster"].nunique()
    logger.info(
        "HDBSCAN found %d cluster(s) with %d reads",
        n_clusters,
        len(clustered),
    )

    breakpoints: list[Breakpoint] = []
    for _cluster_id, group in clustered.groupby("cluster"):
        chrom = group["chrom"].iloc[0]
        breakpoints.append(Breakpoint(
            chrom=chrom,
            position=_refine_breakpoint_position(group),
            pos_min=int(group["pos"].min()),
            pos_max=int(group["pos"].max()),
            supporting_reads=len(group),
            confidence=float(group["probability"].mean()),
            acro_specificity=None,
            read_positions=sorted(group["pos"].astype(int).tolist()),
        ))

    return breakpoints, reads_df
