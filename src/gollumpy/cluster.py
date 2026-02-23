"""Single-pass HDBSCAN clustering of breakpoint positions."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import hdbscan
import numpy as np
import pandas as pd

from gollumpy.models import Breakpoint

if TYPE_CHECKING:
    from gollumpy.config import ClusterParams

logger = logging.getLogger(__name__)


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
    """
    if reads_df.empty or "pos" not in reads_df.columns:
        return [], pd.DataFrame()

    positions = np.array(reads_df["pos"], dtype=float).reshape(-1, 1)

    if len(positions) < cluster_params.min_cluster_size:
        logger.info(
            "Too few reads (%d) for clustering (min_cluster_size=%d)", len(positions), cluster_params.min_cluster_size
        )
        return [], pd.DataFrame()

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=cluster_params.min_cluster_size,
        min_samples=cluster_params.min_samples,
    )
    labels = clusterer.fit_predict(positions)
    probabilities = clusterer.probabilities_

    reads_df = reads_df.copy()
    reads_df["cluster"] = labels
    reads_df["probability"] = probabilities

    # Filter out noise (label == -1)
    clustered = reads_df[reads_df["cluster"] != -1]

    if clustered.empty:
        logger.info("No clusters found by HDBSCAN")
        return [], pd.DataFrame()

    n_clusters = clustered["cluster"].nunique()
    logger.info("HDBSCAN found %d cluster(s) with %d reads", n_clusters, len(clustered))

    breakpoints: list[Breakpoint] = []
    for _cluster_id, group in clustered.groupby("cluster"):
        chrom = group["chrom"].iloc[0]
        breakpoints.append(Breakpoint(
            chrom=chrom,
            position=int(group["pos"].mean()),
            pos_min=int(group["pos"].min()),
            pos_max=int(group["pos"].max()),
            supporting_reads=len(group),
            confidence=float(group["probability"].mean()),
            acro_specificity=None,
        ))

    return breakpoints, reads_df
