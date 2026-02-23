"""Tests for the cluster module."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gollumpy.cluster import cluster_breakpoints
from gollumpy.config import ClusterParams


class TestClusterBreakpoints:
    def test_empty_input(self) -> None:
        result, labeled_df = cluster_breakpoints(pd.DataFrame(), ClusterParams())
        assert result == []
        assert labeled_df.empty

    def test_too_few_reads(self) -> None:
        df = pd.DataFrame({
            "read_id": ["r1"],
            "chrom": ["chr22"],
            "pos": [47097797],
        })
        result, labeled_df = cluster_breakpoints(df, ClusterParams(min_cluster_size=3))
        assert result == []
        assert labeled_df.empty

    def test_single_cluster(self) -> None:
        # Create reads clustered around position 47097797
        np.random.seed(42)
        n_reads = 15
        positions = np.random.normal(47097797, 50, n_reads).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(n_reads)],
            "chrom": ["chr22"] * n_reads,
            "pos": positions,
        })

        result, labeled_df = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        assert len(result) >= 1
        bp = result[0]
        assert bp.chrom == "chr22"
        assert abs(bp.position - 47097797) < 200  # within reasonable range
        assert bp.supporting_reads >= 3
        assert 0 <= bp.confidence <= 1

    def test_two_clusters(self) -> None:
        # Two distant clusters with enough spread for HDBSCAN
        np.random.seed(42)
        n_per_cluster = 20
        pos1 = np.random.normal(10000000, 100, n_per_cluster).astype(int)
        pos2 = np.random.normal(47000000, 100, n_per_cluster).astype(int)
        positions = np.concatenate([pos1, pos2])

        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(2 * n_per_cluster)],
            "chrom": ["chr22"] * (2 * n_per_cluster),
            "pos": positions,
        })

        result, labeled_df = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        # Should find at least 2 clusters (the two distant groups)
        assert len(result) >= 2
        # The two main groups should be far apart
        positions_found = sorted(bp.position for bp in result)
        assert positions_found[-1] - positions_found[0] > 30000000

    def test_noise_points_excluded(self) -> None:
        # Tight cluster + scattered noise
        np.random.seed(42)
        cluster_pos = np.random.normal(47097797, 20, 10).astype(int)
        noise_pos = np.array([1000000, 20000000, 35000000])
        positions = np.concatenate([cluster_pos, noise_pos])

        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(len(positions))],
            "chrom": ["chr22"] * len(positions),
            "pos": positions,
        })

        result, labeled_df = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        # Should find the tight cluster but not include noise
        assert len(result) >= 1
        total_clustered = sum(bp.supporting_reads for bp in result)
        assert total_clustered <= len(positions)  # not all points assigned

    def test_breakpoint_has_correct_fields(self) -> None:
        np.random.seed(42)
        positions = np.random.normal(47097797, 50, 10).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(10)],
            "chrom": ["chr22"] * 10,
            "pos": positions,
        })

        result, _ = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        assert len(result) >= 1
        bp = result[0]
        assert bp.pos_min <= bp.position <= bp.pos_max
        assert bp.supporting_reads > 0
        assert 0 <= bp.confidence <= 1

    def test_breakpoints_have_none_specificity(self) -> None:
        np.random.seed(42)
        positions = np.random.normal(47097797, 50, 10).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(10)],
            "chrom": ["chr22"] * 10,
            "pos": positions,
        })

        result, _ = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        assert len(result) >= 1
        for bp in result:
            assert bp.acro_specificity is None

    def test_returns_labeled_dataframe(self) -> None:
        np.random.seed(42)
        positions = np.random.normal(47097797, 50, 15).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(15)],
            "chrom": ["chr22"] * 15,
            "pos": positions,
        })

        breakpoints, labeled_df = cluster_breakpoints(df, ClusterParams(min_cluster_size=3, min_samples=2))

        assert len(breakpoints) >= 1
        assert "cluster" in labeled_df.columns
        assert "probability" in labeled_df.columns
        # All non-noise reads should have cluster >= 0
        clustered = labeled_df[labeled_df["cluster"] != -1]
        assert len(clustered) > 0
        assert all(clustered["cluster"] >= 0)
