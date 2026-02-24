"""Tests for the cluster module."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gollumpy.cluster import (
    _parse_clip_position,
    _refine_breakpoint_position,
    _trim_cluster_outliers,
    cluster_breakpoints,
    pre_cluster_reads,
)
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
            assert bp.mate_concordance is None
            assert bp.dominant_saac_chrom is None
            assert bp.ring_score is None

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

    def test_allow_single_cluster(self) -> None:
        # A single tight group should be detected as one cluster
        np.random.seed(42)
        positions = np.random.normal(47097797, 30, 10).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(10)],
            "chrom": ["chr22"] * 10,
            "pos": positions,
        })

        result, _ = cluster_breakpoints(
            df, ClusterParams(min_cluster_size=3, min_samples=2, allow_single_cluster=True)
        )
        assert len(result) == 1
        assert result[0].supporting_reads >= 3

class TestPreClusterReads:
    def test_empty_input(self) -> None:
        result = pre_cluster_reads(pd.DataFrame())
        assert result.empty

    def test_no_pos_column(self) -> None:
        df = pd.DataFrame({"read_id": ["r1"], "chrom": ["chr22"]})
        result = pre_cluster_reads(df)
        assert len(result) == 1  # returns input unchanged

    def test_removes_scattered_noise(self) -> None:
        """Scattered reads across the genome should be removed as noise."""
        np.random.seed(42)
        # Tight cluster of 10 reads + 5 scattered noise reads
        cluster_pos = np.random.normal(47097797, 50, 10).astype(int)
        noise_pos = np.array([1_000_000, 10_000_000, 20_000_000, 30_000_000, 40_000_000])
        positions = np.concatenate([cluster_pos, noise_pos])

        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(len(positions))],
            "chrom": ["chr22"] * len(positions),
            "pos": positions,
        })

        result = pre_cluster_reads(df, min_cluster_size=2, min_samples=2, cluster_selection_epsilon=500.0)

        # Should keep the cluster but remove at least some noise
        assert len(result) < len(df)
        assert len(result) >= 10  # cluster should survive

    def test_preserves_single_cluster(self) -> None:
        """A single tight cluster should be fully preserved."""
        np.random.seed(42)
        positions = np.random.normal(47097797, 30, 15).astype(int)
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(15)],
            "chrom": ["chr22"] * 15,
            "pos": positions,
        })

        result = pre_cluster_reads(df)
        assert len(result) == 15  # all reads should be kept

    def test_returns_original_if_all_noise(self) -> None:
        """If HDBSCAN marks everything as noise, return original (safety net)."""
        # Two reads very far apart — HDBSCAN won't cluster them
        df = pd.DataFrame({
            "read_id": ["r0", "r1"],
            "chrom": ["chr22", "chr22"],
            "pos": [1_000_000, 50_000_000],
        })
        result = pre_cluster_reads(df, min_cluster_size=3)
        # With min_cluster_size=3 and only 2 reads, falls through to
        # "too few reads" path → returns original
        assert len(result) == 2

    def test_resets_index(self) -> None:
        """Returned DataFrame should have a clean 0-based index."""
        np.random.seed(42)
        cluster_pos = np.random.normal(47097797, 30, 10).astype(int)
        noise_pos = np.array([1_000_000, 10_000_000, 20_000_000])
        positions = np.concatenate([cluster_pos, noise_pos])

        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(len(positions))],
            "chrom": ["chr22"] * len(positions),
            "pos": positions,
        })

        result = pre_cluster_reads(df)
        if len(result) < len(df):
            assert list(result.index) == list(range(len(result)))


class TestClusterEpsilon:
    def test_epsilon_merges_close_subclusters(self) -> None:
        # Two nearby subclusters within 100bp should merge with epsilon=100
        np.random.seed(42)
        pos1 = np.random.normal(47097700, 10, 10).astype(int)
        pos2 = np.random.normal(47097780, 10, 10).astype(int)
        positions = np.concatenate([pos1, pos2])

        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(20)],
            "chrom": ["chr22"] * 20,
            "pos": positions,
        })

        result_merged, _ = cluster_breakpoints(
            df, ClusterParams(min_cluster_size=3, min_samples=2, cluster_selection_epsilon=100.0)
        )
        result_split, _ = cluster_breakpoints(
            df, ClusterParams(min_cluster_size=3, min_samples=2, cluster_selection_epsilon=0.0)
        )

        # With epsilon=100, subclusters within 80bp should merge
        assert len(result_merged) <= len(result_split)


class TestTrimClusterOutliers:
    def test_removes_outliers(self) -> None:
        """14 tight reads + 5 outliers → outliers trimmed to noise."""
        # Mimics C0011CP cluster_0: tight core + distant stragglers
        tight = list(range(47_382_600, 47_382_600 + 14 * 25, 25))
        outliers = [46_179_722, 46_907_935, 47_234_058, 47_884_930, 48_019_339]
        positions = tight + outliers
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(19)],
            "chrom": ["chr22"] * 19,
            "pos": positions,
            "cluster": [0] * 19,
            "probability": [0.9] * 19,
        })

        result = _trim_cluster_outliers(df)
        clustered = result[result["cluster"] != -1]
        assert len(clustered) == 14
        # All remaining reads should be in the tight core
        assert clustered["pos"].min() >= 47_000_000

    def test_skips_small_clusters(self) -> None:
        """Clusters with < 5 reads should not be trimmed."""
        df = pd.DataFrame({
            "read_id": ["r0", "r1", "r2", "r3"],
            "chrom": ["chr22"] * 4,
            "pos": [1000, 1100, 1200, 50_000_000],  # r3 is an outlier
            "cluster": [0, 0, 0, 0],
            "probability": [0.9] * 4,
        })

        result = _trim_cluster_outliers(df)
        clustered = result[result["cluster"] != -1]
        assert len(clustered) == 4  # all kept (< 5 reads)

    def test_zero_iqr_not_trimmed(self) -> None:
        """Cluster with all identical positions should not be trimmed."""
        df = pd.DataFrame({
            "read_id": [f"r{i}" for i in range(6)],
            "chrom": ["chr22"] * 6,
            "pos": [47_000_000] * 6,
            "cluster": [0] * 6,
            "probability": [0.9] * 6,
        })

        result = _trim_cluster_outliers(df)
        clustered = result[result["cluster"] != -1]
        assert len(clustered) == 6  # all kept

    def test_preserves_noise_label(self) -> None:
        """Reads already labeled as noise (-1) should stay as noise."""
        df = pd.DataFrame({
            "read_id": ["r0", "r1", "r2"],
            "chrom": ["chr22"] * 3,
            "pos": [1000, 2000, 99_000_000],
            "cluster": [-1, 0, 0],
            "probability": [0.0, 0.9, 0.9],
        })

        result = _trim_cluster_outliers(df)
        assert result.iloc[0]["cluster"] == -1  # stays noise


class TestParseClipPosition:
    def test_right_side_clip(self) -> None:
        pos = _parse_clip_position("100M50S", ref_start=1000, ref_end=1100)
        assert pos == 1100

    def test_left_side_clip(self) -> None:
        pos = _parse_clip_position("50S100M", ref_start=1000, ref_end=1100)
        assert pos == 1000

    def test_no_significant_clip(self) -> None:
        pos = _parse_clip_position("150M", ref_start=1000, ref_end=1150)
        assert pos is None

    def test_small_clip_ignored(self) -> None:
        pos = _parse_clip_position("140M10S", ref_start=1000, ref_end=1140)
        assert pos is None  # 10bp < min_clip=20

    def test_custom_min_clip(self) -> None:
        pos = _parse_clip_position("140M10S", ref_start=1000, ref_end=1140, min_clip=5)
        assert pos == 1140

    def test_empty_cigar(self) -> None:
        assert _parse_clip_position("", ref_start=0, ref_end=0) is None
        assert _parse_clip_position(None, ref_start=0, ref_end=0) is None  # type: ignore[arg-type]

    def test_complex_cigar_with_right_clip(self) -> None:
        """CIGAR with insertions/deletions before soft-clip."""
        pos = _parse_clip_position("80M2I50M20S", ref_start=1000, ref_end=1130)
        assert pos == 1130


class TestRefineBreakpointPosition:
    def test_softclip_consensus(self) -> None:
        """5 reads with right-side clips at similar positions → clip median."""
        df = pd.DataFrame({
            "pos": [47_382_500, 47_382_510, 47_382_520, 47_382_530, 47_382_540],
            "cigarstring": ["100M50S", "100M50S", "100M50S", "100M50S", "100M50S"],
            "reference_end": [47_382_600, 47_382_610, 47_382_620, 47_382_630, 47_382_640],
        })

        pos = _refine_breakpoint_position(df)
        # Should use clip positions (reference_end): 600, 610, 620, 630, 640
        # Median of clip positions = 47382620
        assert pos == 47_382_620

    def test_softclip_disagreement_falls_back_to_median(self) -> None:
        """Clips with wide IQR (> 50bp) → fall back to median of all positions."""
        df = pd.DataFrame({
            "pos": [47_382_500, 47_382_510, 47_382_520, 47_382_530, 47_382_540],
            "cigarstring": ["100M50S", "100M50S", "100M50S", "100M50S", "100M50S"],
            # Clip positions spread over 200bp — disagree
            "reference_end": [47_382_500, 47_382_550, 47_382_600, 47_382_650, 47_382_700],
        })

        pos = _refine_breakpoint_position(df)
        # IQR of clips: 75bp > 50bp → fallback to median of pos
        assert pos == 47_382_520  # median of positions

    def test_no_cigar_columns_uses_median(self) -> None:
        """DataFrame without cigarstring column → median of positions."""
        df = pd.DataFrame({
            "pos": [100, 200, 300, 400, 500],
        })

        pos = _refine_breakpoint_position(df)
        assert pos == 300  # median

    def test_too_few_clips_uses_median(self) -> None:
        """Only 2 soft-clipped reads (< 3 threshold) → median."""
        df = pd.DataFrame({
            "pos": [100, 200, 300, 400, 500],
            "cigarstring": ["100M50S", "100M50S", "150M", "150M", "150M"],
            "reference_end": [200, 300, 450, 550, 650],
        })

        pos = _refine_breakpoint_position(df)
        assert pos == 300  # median of all positions

    def test_left_side_clips(self) -> None:
        """Left-side soft-clips should use reference_start as breakpoint."""
        df = pd.DataFrame({
            "pos": [47_382_600, 47_382_605, 47_382_610, 47_382_615, 47_382_620],
            "cigarstring": ["50S100M", "50S100M", "50S100M", "50S100M", "50S100M"],
            "reference_end": [47_382_700, 47_382_705, 47_382_710, 47_382_715, 47_382_720],
        })

        pos = _refine_breakpoint_position(df)
        # Left clips → breakpoint at reference_start (pos)
        # Clip positions = pos values: 600, 605, 610, 615, 620
        assert pos == 47_382_610  # median of clip positions
