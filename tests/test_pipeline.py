"""Tests for the pipeline and report modules."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from gollumpy.config import GollumConfig
from gollumpy.models import Breakpoint
from gollumpy.pipeline import compute_ring_score, run_pipeline
from gollumpy.report import generate_report


def _make_config(tmp_path: Path, **overrides: object) -> GollumConfig:
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    cram = tmp_path / "test.cram"
    cram.touch()
    defaults = {
        "target_chrom": "chr22",
        "mode": "t2t",
        "fasta_t2t": fasta,
        "input_file": cram,
        "output_dir": tmp_path / "out",
        "sample_name": "test_sample",
    }
    defaults.update(overrides)
    return GollumConfig(**defaults)  # type: ignore[arg-type]


class TestGenerateReport:
    def test_with_breakpoints(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        breakpoints = [
            Breakpoint(
                chrom="chr22",
                position=47097797,
                pos_min=47097637,
                pos_max=47097906,
                supporting_reads=11,
                confidence=0.95,
                acro_specificity=9.7,
            ),
        ]
        reads_df = pd.DataFrame({
            "read_id": ["r1"],
            "chrom": ["chr22"],
            "pos": [47097797],
        })

        generate_report(breakpoints, reads_df, config)

        summary_path = config.output_dir / "test_sample.summary.tsv"
        assert summary_path.exists()
        content = summary_path.read_text()
        assert "test_sample" in content
        assert "47097797" in content

        reads_path = config.output_dir / "test_sample.supporting_reads.tsv"
        assert reads_path.exists()

    def test_no_breakpoints(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        generate_report([], pd.DataFrame(), config)

        summary_path = config.output_dir / "test_sample.summary.tsv"
        assert summary_path.exists()
        assert "no ring detected" in summary_path.read_text()

    def test_none_specificity(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        breakpoints = [
            Breakpoint(
                chrom="chr22",
                position=47097797,
                pos_min=47097637,
                pos_max=47097906,
                supporting_reads=5,
                confidence=0.8,
                acro_specificity=None,
            ),
        ]

        generate_report(breakpoints, pd.DataFrame(), config)

        summary_path = config.output_dir / "test_sample.summary.tsv"
        content = summary_path.read_text()
        assert "NA" in content


class TestRunPipeline:
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_no_discordant_reads(self, mock_extract: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        mock_extract.return_value = pd.DataFrame(
            columns=["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"]
        )

        result = run_pipeline(config)
        assert result == []

    @patch("gollumpy.pipeline.pre_cluster_reads")
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_pre_cluster_removes_all_reads(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        tmp_path: Path,
    ) -> None:
        """If pre-clustering removes all reads, pipeline returns empty."""
        config = _make_config(tmp_path)
        reads_df = pd.DataFrame([{
            "read_id": "r1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
        }])
        mock_extract.return_value = reads_df
        mock_precluster.return_value = pd.DataFrame(
            columns=["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"]
        )

        result = run_pipeline(config)
        assert result == []
        mock_precluster.assert_called_once()

    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.pre_cluster_reads", side_effect=lambda df: df)
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_no_alignments(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        tmp_path: Path,
    ) -> None:
        config = _make_config(tmp_path)

        reads_df = pd.DataFrame([{
            "read_id": "r1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
        }])
        mock_extract.return_value = reads_df
        mock_mate_seq.return_value = reads_df.assign(mate_sequence="ACGT")
        mock_align.return_value = pd.DataFrame()

        result = run_pipeline(config)
        assert result == []

    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.pre_cluster_reads", side_effect=lambda df: df)
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_phr_filter_removes_reads_without_target_saac_hit(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Reads with no SAAC hit on the target chrom are filtered out."""
        config = _make_config(tmp_path)  # target_chrom = chr22

        reads_df = pd.DataFrame([
            {"read_id": f"r{i}", "chrom": "chr22", "pos": 47097797 + i, "mapq": 60,
             "mate_chrom": "chr22", "mate_pos": 5000000}
            for i in range(5)
        ])
        mock_extract.return_value = reads_df
        mock_mate_seq.return_value = reads_df.assign(mate_sequence="ACGT" * 30)

        # All reads have SAAC hits on chr15 only, not the target (chr22)
        aligned_df = reads_df.assign(
            mate_sequence="ACGT" * 30,
            best_mlen=148,
            best_divergence=0.01,
            n_saac_hits=1,
            best_align_chrom="chr15",
            has_target_saac_hit=False,  # no hit on target chrom
        )
        mock_align.return_value = aligned_df

        result = run_pipeline(config)
        # All reads filtered by PHR filter → no ring detected
        assert result == []
        assert (config.output_dir / "test_sample.summary.tsv").exists()

    @patch("gollumpy.pipeline.build_aligner")
    @patch("gollumpy.pipeline.compute_mate_concordance")
    @patch("gollumpy.pipeline.cluster_breakpoints")
    @patch("gollumpy.pipeline.compute_acro_specificity")
    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.pre_cluster_reads", side_effect=lambda df: df)
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_full_pipeline_with_detection(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        mock_specificity: MagicMock,
        mock_cluster: MagicMock,
        mock_concordance: MagicMock,
        mock_build_aligner: MagicMock,
        tmp_path: Path,
    ) -> None:
        config = _make_config(tmp_path)

        reads_data = [
            {
                "read_id": f"r{i}",
                "chrom": "chr22",
                "pos": 47097797 + i * 10,
                "mapq": 60,
                "mate_chrom": "chr22",
                "mate_pos": 5000000 + i,
            }
            for i in range(15)
        ]
        reads_df = pd.DataFrame(reads_data)
        mock_extract.return_value = reads_df
        mock_mate_seq.return_value = reads_df.assign(mate_sequence="ACGT" * 30)

        aligned_df = reads_df.assign(
            mate_sequence="ACGT" * 30,
            best_mlen=148,
            best_divergence=0.01,
            n_saac_hits=1,
            best_align_chrom="chr22",
            has_target_saac_hit=True,
        )
        mock_align.return_value = aligned_df

        # cluster_breakpoints returns (breakpoints, labeled_df)
        bp = Breakpoint(
            chrom="chr22",
            position=47097867,
            pos_min=47097797,
            pos_max=47097937,
            supporting_reads=15,
            confidence=0.95,
            acro_specificity=None,
        )
        labeled_df = aligned_df.copy()
        labeled_df["cluster"] = 0
        labeled_df["probability"] = 0.95
        mock_cluster.return_value = ([bp], labeled_df)

        # Per-cluster specificity and concordance
        mock_specificity.return_value = 9.7
        mock_concordance.return_value = (0.95, "chr22")
        mock_build_aligner.return_value = MagicMock()

        result = run_pipeline(config)
        assert len(result) == 1
        assert result[0].chrom == "chr22"
        assert result[0].supporting_reads == 15
        assert result[0].acro_specificity == 9.7
        assert result[0].mate_concordance == 0.95
        assert result[0].dominant_saac_chrom == "chr22"
        assert result[0].ring_score is not None
        assert result[0].ring_score > 0

        # Verify compute_acro_specificity was called once (one cluster)
        mock_specificity.assert_called_once()
        mock_concordance.assert_called_once()

        # Check output files were created
        assert (config.output_dir / "test_sample.summary.tsv").exists()
        assert (config.output_dir / "test_sample.log").exists()

    @patch("gollumpy.pipeline.cluster_breakpoints")
    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.pre_cluster_reads", side_effect=lambda df: df)
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_no_clusters(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        mock_cluster: MagicMock,
        tmp_path: Path,
    ) -> None:
        config = _make_config(tmp_path)

        reads_df = pd.DataFrame([{
            "read_id": "r1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
        }])
        mock_extract.return_value = reads_df
        mock_mate_seq.return_value = reads_df.assign(mate_sequence="ACGT" * 30)
        mock_align.return_value = reads_df.assign(
            mate_sequence="ACGT" * 30,
            best_mlen=148,
            best_divergence=0.01,
            n_saac_hits=1,
            best_align_chrom="chr22",
            has_target_saac_hit=True,
        )
        mock_cluster.return_value = ([], pd.DataFrame())

        result = run_pipeline(config)
        assert result == []
        assert (config.output_dir / "test_sample.summary.tsv").exists()

    @patch("gollumpy.pipeline.build_aligner")
    @patch("gollumpy.pipeline.compute_mate_concordance")
    @patch("gollumpy.pipeline.cluster_breakpoints")
    @patch("gollumpy.pipeline.compute_acro_specificity")
    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.pre_cluster_reads", side_effect=lambda df: df)
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_min_supporting_reads_filter(
        self,
        mock_extract: MagicMock,
        mock_precluster: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        mock_specificity: MagicMock,
        mock_cluster: MagicMock,
        mock_concordance: MagicMock,
        mock_build_aligner: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Clusters with fewer reads than min_supporting_reads are filtered out."""
        from gollumpy.config import ClusterParams

        config = _make_config(tmp_path)
        config.cluster_params = ClusterParams(min_supporting_reads=10)

        reads_data = [
            {
                "read_id": f"r{i}",
                "chrom": "chr22",
                "pos": 47097797 + i * 10,
                "mapq": 60,
                "mate_chrom": "chr22",
                "mate_pos": 5000000 + i,
            }
            for i in range(20)
        ]
        reads_df = pd.DataFrame(reads_data)
        mock_extract.return_value = reads_df
        mock_mate_seq.return_value = reads_df.assign(mate_sequence="ACGT" * 30)

        aligned_df = reads_df.assign(
            mate_sequence="ACGT" * 30,
            best_mlen=148,
            best_divergence=0.01,
            n_saac_hits=1,
            best_align_chrom="chr22",
            has_target_saac_hit=True,
        )
        mock_align.return_value = aligned_df

        # Two clusters: one with 15 reads (passes), one with 5 reads (filtered)
        bp_big = Breakpoint(
            chrom="chr22", position=47097867, pos_min=47097797, pos_max=47097937,
            supporting_reads=15, confidence=0.95, acro_specificity=None,
        )
        bp_small = Breakpoint(
            chrom="chr22", position=10000000, pos_min=9999900, pos_max=10000100,
            supporting_reads=5, confidence=0.80, acro_specificity=None,
        )
        labeled_df = aligned_df.copy()
        labeled_df["cluster"] = [0] * 15 + [1] * 5
        labeled_df["probability"] = 0.9
        mock_cluster.return_value = ([bp_big, bp_small], labeled_df)

        mock_specificity.return_value = 9.7
        mock_concordance.return_value = (0.95, "chr22")
        mock_build_aligner.return_value = MagicMock()

        result = run_pipeline(config)

        # Only the big cluster should survive the filter
        assert len(result) == 1
        assert result[0].supporting_reads == 15

        # compute_acro_specificity should be called only once (for the surviving cluster)
        mock_specificity.assert_called_once()


class TestComputeRingScore:
    def test_perfect_signal(self) -> None:
        """High reads, confidence, specificity, concordance, tight span → high score."""
        score = compute_ring_score(
            supporting_reads=100,
            confidence=1.0,
            acro_specificity=float("inf"),
            mate_concordance=1.0,
            span=200,  # tight cluster
        )
        assert score > 8.0
        assert score <= 10.0

    def test_noise_signal(self) -> None:
        """Low reads, confidence, specificity, concordance, wide span → low score."""
        score = compute_ring_score(
            supporting_reads=3,
            confidence=0.3,
            acro_specificity=1.5,
            mate_concordance=0.3,
            span=5_000_000,  # wide cluster
        )
        assert score < 3.0

    def test_none_specificity(self) -> None:
        """None specificity should contribute 0 to the score."""
        score = compute_ring_score(
            supporting_reads=10,
            confidence=0.8,
            acro_specificity=None,
            mate_concordance=0.9,
            span=500,
        )
        assert score > 0
        # Compare with non-None specificity (should be lower)
        score_with_spec = compute_ring_score(
            supporting_reads=10,
            confidence=0.8,
            acro_specificity=10.0,
            mate_concordance=0.9,
            span=500,
        )
        assert score < score_with_spec

    def test_inf_specificity(self) -> None:
        """Infinite specificity should contribute max (1.0) to the spec component."""
        score_inf = compute_ring_score(
            supporting_reads=10,
            confidence=0.8,
            acro_specificity=float("inf"),
            mate_concordance=0.9,
            span=500,
        )
        score_high = compute_ring_score(
            supporting_reads=10,
            confidence=0.8,
            acro_specificity=20.0,
            mate_concordance=0.9,
            span=500,
        )
        # inf and 20.0 should both cap at 1.0
        assert score_inf == pytest.approx(score_high)

    def test_score_range(self) -> None:
        """Ring score should always be between 0 and 10."""
        score = compute_ring_score(
            supporting_reads=0,
            confidence=0.0,
            acro_specificity=None,
            mate_concordance=0.0,
            span=0,
        )
        assert score >= 0.0
        assert score <= 10.0

    def test_span_penalty(self) -> None:
        """Wide span should produce a lower score than tight span."""
        tight = compute_ring_score(
            supporting_reads=10, confidence=0.8, acro_specificity=10.0,
            mate_concordance=0.9, span=200,
        )
        wide = compute_ring_score(
            supporting_reads=10, confidence=0.8, acro_specificity=10.0,
            mate_concordance=0.9, span=5_000_000,
        )
        assert tight > wide

    def test_concordance_and_span_tied_highest_weight(self) -> None:
        """Mate concordance (0.25) and span_tightness (0.25) are top-weighted components."""
        base = compute_ring_score(
            supporting_reads=10, confidence=0.5, acro_specificity=5.0,
            mate_concordance=0.0, span=500,
        )
        with_concordance = compute_ring_score(
            supporting_reads=10, confidence=0.5, acro_specificity=5.0,
            mate_concordance=1.0, span=500,
        )
        # Concordance adds 0.25 * 1.0 * 10 = 2.5 points
        assert with_concordance - base == pytest.approx(2.5)

    def test_linear_read_signal(self) -> None:
        """read_signal uses linear scaling: 5 reads → 0.33, 10 → 0.67, 15+ → 1.0."""
        score_5 = compute_ring_score(
            supporting_reads=5, confidence=0.5, acro_specificity=10.0,
            mate_concordance=0.5, span=500,
        )
        score_10 = compute_ring_score(
            supporting_reads=10, confidence=0.5, acro_specificity=10.0,
            mate_concordance=0.5, span=500,
        )
        score_15 = compute_ring_score(
            supporting_reads=15, confidence=0.5, acro_specificity=10.0,
            mate_concordance=0.5, span=500,
        )
        score_30 = compute_ring_score(
            supporting_reads=30, confidence=0.5, acro_specificity=10.0,
            mate_concordance=0.5, span=500,
        )
        # Linear scaling: 10 reads should be substantially higher than 5
        assert score_10 - score_5 > 0.5
        # 15 and 30 reads should be capped at the same value
        assert score_15 == pytest.approx(score_30)

    def test_confidence_low_weight(self) -> None:
        """Confidence weight is 0.10 — small clusters with high confidence shouldn't dominate."""
        # Simulate: 5-read noise cluster with confidence=0.98
        noise = compute_ring_score(
            supporting_reads=5, confidence=0.98, acro_specificity=15.0,
            mate_concordance=0.4, span=200,
        )
        # Simulate: 15-read real cluster with confidence=0.55
        real = compute_ring_score(
            supporting_reads=15, confidence=0.55, acro_specificity=15.0,
            mate_concordance=0.7, span=500_000,
        )
        # Real cluster should still outscore noise despite lower confidence
        assert real > noise
