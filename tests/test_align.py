"""Tests for the align module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from gollumpy.align import align_mates, build_aligner, compute_acro_specificity, compute_mate_concordance
from gollumpy.config import GollumConfig


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
        "sample_name": "test",
    }
    defaults.update(overrides)
    return GollumConfig(**defaults)  # type: ignore[arg-type]


def _mock_hit(
    *,
    ctg: str = "chr22",
    r_st: int = 5000000,
    r_en: int = 5000150,
    mapq: int = 0,
    NM: int = 2,
    blen: int = 150,
    mlen: int = 148,
    is_primary: bool = True,
) -> MagicMock:
    hit = MagicMock()
    hit.ctg = ctg
    hit.r_st = r_st
    hit.r_en = r_en
    hit.mapq = mapq
    hit.NM = NM
    hit.blen = blen
    hit.mlen = mlen
    hit.is_primary = is_primary
    return hit


class TestAlignMates:
    def test_empty_input(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        result = align_mates(pd.DataFrame(), config)
        assert result.empty

    def test_no_mate_sequence_column(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        df = pd.DataFrame({"read_id": ["r1"], "pos": [100]})
        result = align_mates(df, config)
        assert result.empty

    @patch("gollumpy.align.mappy.Aligner")
    def test_saac_hit_passes(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Hit in chr22 SAAC region (< 14200000)
        saac_hit = _mock_hit(ctg="chr22", r_st=5000000, r_en=5000150, mlen=148, NM=2, blen=150)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [saac_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert len(result) == 1
        assert result.iloc[0]["read_id"] == "read1"
        assert "best_mlen" in result.columns
        assert "best_align_chrom" in result.columns
        assert "has_target_saac_hit" in result.columns
        assert result.iloc[0]["best_align_chrom"] == "chr22"
        assert bool(result.iloc[0]["has_target_saac_hit"]) is True

    @patch("gollumpy.align.mappy.Aligner")
    def test_non_target_saac_hit_still_in_output(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        """A read with SAAC hit on non-target chrom is kept but has_target_saac_hit=False."""
        config = _make_config(tmp_path)  # target_chrom = chr22

        # Hit in chr14 SAAC region (not the target chr22)
        chr14_hit = _mock_hit(ctg="chr14", r_st=5000000, r_en=5000150, mlen=148, NM=2, blen=150)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [chr14_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr14",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert len(result) == 1
        assert result.iloc[0]["best_align_chrom"] == "chr14"
        assert bool(result.iloc[0]["has_target_saac_hit"]) is False

    @patch("gollumpy.align.mappy.Aligner")
    def test_multi_hit_with_any_target_saac(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        """A read with best hit on chr14 but also a hit on chr22 SAAC gets has_target_saac_hit=True."""
        config = _make_config(tmp_path)  # target_chrom = chr22

        # Best hit on chr14 (higher mlen), secondary hit on chr22
        chr14_hit = _mock_hit(ctg="chr14", r_st=5000000, r_en=5000150, mlen=148, NM=2, blen=150)
        chr22_hit = _mock_hit(ctg="chr22", r_st=5000000, r_en=5000140, mlen=138, NM=3, blen=140)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [chr14_hit, chr22_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr14",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert len(result) == 1
        assert result.iloc[0]["best_align_chrom"] == "chr14"  # best hit is on chr14
        assert bool(result.iloc[0]["has_target_saac_hit"]) is True  # but has ANY hit on chr22

    @patch("gollumpy.align.mappy.Aligner")
    def test_non_saac_hit_filtered(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Hit outside SAAC (position > 14200000 on chr22)
        non_saac_hit = _mock_hit(ctg="chr22", r_st=20000000, r_en=20000150)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [non_saac_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert result.empty

    @patch("gollumpy.align.mappy.Aligner")
    def test_non_acrocentric_hit_filtered(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Hit on non-acrocentric chromosome
        chr1_hit = _mock_hit(ctg="chr1", r_st=5000000, r_en=5000150)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [chr1_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert result.empty

    @patch("gollumpy.align.mappy.Aligner")
    def test_low_quality_hit_filtered(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Hit with low mlen (< 100 default threshold)
        low_qual_hit = _mock_hit(ctg="chr22", r_st=5000000, r_en=5000050, mlen=50, NM=1, blen=50)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [low_qual_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert result.empty

    @patch("gollumpy.align.mappy.Aligner")
    def test_high_divergence_filtered(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Hit with high divergence (NM/blen = 20/150 > 0.05)
        divergent_hit = _mock_hit(ctg="chr22", r_st=5000000, r_en=5000150, mlen=130, NM=20, blen=150)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [divergent_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
            "mate_sequence": "ACGTACGTACGT",
        }])

        result = align_mates(reads_df, config)
        assert result.empty


class TestComputeAcroSpecificity:
    def test_empty_input(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        assert compute_acro_specificity(pd.DataFrame(), config) is None

    @patch("gollumpy.align.mappy.Aligner")
    def test_all_acrocentric(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        acro_hit = _mock_hit(ctg="chr22", mlen=150)
        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [acro_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "mate_sequence": "ACGTACGT",
        }])

        result = compute_acro_specificity(reads_df, config)
        assert result == float("inf")

    @patch("gollumpy.align.mappy.Aligner")
    def test_mixed_hits(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        acro_hit = _mock_hit(ctg="chr22", mlen=150)
        non_acro_hit = _mock_hit(ctg="chr1", mlen=50)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [acro_hit, non_acro_hit]
        mock_aligner_class.return_value = mock_aligner

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "mate_sequence": "ACGTACGT",
        }])

        result = compute_acro_specificity(reads_df, config)
        assert result == pytest.approx(3.0)  # 150 / 50

    @patch("gollumpy.align.mappy.Aligner")
    def test_with_prebuilt_aligner(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        acro_hit = _mock_hit(ctg="chr22", mlen=150)
        non_acro_hit = _mock_hit(ctg="chr1", mlen=50)

        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner.map.return_value = [acro_hit, non_acro_hit]

        reads_df = pd.DataFrame([{"read_id": "read1", "mate_sequence": "ACGTACGT"}])

        result = compute_acro_specificity(reads_df, config, aligner=mock_aligner)
        assert result == pytest.approx(3.0)
        # Aligner constructor should NOT have been called (we passed our own)
        mock_aligner_class.assert_not_called()


class TestBuildAligner:
    @patch("gollumpy.align.mappy.Aligner")
    def test_success(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=True)
        mock_aligner_class.return_value = mock_aligner

        result = build_aligner(config)
        assert result is mock_aligner
        mock_aligner_class.assert_called_once_with(str(config.fasta_t2t), preset="sr", best_n=5)

    @patch("gollumpy.align.mappy.Aligner")
    def test_failure(self, mock_aligner_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        mock_aligner = MagicMock()
        mock_aligner.__bool__ = MagicMock(return_value=False)
        mock_aligner_class.return_value = mock_aligner

        with pytest.raises(RuntimeError, match="Failed to build minimap2 index"):
            build_aligner(config)


class TestComputeMateConcordance:
    def test_empty_input(self) -> None:
        concordance, dominant = compute_mate_concordance(pd.DataFrame())
        assert concordance == 0.0
        assert dominant == ""

    def test_missing_column(self) -> None:
        df = pd.DataFrame({"read_id": ["r1"]})
        concordance, dominant = compute_mate_concordance(df)
        assert concordance == 0.0
        assert dominant == ""

    def test_perfect_concordance(self) -> None:
        df = pd.DataFrame({
            "read_id": ["r1", "r2", "r3", "r4", "r5"],
            "best_align_chrom": ["chr22", "chr22", "chr22", "chr22", "chr22"],
        })
        concordance, dominant = compute_mate_concordance(df)
        assert concordance == 1.0
        assert dominant == "chr22"

    def test_mixed_concordance(self) -> None:
        df = pd.DataFrame({
            "read_id": ["r1", "r2", "r3", "r4", "r5"],
            "best_align_chrom": ["chr22", "chr22", "chr22", "chr14", "chr15"],
        })
        concordance, dominant = compute_mate_concordance(df)
        assert concordance == pytest.approx(0.6)
        assert dominant == "chr22"

    def test_scattered_noise(self) -> None:
        df = pd.DataFrame({
            "read_id": ["r1", "r2", "r3", "r4"],
            "best_align_chrom": ["chr22", "chr14", "chr15", "chr13"],
        })
        concordance, dominant = compute_mate_concordance(df)
        assert concordance == pytest.approx(0.25)

    def test_single_read(self) -> None:
        df = pd.DataFrame({
            "read_id": ["r1"],
            "best_align_chrom": ["chr22"],
        })
        concordance, dominant = compute_mate_concordance(df)
        assert concordance == 1.0
        assert dominant == "chr22"
