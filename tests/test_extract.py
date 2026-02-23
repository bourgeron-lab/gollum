"""Tests for the extract module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from gollumpy.config import GollumConfig
from gollumpy.extract import extract_discordant_reads_t2t, extract_mate_sequences


def _make_config(tmp_path: Path, **overrides: object) -> GollumConfig:
    """Create a GollumConfig with temp files."""
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


def _mock_read(
    *,
    query_name: str = "read1",
    reference_name: str = "chr22",
    reference_start: int = 47000000,
    mapping_quality: int = 60,
    is_proper_pair: bool = False,
    is_unmapped: bool = False,
    mate_is_unmapped: bool = False,
    next_reference_name: str | None = "chr22",
    next_reference_start: int = 5000000,
    query_sequence: str = "ACGTACGTACGTACGT",
) -> MagicMock:
    """Create a mock pysam read."""
    read = MagicMock()
    read.query_name = query_name
    read.reference_name = reference_name
    read.reference_start = reference_start
    read.mapping_quality = mapping_quality
    read.is_proper_pair = is_proper_pair
    read.is_unmapped = is_unmapped
    read.mate_is_unmapped = mate_is_unmapped
    read.next_reference_name = next_reference_name
    read.next_reference_start = next_reference_start
    read.query_sequence = query_sequence
    return read


class TestExtractDiscordantReadsT2T:
    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_discordant_reads(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Mock a discordant read with mate in SAAC
        good_read = _mock_read(
            query_name="good_read",
            reference_start=47097797,
            next_reference_name="chr22",
            next_reference_start=5000000,  # within SAAC (< 14200000)
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [good_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)

        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "good_read"
        assert df.iloc[0]["pos"] == 47097797

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_proper_pairs(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        proper_read = _mock_read(is_proper_pair=True)

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [proper_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_low_mapq(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        low_mq_read = _mock_read(mapping_quality=30)

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [low_mq_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_mate_outside_saac(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Mate position beyond centromere boundary (14200000 for chr22)
        non_saac_read = _mock_read(
            next_reference_name="chr22",
            next_reference_start=20000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [non_saac_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_non_acrocentric_mate(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        # Mate on non-acrocentric chromosome
        non_acro_read = _mock_read(
            next_reference_name="chr1",
            next_reference_start=5000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [non_acro_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_blacklisted_positions(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        # Create blacklist file
        bed = tmp_path / "blacklist.bed"
        bed.write_text("chr22\t47097000\t47098000\n")

        config = _make_config(tmp_path, blacklist_bed=bed)

        # Read position falls within blacklisted region
        blacklisted_read = _mock_read(reference_start=47097500)

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [blacklisted_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_returns_empty_df_with_correct_columns(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = []
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert df.empty
        assert list(df.columns) == ["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"]


class TestExtractMateSequences:
    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_sequences(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        config = _make_config(tmp_path)

        reads_df = pd.DataFrame([{
            "read_id": "read1",
            "chrom": "chr22",
            "pos": 47097797,
            "mapq": 60,
            "mate_chrom": "chr22",
            "mate_pos": 5000000,
        }])

        mate_read = _mock_read(query_name="read1", query_sequence="AAACCCGGGTTT")

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [mate_read]
        mock_bam_class.return_value = mock_bam

        result = extract_mate_sequences(config, reads_df)
        assert len(result) == 1
        assert result.iloc[0]["mate_sequence"] == "AAACCCGGGTTT"

    def test_empty_input(self, tmp_path: Path) -> None:
        config = _make_config(tmp_path)
        empty_df = pd.DataFrame(columns=["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"])
        result = extract_mate_sequences(config, empty_df)
        assert result.empty
        assert "mate_sequence" in result.columns
