"""Tests for the extract module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from gollumpy.config import FilterParams, GollumConfig
from gollumpy.extract import (
    extract_discordant_reads_grch38,
    extract_discordant_reads_t2t,
    extract_mate_sequences,
)


def _make_config(tmp_path: Path, **overrides: object) -> GollumConfig:
    """Create a GollumConfig with temp files."""
    fasta = tmp_path / "ref.fa"
    fasta.touch()
    cram = tmp_path / "test.cram"
    cram.touch()
    defaults: dict[str, object] = {
        "target_chrom": "chr22",
        "mode": "t2t",
        "fasta_t2t": fasta,
        "input_file": cram,
        "output_dir": tmp_path / "out",
        "sample_name": "test",
    }
    defaults.update(overrides)
    # GRCh38 mode requires fasta_grch38
    if defaults["mode"] == "grch38" and "fasta_grch38" not in defaults:
        grch38_fasta = tmp_path / "grch38.fa"
        grch38_fasta.touch()
        defaults["fasta_grch38"] = grch38_fasta
    return GollumConfig(**defaults)  # type: ignore[arg-type]


def _mock_read(
    *,
    query_name: str = "read1",
    reference_name: str = "chr22",
    reference_start: int = 47000000,
    reference_end: int | None = None,
    cigarstring: str = "150M",
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
    read.reference_end = reference_end if reference_end is not None else reference_start + 150
    read.cigarstring = cigarstring
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
        assert df.iloc[0]["cigarstring"] == "150M"
        assert df.iloc[0]["reference_end"] == 47097797 + 150

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
    def test_centromere_buffer_shifts_fetch_start(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Centromere buffer pushes extraction start past pericentromeric zone."""
        config = _make_config(tmp_path, filter_params=FilterParams(centromere_buffer=5_000_000))

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = []
        mock_bam_class.return_value = mock_bam

        extract_discordant_reads_t2t(config)

        # chr22 SAAC ends at 14200000; with 5Mb buffer, fetch should start at 19200000
        mock_bam.fetch.assert_called_once_with("chr22", 19_200_000)

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_zero_buffer_extracts_from_saac_boundary(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """With centromere_buffer=0, extraction starts at SAAC boundary."""
        config = _make_config(tmp_path, filter_params=FilterParams(centromere_buffer=0))

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = []
        mock_bam_class.return_value = mock_bam

        extract_discordant_reads_t2t(config)

        # chr22 SAAC ends at 14200000; with 0 buffer, fetch starts at 14200000
        mock_bam.fetch.assert_called_once_with("chr22", 14_200_000)

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
        expected_cols = [
            "read_id", "chrom", "pos", "cigarstring", "reference_end",
            "mapq", "mate_chrom", "mate_pos",
        ]
        assert list(df.columns) == expected_cols


    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_default_blacklist_filters_known_noise_region(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Default blacklist auto-loads and filters chr22:38861950-38862750."""
        config = _make_config(tmp_path)  # no blacklist_bed → uses default
        assert config.use_default_blacklist is True

        # Read at a position inside the default blacklist (38861950-38862750)
        blacklisted_read = _mock_read(
            query_name="noise_read",
            reference_start=38862100,
            next_reference_name="chr13",
            next_reference_start=5000000,
        )
        # Read at a position outside the blacklist
        good_read = _mock_read(
            query_name="good_read",
            reference_start=47097797,
            next_reference_name="chr22",
            next_reference_start=5000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [blacklisted_read, good_read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "good_read"

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_no_blacklist_skips_filtering(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """With use_default_blacklist=False, no blacklist filtering occurs."""
        config = _make_config(tmp_path, use_default_blacklist=False)

        # Read at a position inside the default blacklist region
        read_in_blacklist = _mock_read(
            query_name="should_pass",
            reference_start=38862100,
            next_reference_name="chr13",
            next_reference_start=5000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read_in_blacklist]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_t2t(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "should_pass"


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


class TestExtractDiscordantReadsGRCh38:
    """Tests for GRCh38 discordant read extraction."""

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_on_acrocentric_short_arm(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on chr21 short arm (pos < 12M) should be kept."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="short_arm_read",
            reference_start=47000000,
            next_reference_name="chr21",
            next_reference_start=8000000,  # chr21 short arm (< 12M)
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "short_arm_read"
        assert df.iloc[0]["mate_chrom"] == "chr21"
        assert df.iloc[0]["mate_pos"] == 8000000

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_on_chr22_short_arm(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on chr22 short arm (pos < 15M) should be kept."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="chr22_short_arm",
            reference_start=47000000,
            next_reference_name="chr22",
            next_reference_start=11834000,  # chr22 short arm (< 15M)
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "chr22_short_arm"

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_on_non_canonical_contig(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on named decoy contig (KMT2C_chr21_*) should be kept."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="decoy_read",
            reference_start=47000000,
            next_reference_name="KMT2C_chr21_7687010_7731520",
            next_reference_start=21335,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "decoy_read"
        assert df.iloc[0]["mate_chrom"] == "KMT2C_chr21_7687010_7731520"

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_unmapped(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Unmapped mate should be kept."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="unmapped_mate",
            reference_start=47000000,
            mate_is_unmapped=True,
            next_reference_name=None,
            next_reference_start=0,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["read_id"] == "unmapped_mate"
        assert df.iloc[0]["mate_chrom"] == "*"
        assert df.iloc[0]["mate_pos"] == 0

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_on_chrUn(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on chrUn_* contig should be kept (regression)."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="chrUn_read",
            reference_start=47000000,
            next_reference_name="chrUn_GL000220v1",
            next_reference_start=133412,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["mate_chrom"] == "chrUn_GL000220v1"

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_extracts_mate_on_random_contig(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on *_random contig should be kept (regression)."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            query_name="random_read",
            reference_start=47000000,
            next_reference_name="chr22_KI270733v1_random",
            next_reference_start=151211,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert len(df) == 1
        assert df.iloc[0]["mate_chrom"] == "chr22_KI270733v1_random"

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_proper_pairs(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Proper pairs should be filtered."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            is_proper_pair=True,
            next_reference_name="chr21",
            next_reference_start=8000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_low_mapq(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Low mapQ R1 reads should be filtered."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            mapping_quality=30,
            next_reference_name="chr21",
            next_reference_start=8000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_mate_on_qarm(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on q-arm of acrocentric chromosome (beyond short arm) should be filtered."""
        config = _make_config(tmp_path, mode="grch38")

        # Mate at chr22:40000000 — well beyond q-arm start (15M), not short arm
        read = _mock_read(
            next_reference_name="chr22",
            next_reference_start=40000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_filters_mate_on_non_acrocentric_canonical(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Mate on canonical non-acrocentric chromosome (e.g. chr1) should be filtered."""
        config = _make_config(tmp_path, mode="grch38")

        read = _mock_read(
            next_reference_name="chr1",
            next_reference_start=50000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_returns_empty_df_with_correct_columns(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Empty result should have correct column schema."""
        config = _make_config(tmp_path, mode="grch38")

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = []
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty
        expected_cols = [
            "read_id", "chrom", "pos", "cigarstring", "reference_end",
            "mapq", "mate_chrom", "mate_pos",
        ]
        assert list(df.columns) == expected_cols

    @patch("gollumpy.extract.pysam.AlignmentFile")
    def test_blacklist_filters_grch38(self, mock_bam_class: MagicMock, tmp_path: Path) -> None:
        """Blacklisted R1 positions should be filtered in GRCh38 mode."""
        bed = tmp_path / "blacklist.bed"
        bed.write_text("chr22\t47000000\t47001000\n")

        config = _make_config(tmp_path, mode="grch38", blacklist_bed=bed)

        read = _mock_read(
            reference_start=47000500,  # inside blacklist
            next_reference_name="chr21",
            next_reference_start=8000000,
        )

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.fetch.return_value = [read]
        mock_bam_class.return_value = mock_bam

        df = extract_discordant_reads_grch38(config)
        assert df.empty
