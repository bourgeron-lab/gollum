"""Tests for config module."""

from pathlib import Path

import pytest

from gollumpy.config import (
    ACROCENTRIC_CHROMS,
    ChromRegion,
    ClusterParams,
    FilterParams,
    GollumConfig,
    is_blacklisted,
    load_acro_regions,
    load_blacklist,
)
from gollumpy.models import Breakpoint, DiscordantRead


class TestChromRegion:
    def test_creation(self) -> None:
        region = ChromRegion(chrom="chr22", saac_start=0, saac_end=14200000)
        assert region.chrom == "chr22"
        assert region.saac_start == 0
        assert region.saac_end == 14200000

    def test_frozen(self) -> None:
        region = ChromRegion(chrom="chr22", saac_start=0, saac_end=14200000)
        with pytest.raises(AttributeError):
            region.saac_end = 999  # type: ignore[misc]


class TestFilterParams:
    def test_defaults(self) -> None:
        params = FilterParams()
        assert params.min_mapq_r1 == 40
        assert params.samtools_mapq == 40
        assert params.min_alignment_score == 80
        assert params.max_divergence == 0.06
        assert params.centromere_buffer == 5_000_000
        assert params.require_target_saac is False

    def test_custom(self) -> None:
        params = FilterParams(min_mapq_r1=30, max_divergence=0.1)
        assert params.min_mapq_r1 == 30
        assert params.max_divergence == 0.1


class TestClusterParams:
    def test_defaults(self) -> None:
        params = ClusterParams()
        assert params.min_cluster_size == 3
        assert params.min_samples == 5
        assert params.cluster_selection_epsilon == 100.0
        assert params.allow_single_cluster is True
        assert params.min_supporting_reads == 3
        assert params.max_cluster_span == 10_000


class TestGollumConfig:
    def test_invalid_chromosome(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        with pytest.raises(ValueError, match="Invalid chromosome"):
            GollumConfig(
                target_chrom="chr1",
                mode="t2t",
                fasta_t2t=fasta,
                input_file=cram,
                output_dir=tmp_path / "out",
                sample_name="test",
            )

    def test_missing_input(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        with pytest.raises(FileNotFoundError, match="Input file not found"):
            GollumConfig(
                target_chrom="chr22",
                mode="t2t",
                fasta_t2t=fasta,
                input_file=tmp_path / "missing.cram",
                output_dir=tmp_path / "out",
                sample_name="test",
            )

    def test_missing_fasta(self, tmp_path: Path) -> None:
        cram = tmp_path / "test.cram"
        cram.touch()
        with pytest.raises(FileNotFoundError, match="T2T reference not found"):
            GollumConfig(
                target_chrom="chr22",
                mode="t2t",
                fasta_t2t=tmp_path / "missing.fa",
                input_file=cram,
                output_dir=tmp_path / "out",
                sample_name="test",
            )

    def test_missing_blacklist(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        with pytest.raises(FileNotFoundError, match="Blacklist BED not found"):
            GollumConfig(
                target_chrom="chr22",
                mode="t2t",
                fasta_t2t=fasta,
                input_file=cram,
                output_dir=tmp_path / "out",
                sample_name="test",
                blacklist_bed=tmp_path / "missing.bed",
            )

    def test_creates_output_dir(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        out = tmp_path / "deep" / "nested" / "output"
        config = GollumConfig(
            target_chrom="chr22",
            mode="t2t",
            fasta_t2t=fasta,
            input_file=cram,
            output_dir=out,
            sample_name="test",
        )
        assert config.output_dir.exists()

    def test_valid_config(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        config = GollumConfig(
            target_chrom="chr22",
            mode="t2t",
            fasta_t2t=fasta,
            input_file=cram,
            output_dir=tmp_path / "out",
            sample_name="test_sample",
        )
        assert config.target_chrom == "chr22"
        assert config.mode == "t2t"
        assert config.sample_name == "test_sample"

    def test_reference_fasta_t2t_mode(self, tmp_path: Path) -> None:
        fasta = tmp_path / "t2t.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        config = GollumConfig(
            target_chrom="chr22",
            mode="t2t",
            fasta_t2t=fasta,
            input_file=cram,
            output_dir=tmp_path / "out",
            sample_name="test",
        )
        assert config.reference_fasta == fasta

    def test_reference_fasta_grch38_mode(self, tmp_path: Path) -> None:
        t2t_fasta = tmp_path / "t2t.fa"
        t2t_fasta.touch()
        grch38_fasta = tmp_path / "grch38.fa"
        grch38_fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        config = GollumConfig(
            target_chrom="chr22",
            mode="grch38",
            fasta_t2t=t2t_fasta,
            input_file=cram,
            output_dir=tmp_path / "out",
            sample_name="test",
            fasta_grch38=grch38_fasta,
        )
        assert config.reference_fasta == grch38_fasta

    def test_grch38_mode_missing_grch38_fasta(self, tmp_path: Path) -> None:
        fasta = tmp_path / "t2t.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        with pytest.raises(ValueError, match="GRCh38 reference"):
            GollumConfig(
                target_chrom="chr22",
                mode="grch38",
                fasta_t2t=fasta,
                input_file=cram,
                output_dir=tmp_path / "out",
                sample_name="test",
            )

    def test_grch38_mode_missing_grch38_fasta_file(self, tmp_path: Path) -> None:
        fasta = tmp_path / "t2t.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        with pytest.raises(FileNotFoundError, match="GRCh38 reference not found"):
            GollumConfig(
                target_chrom="chr22",
                mode="grch38",
                fasta_t2t=fasta,
                input_file=cram,
                output_dir=tmp_path / "out",
                sample_name="test",
                fasta_grch38=tmp_path / "missing.fa",
            )


class TestLoadAcroRegions:
    def test_loads_all_chromosomes(self) -> None:
        regions = load_acro_regions()
        assert set(regions.keys()) == set(ACROCENTRIC_CHROMS)

    def test_chr22_values(self) -> None:
        regions = load_acro_regions()
        assert regions["chr22"].saac_start == 0
        assert regions["chr22"].saac_end == 14200000

    def test_all_regions_have_valid_boundaries(self) -> None:
        regions = load_acro_regions()
        for chrom, region in regions.items():
            assert region.chrom == chrom
            assert region.saac_start >= 0
            assert region.saac_end > region.saac_start


class TestBlacklist:
    def test_load_blacklist(self, tmp_path: Path) -> None:
        bed = tmp_path / "blacklist.bed"
        bed.write_text("chr22\t1000\t2000\nchr22\t5000\t6000\n")
        regions = load_blacklist(bed)
        assert len(regions) == 2
        assert regions[0] == ("chr22", 1000, 2000)
        assert regions[1] == ("chr22", 5000, 6000)

    def test_load_blacklist_with_comments(self, tmp_path: Path) -> None:
        bed = tmp_path / "blacklist.bed"
        bed.write_text("# comment\nchr22\t1000\t2000\n\nchr22\t5000\t6000\n")
        regions = load_blacklist(bed)
        assert len(regions) == 2

    def test_is_blacklisted_true(self) -> None:
        blacklist = [("chr22", 1000, 2000), ("chr22", 5000, 6000)]
        assert is_blacklisted("chr22", 1500, blacklist)
        assert is_blacklisted("chr22", 1000, blacklist)  # inclusive start
        assert is_blacklisted("chr22", 5999, blacklist)

    def test_is_blacklisted_false(self) -> None:
        blacklist = [("chr22", 1000, 2000), ("chr22", 5000, 6000)]
        assert not is_blacklisted("chr22", 999, blacklist)
        assert not is_blacklisted("chr22", 2000, blacklist)  # exclusive end
        assert not is_blacklisted("chr22", 3000, blacklist)
        assert not is_blacklisted("chr13", 1500, blacklist)  # different chrom

    def test_is_blacklisted_empty(self) -> None:
        assert not is_blacklisted("chr22", 1500, [])


class TestModels:
    def test_discordant_read(self) -> None:
        read = DiscordantRead(
            read_id="read1",
            chrom="chr22",
            pos=47097797,
            mapq=60,
            mate_chrom="chr22",
            mate_pos=5000000,
            mate_sequence="ACGTACGT",
        )
        assert read.read_id == "read1"
        assert read.pos == 47097797

    def test_breakpoint(self) -> None:
        bp = Breakpoint(
            chrom="chr22",
            position=47097797,
            pos_min=47097637,
            pos_max=47097906,
            supporting_reads=11,
            confidence=0.95,
            acro_specificity=9.7,
        )
        assert bp.supporting_reads == 11
        assert bp.confidence == 0.95

    def test_breakpoint_new_fields_default_none(self) -> None:
        bp = Breakpoint(
            chrom="chr22",
            position=47097797,
            pos_min=47097637,
            pos_max=47097906,
            supporting_reads=11,
            confidence=0.95,
            acro_specificity=9.7,
        )
        assert bp.mate_concordance is None
        assert bp.dominant_saac_chrom is None
        assert bp.ring_score is None

    def test_breakpoint_with_all_fields(self) -> None:
        bp = Breakpoint(
            chrom="chr22",
            position=47097797,
            pos_min=47097637,
            pos_max=47097906,
            supporting_reads=11,
            confidence=0.95,
            acro_specificity=9.7,
            mate_concordance=0.91,
            dominant_saac_chrom="chr22",
            ring_score=7.5,
        )
        assert bp.mate_concordance == 0.91
        assert bp.dominant_saac_chrom == "chr22"
        assert bp.ring_score == 7.5
