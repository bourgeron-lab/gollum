"""Tests for the pipeline and report modules."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from gollumpy.config import GollumConfig
from gollumpy.models import Breakpoint
from gollumpy.pipeline import run_pipeline
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

    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_no_alignments(
        self,
        mock_extract: MagicMock,
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

    @patch("gollumpy.pipeline.build_aligner")
    @patch("gollumpy.pipeline.cluster_breakpoints")
    @patch("gollumpy.pipeline.compute_acro_specificity")
    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_full_pipeline_with_detection(
        self,
        mock_extract: MagicMock,
        mock_mate_seq: MagicMock,
        mock_align: MagicMock,
        mock_specificity: MagicMock,
        mock_cluster: MagicMock,
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

        # Per-cluster specificity
        mock_specificity.return_value = 9.7
        mock_build_aligner.return_value = MagicMock()

        result = run_pipeline(config)
        assert len(result) == 1
        assert result[0].chrom == "chr22"
        assert result[0].supporting_reads == 15
        assert result[0].acro_specificity == 9.7

        # Verify compute_acro_specificity was called once (one cluster)
        mock_specificity.assert_called_once()

        # Check output files were created
        assert (config.output_dir / "test_sample.summary.tsv").exists()
        assert (config.output_dir / "test_sample.log").exists()

    @patch("gollumpy.pipeline.cluster_breakpoints")
    @patch("gollumpy.pipeline.align_mates")
    @patch("gollumpy.pipeline.extract_mate_sequences")
    @patch("gollumpy.pipeline.extract_discordant_reads_t2t")
    def test_no_clusters(
        self,
        mock_extract: MagicMock,
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
        )
        mock_cluster.return_value = ([], pd.DataFrame())

        result = run_pipeline(config)
        assert result == []
        assert (config.output_dir / "test_sample.summary.tsv").exists()
