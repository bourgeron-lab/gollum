"""Tests for the CLI module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from gollumpy.cli import main


@patch("gollumpy.cli.run_pipeline")
class TestCLI:
    def test_help(self, mock_pipeline) -> None:  # noqa: ANN001
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "Gollum" in result.output

    def test_t2t_help(self, mock_pipeline) -> None:  # noqa: ANN001
        runner = CliRunner()
        result = runner.invoke(main, ["t2t", "--help"])
        assert result.exit_code == 0
        assert "--fasta" in result.output
        assert "--chr" in result.output
        assert "--cluster-epsilon" in result.output
        assert "--allow-single" in result.output
        assert "--min-supporting-reads" in result.output
        assert "--max-cluster-span" in result.output
        assert "--min-cluster-span" in result.output
        assert "--require-target-saac" in result.output

    def test_grch38_help(self, mock_pipeline) -> None:  # noqa: ANN001
        runner = CliRunner()
        result = runner.invoke(main, ["grch38", "--help"])
        assert result.exit_code == 0
        assert "GRCh38" in result.output
        assert "--ref" in result.output

    def test_t2t_basic(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
        ])

        assert result.exit_code == 0
        mock_pipeline.assert_called_once()
        config = mock_pipeline.call_args[0][0]
        assert config.mode == "t2t"
        assert config.target_chrom == "chr22"
        assert config.sample_name == "sample"
        assert config.filter_params.min_alignment_score == 80
        assert config.filter_params.max_divergence == 0.06
        assert config.filter_params.centromere_buffer == 5_000_000
        assert config.cluster_params.min_cluster_size == 5
        assert config.cluster_params.min_samples == 3
        assert config.cluster_params.cluster_selection_epsilon == 100.0
        assert config.cluster_params.allow_single_cluster is True
        assert config.cluster_params.min_supporting_reads == 5
        assert config.cluster_params.max_cluster_span == 10_000
        assert config.cluster_params.min_cluster_span == 50
        assert config.filter_params.require_target_saac is False

    def test_t2t_custom_options(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--chr", "chr13",
            "-s", "my_sample",
            "--min-alignment-score", "150",
            "--max-divergence", "0.03",
            "--min-cluster-size", "5",
            "--min-samples", "3",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.target_chrom == "chr13"
        assert config.sample_name == "my_sample"
        assert config.filter_params.min_alignment_score == 150
        assert config.filter_params.max_divergence == 0.03
        assert config.cluster_params.min_cluster_size == 5
        assert config.cluster_params.min_samples == 3

    def test_grch38_basic(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        t2t_fasta = tmp_path / "t2t.fa"
        t2t_fasta.touch()
        grch38_fasta = tmp_path / "grch38.fa"
        grch38_fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "grch38", str(cram),
            "-f", str(t2t_fasta),
            "--ref", str(grch38_fasta),
            "-o", str(out),
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.mode == "grch38"
        assert config.fasta_t2t == t2t_fasta
        assert config.fasta_grch38 == grch38_fasta
        assert config.reference_fasta == grch38_fasta

    def test_grch38_missing_ref(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "t2t.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()

        runner = CliRunner()
        result = runner.invoke(main, [
            "grch38", str(cram),
            "-f", str(fasta),
            "-o", str(tmp_path / "out"),
            # missing --ref
        ])

        assert result.exit_code != 0

    def test_missing_fasta(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        cram = tmp_path / "test.cram"
        cram.touch()

        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(tmp_path / "missing.fa"),
            "-o", str(tmp_path / "out"),
        ])

        assert result.exit_code != 0

    def test_invalid_chromosome(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "test.cram"
        cram.touch()

        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(tmp_path / "out"),
            "--chr", "chr1",
        ])

        assert result.exit_code != 0

    def test_min_supporting_reads_option(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--min-supporting-reads", "10",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.cluster_params.min_supporting_reads == 10

    def test_centromere_buffer_option(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--centromere-buffer", "3000000",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.filter_params.centromere_buffer == 3_000_000

    def test_require_target_saac_flag(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--require-target-saac",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.filter_params.require_target_saac is True

    def test_max_cluster_span_option(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--max-cluster-span", "5000",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.cluster_params.max_cluster_span == 5000

    def test_min_cluster_span_option(self, mock_pipeline, tmp_path: Path) -> None:  # noqa: ANN001
        fasta = tmp_path / "ref.fa"
        fasta.touch()
        cram = tmp_path / "sample.cram"
        cram.touch()
        out = tmp_path / "out"

        mock_pipeline.return_value = []
        runner = CliRunner()
        result = runner.invoke(main, [
            "t2t", str(cram),
            "-f", str(fasta),
            "-o", str(out),
            "--min-cluster-span", "100",
        ])

        assert result.exit_code == 0
        config = mock_pipeline.call_args[0][0]
        assert config.cluster_params.min_cluster_span == 100

    def test_version(self, mock_pipeline) -> None:  # noqa: ANN001
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "2.0.0" in result.output
