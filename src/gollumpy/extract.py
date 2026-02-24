"""Discordant read extraction from BAM/CRAM files using pysam."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd
import pysam

from gollumpy.config import ACROCENTRIC_CHROMS, is_blacklisted, load_acro_regions, load_blacklist

if TYPE_CHECKING:
    from gollumpy.config import GollumConfig

logger = logging.getLogger(__name__)


def extract_discordant_reads_t2t(config: GollumConfig) -> pd.DataFrame:
    """Extract discordant reads from a T2T-aligned BAM/CRAM.

    Reads are extracted from the target chromosome q-arm region.
    Keeps reads where R1 anchors on q-arm (mapQ >= threshold) and
    R2 mate maps to any acrocentric SAAC region.
    """
    acro_regions = load_acro_regions()
    target_region = acro_regions[config.target_chrom]

    blacklist: list[tuple[str, int, int]] = []
    if config.blacklist_bed is not None:
        blacklist = load_blacklist(config.blacklist_bed)

    reads: list[dict[str, str | int]] = []

    with pysam.AlignmentFile(str(config.input_file), reference_filename=str(config.reference_fasta)) as bam:
        # Fetch from q-arm: skip pericentromeric zone (buffer) to avoid artifacts
        fetch_start = target_region.saac_end + config.filter_params.centromere_buffer
        buf = config.filter_params.centromere_buffer
        logger.info("Extracting from %s:%d (buffer=%d)", config.target_chrom, fetch_start, buf)
        for read in bam.fetch(config.target_chrom, fetch_start):
            # Skip proper pairs — we want discordant reads
            if read.is_proper_pair:
                continue

            # Skip unmapped reads or reads without mate info
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            # Mapping quality filter (samtools -q 40 equivalent)
            if read.mapping_quality < config.filter_params.samtools_mapq:
                continue

            # Stricter R1 mapQ filter
            if read.mapping_quality < config.filter_params.min_mapq_r1:
                continue

            # Get mate chromosome
            mate_chrom = read.next_reference_name
            if mate_chrom is None:
                continue

            mate_pos = read.next_reference_start

            # Check if mate maps to any acrocentric SAAC region
            if mate_chrom not in ACROCENTRIC_CHROMS:
                continue

            mate_region = acro_regions.get(mate_chrom)
            if mate_region is None or mate_pos > mate_region.saac_end:
                continue

            # Blacklist filter on R1 position
            if blacklist and is_blacklisted(config.target_chrom, read.reference_start, blacklist):
                continue

            # Get mate sequence from the read's query sequence
            # Note: in paired-end data, pysam gives us the read's own sequence,
            # not the mate's. We store the read's sequence for now and will
            # extract mate sequences separately.
            reads.append({
                "read_id": read.query_name,
                "chrom": read.reference_name,
                "pos": read.reference_start,
                "mapq": read.mapping_quality,
                "mate_chrom": mate_chrom,
                "mate_pos": mate_pos,
            })

    logger.info("Discordant reads extracted: %d", len(reads))

    if not reads:
        return pd.DataFrame(columns=["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"])

    return pd.DataFrame(reads)


def extract_mate_sequences(
    config: GollumConfig,
    reads_df: pd.DataFrame,
) -> pd.DataFrame:
    """Extract mate sequences for discordant reads from BAM/CRAM.

    For each discordant read, fetches the mate's sequence using pysam.
    Returns the input DataFrame with an added 'mate_sequence' column.
    """
    if reads_df.empty:
        reads_df["mate_sequence"] = pd.Series(dtype=str)
        return reads_df

    mate_sequences: dict[str, str] = {}

    with pysam.AlignmentFile(str(config.input_file), reference_filename=str(config.reference_fasta)) as bam:
        # Fetch mates from their mapped positions
        for _, row in reads_df.iterrows():
            mate_chrom = row["mate_chrom"]
            mate_pos = int(row["mate_pos"])
            read_id = row["read_id"]

            if read_id in mate_sequences:
                continue

            try:
                for read in bam.fetch(mate_chrom, mate_pos, mate_pos + 1):
                    if read.query_name == read_id and read.query_sequence is not None:
                        mate_sequences[read_id] = read.query_sequence
                        break
            except ValueError:
                logger.warning("Could not fetch mate for %s at %s:%d", read_id, mate_chrom, mate_pos)
                continue

    reads_df = reads_df.copy()
    reads_df["mate_sequence"] = reads_df["read_id"].map(mate_sequences)

    # Drop reads where we couldn't get mate sequence
    before = len(reads_df)
    reads_df = reads_df.dropna(subset=["mate_sequence"])
    dropped = before - len(reads_df)
    if dropped > 0:
        logger.warning("Dropped %d reads without mate sequence", dropped)

    logger.info("Reads with mate sequences: %d", len(reads_df))
    return reads_df


def extract_discordant_reads_grch38(config: GollumConfig) -> pd.DataFrame:
    """Extract discordant reads from a GRCh38-aligned BAM/CRAM.

    In GRCh38, SAAC sequences are absent, so R2 reads from ring breakpoints will be:
    - Unmapped
    - Mapped to chrUn_* or alt contigs with low quality
    - Mapped elsewhere with poor alignment (mapQ < 10)
    """
    blacklist: list[tuple[str, int, int]] = []
    if config.blacklist_bed is not None:
        blacklist = load_blacklist(config.blacklist_bed)

    # GRCh38 approximate centromere boundaries for q-arm extraction
    # These are rough boundaries — reads will be re-aligned to T2T anyway
    grch38_qarm_starts: dict[str, int] = {
        "chr13": 17700000,
        "chr14": 17200000,
        "chr15": 19000000,
        "chr21": 12000000,
        "chr22": 15000000,
    }

    qarm_start = grch38_qarm_starts.get(config.target_chrom)
    if qarm_start is None:
        msg = f"No GRCh38 q-arm start defined for {config.target_chrom}"
        raise ValueError(msg)

    reads: list[dict[str, str | int]] = []

    with pysam.AlignmentFile(str(config.input_file), reference_filename=str(config.reference_fasta)) as bam:
        fetch_start = qarm_start + config.filter_params.centromere_buffer
        buf = config.filter_params.centromere_buffer
        logger.info("Extracting from %s:%d (buffer=%d)", config.target_chrom, fetch_start, buf)
        for read in bam.fetch(config.target_chrom, fetch_start):
            if read.is_proper_pair:
                continue
            if read.is_unmapped:
                continue
            if read.mapping_quality < config.filter_params.min_mapq_r1:
                continue

            # In GRCh38 mode, keep reads where mate is:
            # 1. Unmapped
            # 2. Low mapQ (< 10)
            # 3. On chrUn_*, *_random, or alt contigs
            mate_chrom = read.next_reference_name
            mate_is_candidate = read.mate_is_unmapped or (
                mate_chrom is not None
                and (mate_chrom.startswith("chrUn_") or "_random" in mate_chrom or "_alt" in mate_chrom)
            )

            if not mate_is_candidate:
                continue

            if blacklist and is_blacklisted(config.target_chrom, read.reference_start, blacklist):
                continue

            reads.append({
                "read_id": read.query_name,
                "chrom": read.reference_name,
                "pos": read.reference_start,
                "mapq": read.mapping_quality,
                "mate_chrom": mate_chrom if mate_chrom is not None else "*",
                "mate_pos": read.next_reference_start if not read.mate_is_unmapped else 0,
            })

    logger.info("GRCh38 discordant reads extracted: %d", len(reads))

    if not reads:
        return pd.DataFrame(columns=["read_id", "chrom", "pos", "mapq", "mate_chrom", "mate_pos"])

    return pd.DataFrame(reads)
