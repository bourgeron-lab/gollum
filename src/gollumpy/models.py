"""Data models for gollum pipeline."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class DiscordantRead:
    """A discordant read pair where R1 anchors on q-arm and R2 maps to SAAC."""

    read_id: str
    chrom: str
    pos: int
    mapq: int
    mate_chrom: str
    mate_pos: int
    mate_sequence: str


@dataclass
class Breakpoint:
    """A detected ring chromosome breakpoint supported by clustered discordant reads."""

    chrom: str
    position: int
    pos_min: int
    pos_max: int
    supporting_reads: int
    confidence: float
    acro_specificity: float | None
    mate_concordance: float | None = None
    dominant_saac_chrom: str | None = None
    ring_score: float | None = None
