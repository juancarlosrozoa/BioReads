"""
Platform detector.

Infers the sequencing platform and experiment type from FASTQ read headers
or from user-provided hints.

Illumina headers:  @A00123:123:HHMMTDRXX:1:1101:1234:1234 1:N:0:ATCG
Nanopore headers:  @<uuid> runid=... read=... ch=... start_time=...
PacBio headers:    @m64012_190920_173625/1/ccs
Ion Torrent:       @<flow_cell>:<x>:<y>

Usage:
    from bioreads.core.detector import detect_platform

    platform, confidence = detect_platform("reads.fastq")
    # → ("illumina", "high")
"""

from __future__ import annotations

import re
from pathlib import Path
from Bio import SeqIO


# Header patterns per platform
_PATTERNS = {
    "illumina":   re.compile(r"^@[A-Z0-9]+:\d+:[A-Z0-9]+:\d+:\d+:\d+:\d+"),
    "nanopore":   re.compile(r"runid=|read=\d+.*ch=\d+|@[0-9a-f-]{36}"),
    "pacbio":     re.compile(r"^@m\d+_\d+_\d+/\d+/(ccs|subreads?)"),
    "iontorrent": re.compile(r"^@[A-Z0-9]+:\d+:\d+$"),
    "bgi":        re.compile(r"^@[A-Z0-9]+L\d+C\d+R\d+/"),
}

# Read length thresholds
_SHORT_READ_MAX = 600   # bp — Illumina / Ion Torrent
_LONG_READ_MIN  = 1000  # bp — Nanopore / PacBio


def detect_platform(
    fastq_file: str | Path,
    n_reads: int = 500,
) -> tuple[str, str]:
    """Infer sequencing platform from a FASTQ file.

    Reads up to `n_reads` records and scores header patterns + read lengths.

    Args:
        fastq_file: Path to a FASTQ file (.fastq / .fq / .gz).
        n_reads:    Number of reads to sample (default: 500).

    Returns:
        Tuple of (platform, confidence) where platform is one of:
        "illumina", "nanopore", "pacbio", "iontorrent", "bgi", "unknown"
        and confidence is "high", "medium", or "low".
    """
    path = Path(fastq_file)
    fmt  = "fastq"

    scores: dict[str, int] = {k: 0 for k in _PATTERNS}
    lengths: list[int] = []
    total = 0

    try:
        for record in SeqIO.parse(path, fmt):
            header = record.description
            lengths.append(len(record.seq))
            for platform, pattern in _PATTERNS.items():
                if pattern.search(header):
                    scores[platform] += 1
            total += 1
            if total >= n_reads:
                break
    except Exception:
        return "unknown", "low"

    if total == 0:
        return "unknown", "low"

    # Pick platform with highest score
    best = max(scores, key=scores.__getitem__)
    best_score = scores[best]

    # Fallback: use read length if no header matched
    if best_score == 0 and lengths:
        mean_len = sum(lengths) / len(lengths)
        if mean_len >= _LONG_READ_MIN:
            return "nanopore", "low"
        return "illumina", "low"

    confidence = "high" if best_score / total > 0.8 else "medium" if best_score / total > 0.4 else "low"
    return best, confidence


def suggest_aligner(platform: str, experiment: str = "dna") -> dict:
    """Suggest the best aligner for a platform + experiment combination.

    Args:
        platform:   Platform string from detect_platform().
        experiment: "dna", "rna", or "amplicon".

    Returns:
        Dict with keys: aligner, preset, level, fallback.
        - aligner:  recommended external binary
        - preset:   aligner-specific preset flag
        - level:    1 (pip-only) or 2 (external binary required)
        - fallback: pip-installable fallback (mappy)
    """
    matrix = {
        ("illumina",   "dna"):      {"aligner": "bwa-mem2",  "preset": None,          "level": 2, "fallback": None},
        ("illumina",   "rna"):      {"aligner": "STAR",      "preset": None,          "level": 2, "fallback": None},
        ("illumina",   "amplicon"): {"aligner": "bwa-mem2",  "preset": None,          "level": 2, "fallback": None},
        ("nanopore",   "dna"):      {"aligner": "minimap2",  "preset": "map-ont",     "level": 2, "fallback": "mappy"},
        ("nanopore",   "rna"):      {"aligner": "minimap2",  "preset": "splice",      "level": 2, "fallback": "mappy"},
        ("pacbio",     "dna"):      {"aligner": "minimap2",  "preset": "map-hifi",    "level": 2, "fallback": "mappy"},
        ("pacbio",     "rna"):      {"aligner": "minimap2",  "preset": "splice:hq",   "level": 2, "fallback": "mappy"},
        ("iontorrent", "dna"):      {"aligner": "bwa-mem2",  "preset": None,          "level": 2, "fallback": None},
        ("bgi",        "dna"):      {"aligner": "bwa-mem2",  "preset": None,          "level": 2, "fallback": None},
        ("bgi",        "rna"):      {"aligner": "STAR",      "preset": None,          "level": 2, "fallback": None},
    }
    key = (platform, experiment)
    return matrix.get(key, {"aligner": "minimap2", "preset": "map-ont", "level": 2, "fallback": "mappy"})
