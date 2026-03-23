"""
QC Engine.

Computes quality metrics for FASTQ files from any sequencing platform.
Works with pip-only dependencies (biopython).

Metrics computed:
  - Total reads and bases
  - Read length distribution (min, max, mean, N50)
  - Per-base quality scores (mean Phred)
  - % bases above Q20 and Q30
  - GC content
  - Adapter contamination hints (Illumina TruSeq / Nextera)
  - Duplicate rate estimate (first-read sampling)

Usage:
    from bioreads.core.qc import QCEngine

    engine = QCEngine()
    result = engine.run("reads.fastq", platform="illumina")
    print(result.summary())
    result.save_html("qc_report.html")
"""

from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path
from statistics import mean, median
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bioreads.core.detector import detect_platform


# Common Illumina adapter sequences (first 12 bp)
_ADAPTER_HINTS = {
    "TruSeq Universal": "AGATCGGAAGAGC",
    "Nextera":          "CTGTCTCTTATAC",
    "TruSeq Read2":     "AGATCGGAAGAGC",
}
_ADAPTER_CHECK_LEN = 13
_SAMPLE_READS      = 200_000   # reads to sample for QC


class QCResult:
    """Holds QC metrics for one FASTQ file."""

    def __init__(self, file: str, platform: str, platform_confidence: str):
        self.file                = file
        self.platform            = platform
        self.platform_confidence = platform_confidence

        self.total_reads:  int   = 0
        self.total_bases:  int   = 0
        self.min_length:   int   = 0
        self.max_length:   int   = 0
        self.mean_length:  float = 0.0
        self.median_length:float = 0.0
        self.n50:          int   = 0
        self.mean_quality: float = 0.0
        self.pct_q20:      float = 0.0
        self.pct_q30:      float = 0.0
        self.gc_content:   float = 0.0
        self.duplicate_rate: float = 0.0
        self.adapter_hits: dict[str, int] = {}
        self.length_dist:  dict[int, int] = {}   # binned

    def summary(self) -> str:
        lines = [
            f"File             : {self.file}",
            f"Platform         : {self.platform} (confidence: {self.platform_confidence})",
            f"",
            f"--- Read statistics ---",
            f"  Total reads    : {self.total_reads:,}",
            f"  Total bases    : {self.total_bases:,}",
            f"  Length min/max : {self.min_length} / {self.max_length} bp",
            f"  Length mean    : {self.mean_length:.1f} bp",
            f"  Length median  : {self.median_length:.1f} bp",
            f"  N50            : {self.n50:,} bp",
            f"",
            f"--- Quality ---",
            f"  Mean Phred Q   : {self.mean_quality:.1f}",
            f"  % bases >= Q20 : {self.pct_q20:.1f}%",
            f"  % bases >= Q30 : {self.pct_q30:.1f}%",
            f"  GC content     : {self.gc_content:.1f}%",
            f"  Duplicate rate : {self.duplicate_rate:.1f}%",
        ]
        if self.adapter_hits:
            lines += ["", "--- Adapter contamination ---"]
            for name, count in self.adapter_hits.items():
                pct = count / self.total_reads * 100 if self.total_reads else 0
                lines.append(f"  {name:<22}: {count:,} reads ({pct:.1f}%)")
        return "\n".join(lines)


class QCEngine:
    """Run quality control on a FASTQ file."""

    def run(
        self,
        fastq_file: str | Path,
        platform: str | None = None,
        max_reads: int = _SAMPLE_READS,
    ) -> QCResult:
        """Compute QC metrics for a FASTQ file.

        Args:
            fastq_file: Path to FASTQ file (.fastq, .fq, or .gz).
            platform:   Override platform detection ("illumina", "nanopore",
                        "pacbio", "iontorrent", "bgi"). If None, auto-detects.
            max_reads:  Maximum reads to process (default: 200,000).

        Returns:
            QCResult with all metrics populated.
        """
        path = Path(fastq_file)

        # Platform detection
        if platform:
            detected, confidence = platform, "user"
        else:
            detected, confidence = detect_platform(path)

        result = QCResult(str(path), detected, confidence)

        lengths:   list[int]   = []
        qualities: list[float] = []
        gc_counts: list[float] = []
        q20_bases = 0
        q30_bases = 0
        total_bases = 0
        first_seqs:  list[str] = []   # for duplicate estimation
        adapter_hits: Counter = Counter()

        for record in self._iter_records(path, max_reads):
            seq = str(record.seq).upper()
            ln  = len(seq)
            if ln == 0:
                continue

            lengths.append(ln)
            total_bases += ln

            # Quality
            quals = record.letter_annotations.get("phred_quality", [])
            if quals:
                mean_q = sum(quals) / len(quals)
                qualities.append(mean_q)
                q20_bases += sum(1 for q in quals if q >= 20)
                q30_bases += sum(1 for q in quals if q >= 30)

            # GC
            g = seq.count("G")
            c = seq.count("C")
            gc_counts.append((g + c) / ln * 100)

            # Adapter hints (check last 20 bp for Illumina)
            tail = seq[-20:] if ln >= 20 else seq
            for name, adapter in _ADAPTER_HINTS.items():
                if adapter[:_ADAPTER_CHECK_LEN] in tail:
                    adapter_hits[name] += 1

            # Duplicate estimation (store first 50 bp of first 50k reads)
            if len(first_seqs) < 50_000:
                first_seqs.append(seq[:50])

        n = len(lengths)
        if n == 0:
            return result

        result.total_reads   = n
        result.total_bases   = total_bases
        result.min_length    = min(lengths)
        result.max_length    = max(lengths)
        result.mean_length   = sum(lengths) / n
        result.median_length = float(sorted(lengths)[n // 2])
        result.n50           = self._n50(lengths)
        result.mean_quality  = mean(qualities) if qualities else 0.0
        result.pct_q20       = q20_bases / total_bases * 100 if total_bases else 0.0
        result.pct_q30       = q30_bases / total_bases * 100 if total_bases else 0.0
        result.gc_content    = sum(gc_counts) / n
        result.adapter_hits  = dict(adapter_hits)

        # Duplicate rate estimate
        unique = len(set(first_seqs))
        result.duplicate_rate = (1 - unique / len(first_seqs)) * 100 if first_seqs else 0.0

        # Length distribution (binned)
        bin_size = max(1, (result.max_length - result.min_length) // 20 or 1)
        dist: Counter = Counter()
        for ln in lengths:
            bucket = (ln // bin_size) * bin_size
            dist[bucket] += 1
        result.length_dist = dict(sorted(dist.items()))

        return result

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _iter_records(self, path: Path, max_reads: int) -> Iterator[SeqRecord]:
        opener = gzip.open if path.suffix == ".gz" else open
        fmt    = "fastq"
        with opener(path, "rt") as fh:
            count = 0
            for record in SeqIO.parse(fh, fmt):
                yield record
                count += 1
                if count >= max_reads:
                    break

    @staticmethod
    def _n50(lengths: list[int]) -> int:
        sorted_l = sorted(lengths, reverse=True)
        total    = sum(sorted_l)
        running  = 0
        for ln in sorted_l:
            running += ln
            if running >= total / 2:
                return ln
        return 0
