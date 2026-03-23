"""
Alignment Engine.

Selects and runs the appropriate aligner for a given platform and experiment.

Level 1 (pip-only):  mappy — Python bindings for minimap2
Level 2 (external):  BWA-MEM2, STAR, minimap2 binary via subprocess

The engine auto-detects which level is available and uses the best option.

Usage:
    from bioreads.core.aligner import AlignmentEngine

    engine = AlignmentEngine.auto(platform="nanopore", experiment="dna")
    result = engine.align(
        reads="reads.fastq",
        reference="genome.fa",
        output="aligned.bam",
    )
    print(result.mapping_rate)
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from bioreads.core.detector import suggest_aligner


@dataclass
class AlignmentResult:
    """Statistics from an alignment run."""
    output_bam:    str
    aligner_used:  str
    level:         int
    total_reads:   int   = 0
    mapped_reads:  int   = 0
    mapping_rate:  float = 0.0
    multimappers:  int   = 0
    log:           str   = ""

    def summary(self) -> str:
        return (
            f"Aligner       : {self.aligner_used} (level {self.level})\n"
            f"Output        : {self.output_bam}\n"
            f"Total reads   : {self.total_reads:,}\n"
            f"Mapped reads  : {self.mapped_reads:,}\n"
            f"Mapping rate  : {self.mapping_rate:.1f}%\n"
            f"Multimappers  : {self.multimappers:,}\n"
        )


class AlignmentEngine:
    """Run alignment using the best available tool."""

    def __init__(self, platform: str, experiment: str = "dna"):
        self.platform   = platform
        self.experiment = experiment
        self._suggestion = suggest_aligner(platform, experiment)

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def auto(cls, platform: str, experiment: str = "dna") -> "AlignmentEngine":
        """Create an engine with automatic aligner selection.

        Args:
            platform:   "illumina", "nanopore", "pacbio", "iontorrent", "bgi".
            experiment: "dna", "rna", or "amplicon".

        Returns:
            Configured AlignmentEngine.
        """
        return cls(platform=platform, experiment=experiment)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    @property
    def available_level(self) -> int:
        """Return the highest operation level available (1 or 2)."""
        aligner = self._suggestion["aligner"]
        if shutil.which(aligner):
            return 2
        if self._suggestion.get("fallback") == "mappy":
            try:
                import mappy  # noqa: F401
                return 1
            except ImportError:
                pass
        return 0

    @property
    def aligner_info(self) -> dict:
        """Return the suggested aligner info for this platform/experiment."""
        return self._suggestion

    def align(
        self,
        reads: str | Path,
        reference: str | Path,
        output: str | Path = "aligned.bam",
        reads2: str | Path | None = None,
        threads: int = 4,
    ) -> AlignmentResult:
        """Align reads to a reference genome.

        Automatically selects the best available aligner:
        - Level 2: external binary (BWA-MEM2 / STAR / minimap2)
        - Level 1: mappy (Python minimap2 bindings)

        Args:
            reads:     Path to FASTQ file (R1 for paired-end).
            reference: Path to reference genome FASTA or index.
            output:    Output BAM file path.
            reads2:    Path to R2 FASTQ for paired-end (Illumina).
            threads:   Number of threads for the aligner.

        Returns:
            AlignmentResult with mapping statistics.
        """
        level = self.available_level
        if level == 0:
            raise RuntimeError(
                f"No aligner available for {self.platform}/{self.experiment}.\n"
                f"  Recommended: install {self._suggestion['aligner']} via bioconda or WSL2.\n"
                f"  Fallback:    pip install mappy (long reads only)."
            )
        if level == 2:
            return self._align_external(reads, reference, output, reads2, threads)
        return self._align_mappy(reads, reference, output)

    def check_installation(self) -> dict:
        """Check which aligners are installed and their versions.

        Returns:
            Dict mapping aligner name to version string or None if not found.
        """
        aligners = ["bwa-mem2", "bwa", "STAR", "minimap2"]
        result   = {}
        for a in aligners:
            path = shutil.which(a)
            if path:
                try:
                    out = subprocess.run(
                        [a, "--version"], capture_output=True, text=True, timeout=5
                    )
                    version = (out.stdout or out.stderr).split("\n")[0].strip()
                except Exception:
                    version = "installed"
                result[a] = version
            else:
                result[a] = None

        # mappy
        try:
            import mappy
            result["mappy"] = getattr(mappy, "__version__", "installed")
        except ImportError:
            result["mappy"] = None

        return result

    # ------------------------------------------------------------------
    # Private — external aligners
    # ------------------------------------------------------------------

    def _align_external(
        self,
        reads: str | Path,
        reference: str | Path,
        output: str | Path,
        reads2: str | Path | None,
        threads: int,
    ) -> AlignmentResult:
        aligner = self._suggestion["aligner"]
        preset  = self._suggestion.get("preset")
        out     = Path(output)

        if aligner == "minimap2":
            cmd = ["minimap2", "-a", "-t", str(threads)]
            if preset:
                cmd += ["-x", preset]
            cmd += [str(reference), str(reads)]
            if reads2:
                cmd.append(str(reads2))

        elif aligner == "bwa-mem2":
            cmd = ["bwa-mem2", "mem", "-t", str(threads), str(reference), str(reads)]
            if reads2:
                cmd.append(str(reads2))

        elif aligner == "STAR":
            cmd = [
                "STAR", "--runThreadN", str(threads),
                "--genomeDir", str(reference),
                "--readFilesIn", str(reads),
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outFileNamePrefix", str(out.parent / out.stem) + "_",
            ]
            if reads2:
                cmd.append(str(reads2))
        else:
            raise ValueError(f"Unknown aligner: {aligner}")

        log_lines = []
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600
            )
            log_lines = proc.stderr.splitlines()
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr[:500])
        except FileNotFoundError:
            raise RuntimeError(f"Aligner '{aligner}' not found in PATH.")

        result = AlignmentResult(
            output_bam=str(out),
            aligner_used=aligner,
            level=2,
            log="\n".join(log_lines),
        )
        self._parse_mapping_stats(log_lines, result, aligner)
        return result

    # ------------------------------------------------------------------
    # Private — mappy fallback (Level 1)
    # ------------------------------------------------------------------

    def _align_mappy(
        self,
        reads: str | Path,
        reference: str | Path,
        output: str | Path,
    ) -> AlignmentResult:
        import mappy
        from Bio import SeqIO

        preset = self._suggestion.get("preset", "map-ont")
        aligner = mappy.Aligner(str(reference), preset=preset)
        if not aligner:
            raise RuntimeError(f"mappy: could not load reference '{reference}'")

        total = mapped = 0
        out   = Path(output)

        with open(out, "w") as sam:
            # Write minimal SAM header
            sam.write("@HD\tVN:1.6\tSO:unsorted\n")
            sam.write(f"@PG\tID:mappy\tPN:mappy\n")

            for record in SeqIO.parse(str(reads), "fastq"):
                total += 1
                seq   = str(record.seq)
                hits  = list(aligner.map(seq))
                if hits:
                    mapped += 1
                    h = hits[0]
                    flag = 0 if h.strand == 1 else 16
                    cigar = h.cigar_str
                    sam.write(
                        f"{record.id}\t{flag}\t{h.ctg}\t{h.r_st+1}\t{h.mapq}"
                        f"\t{cigar}\t*\t0\t0\t{seq}\t*\n"
                    )
                else:
                    sam.write(
                        f"{record.id}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t*\n"
                    )

        rate = mapped / total * 100 if total else 0.0
        return AlignmentResult(
            output_bam=str(out),
            aligner_used="mappy",
            level=1,
            total_reads=total,
            mapped_reads=mapped,
            mapping_rate=rate,
        )

    # ------------------------------------------------------------------
    # Private — parse aligner logs
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_mapping_stats(lines: list[str], result: AlignmentResult, aligner: str):
        import re
        for line in lines:
            m = re.search(r"(\d+)\s+\+\s+\d+\s+in total", line)
            if m:
                result.total_reads = int(m.group(1))
            m = re.search(r"(\d+)\s+\+\s+\d+\s+mapped", line)
            if m:
                result.mapped_reads = int(m.group(1))
        if result.total_reads:
            result.mapping_rate = result.mapped_reads / result.total_reads * 100
