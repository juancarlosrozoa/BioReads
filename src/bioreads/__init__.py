"""BioReads — NGS quality control and alignment for multiple sequencing platforms."""

from bioreads.core.qc import QCEngine
from bioreads.core.aligner import AlignmentEngine

__all__ = ["QCEngine", "AlignmentEngine"]
__version__ = "0.1.0"
