# BioReads

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Platforms](https://img.shields.io/badge/Platforms-Illumina%20|%20Nanopore%20|%20PacBio%20|%20IonTorrent-green)](#supported-platforms)

**BioReads** is an open-source Python tool for NGS quality control and alignment across multiple sequencing platforms.

It automatically detects the sequencing platform from your data, selects the appropriate aligner, and generates a detailed QC report — without requiring you to learn a different tool for each platform.

### Supported platforms

| Platform | Type | Aligner |
|---|---|---|
| Illumina (HiSeq, NovaSeq, MiSeq) | Short reads | BWA-MEM2 / STAR |
| Oxford Nanopore (ONT) | Long reads | minimap2 / mappy |
| PacBio HiFi (CCS) | Long reads | minimap2 / mappy |
| Ion Torrent | Short reads | BWA-MEM2 |
| BGI / DNBSEQ | Short reads | BWA-MEM2 |

### Two operation levels

| Level | Requirements | Capabilities |
|---|---|---|
| **Level 1** — pip only | `mappy`, `biopython` | FASTQ QC, lightweight alignment via mappy |
| **Level 2** — bioconda / WSL2 | BWA-MEM2, STAR, minimap2 | Full alignment, splice-aware RNA-seq, BAM output |

BioReads detects which level is available and works with what it finds.

### Three interfaces

```bash
# GUI
bioreads gui

# CLI
bioreads qc reads.fastq
bioreads align --platform illumina --experiment rna reads_R1.fastq reads_R2.fastq --ref genome.fa

# Python library
from bioreads import QCEngine, AlignmentEngine
```
