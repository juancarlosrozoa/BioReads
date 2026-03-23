"""BioReads CLI."""

import click
from pathlib import Path


@click.group()
@click.version_option(package_name="bioreads")
def cli():
    """BioReads — NGS quality control and alignment."""


@cli.command("qc")
@click.argument("fastq", type=click.Path(exists=True, dir_okay=False))
@click.option("--platform", "-p", default=None,
              type=click.Choice(["illumina", "nanopore", "pacbio", "iontorrent", "bgi"]),
              help="Sequencing platform (auto-detected if omitted).")
@click.option("--max-reads", default=200_000, show_default=True,
              help="Maximum reads to process.")
@click.option("--output", "-o", default=None, help="Save text report to file.")
def qc(fastq, platform, max_reads, output):
    """Run quality control on a FASTQ file."""
    from bioreads.core.qc import QCEngine
    result = QCEngine().run(fastq, platform=platform, max_reads=max_reads)
    click.echo(result.summary())
    if output:
        Path(output).write_text(result.summary(), encoding="utf-8")
        click.echo(f"\nReport saved to {output}")


@cli.command("align")
@click.argument("reads", type=click.Path(exists=True, dir_okay=False))
@click.argument("reference", type=click.Path(exists=True, dir_okay=False))
@click.option("--platform", "-p", required=True,
              type=click.Choice(["illumina", "nanopore", "pacbio", "iontorrent", "bgi"]),
              help="Sequencing platform.")
@click.option("--experiment", "-e", default="dna", show_default=True,
              type=click.Choice(["dna", "rna", "amplicon"]),
              help="Experiment type.")
@click.option("--reads2", default=None, type=click.Path(exists=True),
              help="R2 FASTQ for paired-end Illumina.")
@click.option("--output", "-o", default="aligned.bam", show_default=True,
              help="Output BAM/SAM file.")
@click.option("--threads", "-t", default=4, show_default=True,
              help="Number of threads.")
def align(reads, reference, platform, experiment, reads2, output, threads):
    """Align reads to a reference genome."""
    from bioreads.core.aligner import AlignmentEngine
    engine = AlignmentEngine.auto(platform=platform, experiment=experiment)

    info   = engine.aligner_info
    level  = engine.available_level

    click.echo(f"Platform    : {platform} / {experiment}")
    click.echo(f"Aligner     : {info['aligner']}"
               + (f" -x {info['preset']}" if info.get('preset') else ""))
    click.echo(f"Level       : {level} ({'external binary' if level == 2 else 'mappy' if level == 1 else 'NOT AVAILABLE'})")

    if level == 0:
        click.echo(f"\nNo aligner available. Install {info['aligner']} via bioconda or WSL2.", err=True)
        raise SystemExit(1)

    result = engine.align(reads, reference, output=output, reads2=reads2, threads=threads)
    click.echo("\n" + result.summary())


@cli.command("check-tools")
def check_tools():
    """Check which aligners are installed."""
    from bioreads.core.aligner import AlignmentEngine
    engine = AlignmentEngine("illumina")
    tools  = engine.check_installation()
    click.echo("Aligner availability:\n")
    for tool, version in tools.items():
        status = f"OK   {version}" if version else "--   not found"
        click.echo(f"  {tool:<12} {status}")


@cli.command("detect")
@click.argument("fastq", type=click.Path(exists=True, dir_okay=False))
def detect(fastq):
    """Detect the sequencing platform of a FASTQ file."""
    from bioreads.core.detector import detect_platform, suggest_aligner
    platform, confidence = detect_platform(fastq)
    suggestion = suggest_aligner(platform)
    click.echo(f"Platform   : {platform}")
    click.echo(f"Confidence : {confidence}")
    click.echo(f"Aligner    : {suggestion['aligner']}"
               + (f" -x {suggestion['preset']}" if suggestion.get('preset') else ""))
    click.echo(f"Level      : {suggestion['level']}")
    if suggestion.get("fallback"):
        click.echo(f"Fallback   : {suggestion['fallback']} (pip install mappy)")


@cli.command("gui")
def launch_gui():
    """Launch the graphical interface."""
    from bioreads.gui.app import main
    main()
