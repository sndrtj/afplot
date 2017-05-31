"""
afplot.cli
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

import click


@click.group()
def cli_whole_genome(**kwargs):
    """Create whole-genome plots for one or multiple VCFs."""
    pass


@cli_whole_genome.command()
def whole_genome_histogram(**kwargs):
    """Create histograms over every chromosome."""
    pass


@cli_whole_genome.command()
def whole_genome_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every chromosome."""
    pass


@cli_whole_genome.command()
def whole_genome_distance(**kwargs):
    """Create scatter plot distance to theoretical AF over very chromosome."""
    pass


@click.group()
def cli_regions(**kwargs):
    """Create plots for regions of interest for one VCF."""
    pass


@cli_regions.command()
def region_histogram(**kwargs):
    """Create histograms of allele frequencies over every region."""
    pass


@cli_regions.command()
def region_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every region."""
    pass


@cli_regions.command()
def region_distance(**kwargs):
    """Create scatter plot ofdistance to theoretical AF over every region."""
    pass


cli = click.CommandCollection(sources=[cli_whole_genome, cli_regions])

if __name__ == "__main__":
    cli()




