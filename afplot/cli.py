"""
afplot.cli
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

import click


@click.group(short_help="Whole-genome plots")
def cli_whole_genome(**kwargs):
    """Create whole-genome plots for one or multiple VCFs."""
    pass


@click.command(short_help="Whole-genome histogram")
def whole_genome_histogram(**kwargs):
    """Create histograms over every chromosome."""
    pass


@click.command(short_help="Whole-genome scatter plot")
def whole_genome_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every chromosome."""
    pass


@click.command(short_help="Whole-genome distance plot")
def whole_genome_distance(**kwargs):
    """Create scatter plot distance to theoretical AF over very chromosome."""
    pass


@click.group(short_help="Region plots")
def cli_regions(**kwargs):
    """Create plots for regions of interest for one VCF."""
    pass


@click.command(short_help="Region histogram")
def region_histogram(**kwargs):
    """Create histograms of allele frequencies over every region."""
    pass


@click.command(short_help="Region scatter plot")
def region_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every region."""
    pass


@click.command(short_help="Region distance plot")
def region_distance(**kwargs):
    """Create scatter plot ofdistance to theoretical AF over every region."""
    pass


@click.group()
def cli():
    """
    Plot allele frequencies in VCF files. 
    
    \b
    Two basic modes exist:
      - regions: Plot histogram, scatter or distance plots per
        user-specified region.
      - whole-genome: Plot histogram, scatter or distance plots over the
        entire genome.
    """
    pass


def main():
    cli_regions.add_command(region_histogram, "histogram")
    cli_regions.add_command(region_scatter, "scatter")
    cli_regions.add_command(region_distance, "distance")
    cli_whole_genome.add_command(whole_genome_histogram, "histogram")
    cli_whole_genome.add_command(whole_genome_scatter, "scatter")
    cli_whole_genome.add_command(whole_genome_distance, "distance")
    cli.add_command(cli_regions, "regions")
    cli.add_command(cli_whole_genome, "whole-genome")
    cli()

if __name__ == "__main__":
    main()




