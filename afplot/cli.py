"""
afplot.cli
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

import re

import click

from .utils import Region


def validate_region_str(ctx, param, value):
    region_re = re.compile('^([\w\d]+):(\d+)-(\d+)$')
    match = region_re.match(value)
    if match is not None:
        return Region(match.group(1),
                      int(match.group(2)),
                      int(match.group(3)))
    else:
        raise click.BadParameter('{0} is not a '
                                 'valid region string'.format(value))


shared_options_regions = [
    click.option("--vcf",
                 "-v",
                 type=click.Path(exists=True),
                 required=True,
                 help="Path to input VCF file"),
    click.option("--output-dir",
                 "-o",
                 type=click.Path(exists=True),
                 required=True,
                 help="Path to output directory"),
    click.option("--region-file",
                 "-L",
                 type=click.Path(exists=True),
                 help="Path to region file"),
    click.option("--region",
                 "-R",
                 callback=validate_region_str,
                 help="Region string. Must be of format <contig:start-end>"),
    click.option("--margin",
                 "-m",
                 type=int,
                 help="Margin around regions to plot",
                 default=0)
]


shared_options_genome = [
    click.option("--vcf",
                 "-v",
                 type=click.Path(exists=True),
                 required=True,
                 multiple=True,
                 help="Path(s) to input VCF file(s)"),
    click.option("--label",
                 "-l",
                 type=str,
                 required=True,
                 multiple=True,
                 help="Label(s) to VCF file(s)"),
    click.option("--sample",
                 "-s",
                 type=str,
                 multiple=True,
                 help="Sample name(s) of VCF file(s). "
                      "If not given, will use fist sample in each VCF File"),
    click.option("--exclude-pattern",
                 "-e",
                 type=str,
                 multiple=True,
                 help="Regex pattern(s) to exclude from contig list"),
    click.option("--output",
                 "-o",
                 type=click.Path(exists=False),
                 required=True,
                 help="Path to output file")
]


def generic_option(options):
    """
    Decorator to add generic options to Click CLI's
    The group parent should NOT be decorated with this decorator
    :param options: list of click.option
    :return: decorated function
    """
    def __generic_option(func):
        for option in reversed(options):
            func = option(func)
        return func
    return __generic_option


@click.group(short_help="Whole-genome plots")
def cli_whole_genome(**kwargs):
    """Create whole-genome plots for one or multiple VCFs."""
    pass


@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome histogram")
def whole_genome_histogram(**kwargs):
    """Create histograms over every chromosome."""
    pass


@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome scatter plot")
def whole_genome_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every chromosome."""
    pass


@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome distance plot")
def whole_genome_distance(**kwargs):
    """Create scatter plot distance to theoretical AF over very chromosome."""
    pass


@click.group(short_help="Region plots")
def cli_regions(**kwargs):
    """Create plots for regions of interest for one VCF."""
    pass


@generic_option(shared_options_regions)
@click.command(short_help="Region histogram")
def region_histogram(**kwargs):
    """Create histograms of allele frequencies over every region."""
    pass


@generic_option(shared_options_regions)
@click.command(short_help="Region scatter plot")
def region_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every region."""
    pass


@generic_option(shared_options_regions)
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
