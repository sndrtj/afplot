"""
afplot.cli
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

import re

import click
import seaborn as sns
import vcf

from .utils import Region, get_contigs, bed_reader
from .whole_genome import histogram_main, scatter_main, distance_main
from .region import region_histogram_main, \
    region_scatter_main, region_distance_main


def validate_region_str(ctx, param, value):
    if value is None:
        return None
    region_re = re.compile('^([\w\d]+):(\d+)-(\d+)$')
    match = region_re.match(value)
    if match is not None:
        return Region(match.group(1),
                      int(match.group(2)),
                      int(match.group(3)))
    else:
        raise click.BadParameter('{0} is not a '
                                 'valid region string'.format(value))


shared_options_all = [
    click.option("--dpi",
                 type=int,
                 help="DPI for output PNGs (default: 300)",
                 default=300),
    click.option("--color-palette",
                 type=str,
                 help="The name of a color palette "
                      "to pass to seaborn.set_palette")
]


shared_options_regions = shared_options_all + [
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
    click.option("--name",
                 "-n",
                 type=str,
                 help="Optional title for plot"),
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


shared_options_genome = shared_options_all + [
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


def _setup_genome_values(**kwargs):
    """Setup values used for whole-genome plotting."""
    readers = [vcf.Reader(filename=x) for x in kwargs.get("vcf", [])]
    contigs = get_contigs(readers, kwargs.get("exclude-pattern", []))
    if len(kwargs.get('sample', [])) == 0:
        samples = [x.samples[0] for x in readers]
    else:
        samples = kwargs.get('sample', [])
    if kwargs.get('color-palette') is not None:
        if len(samples) == 1:
            sns.set_palette(kwargs.get('color-palette'), 4)
        else:
            sns.set_palette(kwargs.get('color-palette'), len(samples))
    return readers, contigs, samples


def _setup_region_values(**kwargs):
    """Setup values for region plotting."""
    reader = vcf.Reader(filename=kwargs.get("vcf"))
    region = kwargs.get("region")
    region_file = kwargs.get("region_file")
    margin = kwargs.get("margin", 0)
    if region is not None:
        nrs = [Region(region.chr, int(region.start)-margin,
                      int(region.end)+margin)]
    elif region_file is not None:
        nrs = bed_reader(region_file, margin)
    else:
        nrs = []
    return reader, nrs


@click.group(short_help="Whole-genome plots")
def cli_whole_genome(**kwargs):
    """
    Create whole-genome plots for one or multiple VCFs.

    If only one VCF is supplied, plots will be colored
    on call type (het/hom_ref/hom_alt).
    If multiple VCF files are supplied, plots will be colored per file/label.
    Only *one* sample per VCF file can be plotted.

    Your VCF file *MUST* contain an AD column in the FORMAT field.
    Your VCF file *MUST* have contig names and lengths placed in the header.
    Your VCF file *MUST* be indexed with tabix.

    VCF files preferably have the same contigs,
    i.e. they are produced with the same reference.
    If this is not the case, this script will select the vcf file with the
    largest number of contigs.

    You may exclude contigs by supplying a regex pattern to the -e parameter.
    This parameter may be repeated.
    """
    pass


@click.option("--kde-only", "-k", is_flag=True,
              help="Only show kernel density plot")
@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome histogram")
def whole_genome_histogram(**kwargs):
    """Create histograms over every chromosome."""
    readers, contigs, samples = _setup_genome_values(**kwargs)
    labels = kwargs.get('label', [])
    dpi = kwargs.get('dpi', None)
    kde = kwargs.get('kde-only', False)
    output = kwargs.get('output')
    if dpi is None:
        histogram_main(readers, labels, samples,
                       contigs, output, kde_only=kde)
    else:
        histogram_main(readers, labels, samples,
                       contigs, output, kde_only=kde, dpi=dpi)


@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome scatter plot")
def whole_genome_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every chromosome."""
    readers, contigs, samples = _setup_genome_values(**kwargs)
    labels = kwargs.get('label', [])
    dpi = kwargs.get('dpi', None)
    output = kwargs.get('output')
    if dpi is None:
        scatter_main(readers, labels, samples, contigs, output)
    else:
        scatter_main(readers, labels, samples, contigs, output, dpi=dpi)


@generic_option(shared_options_genome)
@click.command(short_help="Whole-genome distance plot")
def whole_genome_distance(**kwargs):
    """Create scatter plot distance to theoretical AF over very chromosome."""
    readers, contigs, samples = _setup_genome_values(**kwargs)
    labels = kwargs.get('label', [])
    dpi = kwargs.get('dpi', None)
    output = kwargs.get('output')
    if dpi is None:
        distance_main(readers, labels, samples, contigs, output)
    else:
        distance_main(readers, labels, samples, contigs, output, dpi=dpi)


@click.group(short_help="Region plots")
def cli_regions(**kwargs):
    """
    Create plots for regions of interest for one VCF.
    
    Plots will be colored on call type (het/hom_alt/hom_ref).
    
    Your VCF file *MUST* contain an AD column in the FORMAT field.
    Your VCF file *MUST* have contig names and lengths placed in the header.
    Your VCF file *MUST* be indexed with tabix.
    """
    pass


@click.option("--kde-only", "-k", is_flag=True,
              help="Only show kernel density plot")
@generic_option(shared_options_regions)
@click.command(short_help="Region histogram")
def region_histogram(**kwargs):
    """Create histograms of allele frequencies over every region."""
    reader, regions = _setup_region_values(**kwargs)
    region_histogram_main(
        reader,
        kwargs.get("output_dir"),
        regions,
        kwargs.get("name"),
        kwargs.get("dpi"),
        kwargs.get("kde-only")
    )


@generic_option(shared_options_regions)
@click.command(short_help="Region scatter plot")
def region_scatter(**kwargs):
    """Create scatter plot of allele frequencies over every region."""
    reader, regions = _setup_region_values(**kwargs)
    region_scatter_main(
        reader,
        kwargs.get("output_dir"),
        regions,
        kwargs.get("name"),
        kwargs.get("dpi")
    )


@generic_option(shared_options_regions)
@click.command(short_help="Region distance plot")
def region_distance(**kwargs):
    """Create scatter plot of distance to theoretical AF over every region."""
    reader, regions = _setup_region_values(**kwargs)
    region_distance_main(
        reader,
        kwargs.get("output_dir"),
        regions,
        kwargs.get("name"),
        kwargs.get("dpi")
    )


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
