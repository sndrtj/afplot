"""
afplot.afplot
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""

from __future__ import print_function
import argparse
import re
import seaborn as sns
import vcf

from .utils import _is_vcf_version_at_least_0_6_8, \
    NEW_VCF, get_longest_contig_list
from .variation import get_all_allele_freqs, \
    get_variant_type, get_distance_to_exp
from .whole_genome import get_array_for_chrom_all, build_dataframe, \
    clean_df, scatter_main, histogram_main, distance_main


def main():
    desc = """
    Create scatter plots or histogram of allele frequencies in vcf files.
    If only one VCF is supplied, plots will be colored
    on call type (het/hom_ref/hom_alt).
    If multiple VCF files are supplied, plots will be colored per file/label.
    Only *one* sample per VCF file can be plotted.

    Your VCF file *MUST* contain an AD column in the FORMAT field.
    Your VCF file *MUST* have contig names and lengths placed in the header.
    Your VCF file *MUST* be indexed with tabix.

    VCF files preferably have the same contigs,
    i.e. they are produced with the same reference.
    If this is not the case, this script will select the vcf file
    with the largest number of contigs.

    You may exclude contigs by supplying a regex pattern to the -e parameter.
    This parameter may be repeated.
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--vcf", type=str,
                        help="Input vcf file(s)",
                        required=True,
                        action="append")
    parser.add_argument("-l", "--label", type=str,
                        help="Labels to vcf file(s)",
                        required=True,
                        action="append")
    parser.add_argument("-s", "--sample", type=str,
                        help="Sample identifiers (1 per vcf). "
                             "Uses first sample in vcf by default",
                        action="append")
    parser.add_argument("-o", "--output",
                        type=str,
                        help="Path to output png",
                        required=True)
    parser.add_argument("--dpi", type=int,
                        default=300,
                        help="DPI for output png (default: 300)")
    parser.add_argument("-k", "--kde-only",
                        help="Only show kernel density plot on histogram",
                        action="store_true")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--scatter",
                       help="Make scatter plot of AFs per chromosome",
                       action="store_true")
    group.add_argument("--histogram",
                       help="Make histogram of AFs per chromosome",
                       action="store_true")
    group.add_argument("--distance",
                       help="Create scatter plot of distances to expected AFs",
                       action="store_true")

    parser.add_argument("-e", "--exclude-pattern",
                        type=str,
                        action="append",
                        default=[],
                        help="Regex pattern to exclude from contig list")
    parser.add_argument("--color-palette",
                        type=str,
                        help="The name of a color palette to pass to "
                             "seaborn.set_palette")

    args = parser.parse_args()

    readers = [vcf.Reader(filename=x) for x in args.vcf]
    if not args.sample:
        samples = [x.samples[0] for x in readers]
    else:
        samples = args.sample

    contigs = get_longest_contig_list(readers)
    for pattern in args.exclude_pattern:
        regex = re.compile(pattern)
        contigs = [x for x in contigs if not regex.match(x)]

    if args.color_palette:
        if len(samples) == 1:
            sns.set_palette(args.color_palette, 4)
        else:
            sns.set_palette(args.color_palette, len(samples))

    if args.scatter:
        scatter_main(readers, args.label, samples,
                     contigs, args.output, args.dpi)
    elif args.histogram:
        histogram_main(readers, args.label, samples,
                       contigs, args.output, args.dpi, args.kde_only)
    elif args.distance:
        distance_main(readers, args.label, samples,
                      contigs, args.output, args.dpi)


if __name__ == '__main__':
    main()
