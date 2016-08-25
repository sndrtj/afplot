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
import sys
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import progressbar
import seaborn as sns
import vcf


def _is_vcf_version_at_least_0_6_8():
    """
    The behaviour of vcfReader.fetch changed significantly from version 0.6.8 onwards
    :return: boolean
    """
    major, minor, patch = vcf.VERSION.split(".")
    if int(major) == 0 and int(minor) == 6 and int(patch) >= 8:
        return True
    elif int(major) == 0 and int(minor) > 6:
        return True
    elif int(major) > 0:
        return True
    return False


NEW_VCF = _is_vcf_version_at_least_0_6_8()


def get_longest_contig_list(readers):
    """
    Get the largest list of contig names from a list of readers
    :param readers: list of vcf readers
    :return: list of contig names
    """
    sorted_readers = sorted(readers, key=lambda x: len(x.contigs))
    return sorted_readers[0].contigs.keys()


def get_all_allele_freqs(record, sample_name):
    fmt = record.genotype(sample_name)
    if not hasattr(fmt.data, 'AD'):
        return []
    ad = fmt.data.AD
    if not ad:
        return []
    if len(ad) == 0:
        return []
    if sum(ad) == 0:
        return [0.0 for _ in ad]
    return [float(x)/sum(ad) for x in ad]


def get_variant_type(record, sample_name):
    fmt = record.genotype(sample_name)
    if not fmt.called:
        return "no_call"
    elif hasattr(fmt.data, "GQ") and fmt.data.GQ == 0:
        return "no_call"
    elif not fmt.is_variant:
        return "hom_ref"
    elif fmt.is_het:
        return "het"
    else:
        return "hom_alt"


def get_distance_to_exp(record, sample_name):
    """
    Get distance to expected theoretical allele frequencies
    This assumes the AD to field to conform to GATK spec
    i.e. the number of AD values EXACTLY matches the number of alleles (INCLUDING ref allele)
    :param record: VCF record
    :param sample_name: sample name
    :return: list of distances
    """
    freqs = get_all_allele_freqs(record, sample_name)
    if len(freqs) == 0:
        return []
    rtype = get_variant_type(record, sample_name)
    fmt = record.genotype(sample_name)
    assert len(freqs) == len(record.alleles)
    if rtype == "no_call":
        return [0 for _ in freqs]
    elif rtype == "hom_ref":
        # freq of ref allele should be 1.0, freq of all other alleles should be 0.0
        return [1 - freqs[0]] + freqs[1:]
    elif rtype == "hom_alt":
        # affected allele should be 1.0, all other alleles should be 0.0
        if fmt.phased:
            idx_affected_allele = int(fmt.data.GT.split("|")[0])
        else:
            idx_affected_allele = int(fmt.data.GT.split("/")[0])
        distances = []
        for i, f in enumerate(freqs):
            if i == idx_affected_allele:
                distances.append(1 - f)
            else:
                distances.append(f)
        return distances
    elif rtype == "het":
        # freq of affected alleles should be 0.5, all other alleles should be 0.0
        if fmt.phased:
            idx_affected_alleles = [int(x) for x in fmt.data.GT.split("|")]
        else:
            idx_affected_alleles = [int(x) for x in fmt.data.GT.split("/")]
        distances = []
        for i, f in enumerate(freqs):
            if i in idx_affected_alleles:
                distances.append(abs(0.5 - f))
            else:
                distances.append(f)
        return distances
    else:
        raise NotImplementedError


def get_array_for_chrom_all(reader, chromosome, label=None, sample=None):
    """
    Get MAF array for a contig from a reader
    :param reader: vcf reader object (must be tabixxed)
    :param chromosome: contig name
    :return: 4d-array of POS:AF:TYPE:DISTANCE
    """

    maf = []
    l = reader.contigs.get(chromosome).length
    if not sample:
        sample = reader.samples[0]
    with progressbar.ProgressBar(max_value=l, redirect_stdout=True) as bar:
        try:
            if NEW_VCF:
                iterator = reader.fetch(chromosome, 0)
            else:
                iterator = reader.fetch(chromosome, 1, l)
        except ValueError:
            return np.array(maf)
        for record in iterator:
            bar.update(record.POS)
            ad = get_all_allele_freqs(record, sample)
            distances = get_distance_to_exp(record, sample)
            if len(ad) == 0:
                continue
            for freq, dist in zip(ad, distances):
                if not label:
                    maf.append([record.POS, freq, get_variant_type(record, sample), dist])
                else:
                    maf.append([record.POS, freq, label, dist])
    return np.array(maf)


def build_dataframe(readers, labels, samples, contigs):
    assert len(readers) == len(labels) and len(readers) == len(samples)
    the_dict = OrderedDict()
    for r, l, s in zip(readers, labels, samples):
        sample_dict = OrderedDict()
        for chrom in contigs:
            message = "Processing chromosome {0} for sample {1}".format(chrom, s)
            print(message, file=sys.stderr)
            if len(readers) == 1:
               arr = get_array_for_chrom_all(r, chrom)
            else:
                arr = get_array_for_chrom_all(r, chrom, l, s)
            message = "{0} data points processed".format(len(arr))
            print(message, file=sys.stderr)
            if len(arr) == 0:
                continue
            df = pd.DataFrame(
                {"pos": [int(x) for x in arr[:, 0]],
                 "af": [float(x) for x in arr[:, 1]],
                 "label": arr[:, 2],
                 "distance": [float(x) for x in arr[:, 3]],
                 "chromosome": [chrom for _ in range(len(arr[:, 2]))]
                 }
            )
            sample_dict[chrom] = df
        sample_df = pd.concat(sample_dict.values())
        the_dict[s] = sample_df
    return pd.concat(the_dict.values())


def clean_df(df, contigs, column="af"):
    """
    Clean dataframe so that it removes categories where all values of column are 0
    :param df:
    :return: cleaned df
    """
    labels = set(df.label)
    tmp_dfs = []
    for chrom in contigs:
        t = df[df.chromosome == chrom]
        for l in labels:
            tt = t[t.label == l]
            col = getattr(tt, column)
            if not all((x == 0 for x in col)):
                tmp_dfs.append(tt)
    return pd.concat(tmp_dfs)


def scatter_main(readers, labels, samples, contigs, png, dpi=300):
    df = build_dataframe(readers, labels, samples, contigs)
    f = sns.lmplot("pos", "af", df, col="chromosome",
                   col_wrap=4, fit_reg=False,
                   hue="label", scatter_kws={"alpha": 0.3}, aspect=3)

    for i, x in enumerate(f.axes):
        x.set_xlim(0, )
    plt.savefig(png, dpi=dpi)


def histogram_main(readers, labels, samples, contigs, png, dpi=300, kde_only=False):
    df = build_dataframe(readers, labels, samples, contigs)
    df = clean_df(df, contigs)
    g = sns.FacetGrid(df, col="chromosome", hue="label", aspect=1, col_wrap=4, sharey=False)
    if kde_only:
        g = (g.map(sns.distplot, "af", hist=False).add_legend())
    else:
        g = (g.map(sns.distplot, "af").add_legend())
    for x in g.axes:
        if x.get_ylim()[1] > 10:
            x.set_ylim(0, 10)
        x.set_xlim(-0.5, 1.5)
    plt.savefig(png, dpi=dpi)


def distance_main(readers, labels, samples, contigs, png, dpi=300):
    df = build_dataframe(readers, labels, samples, contigs)
    f = sns.lmplot("pos", "distance", df, col="chromosome",
                   col_wrap=4, fit_reg=False,
                   hue="label", scatter_kws={"alpha": 0.3}, aspect=3)
    for i, x in enumerate(f.axes):
        x.set_xlim(0, )
        x.set_ylim(0, 0.5)
    plt.savefig(png, dpi=dpi)


def main():
    desc = """
    Create scatter plots or histogram of allele frequencies in vcf files.
    If only one VCF is supplied, plots will be colored on call type (het/hom_ref/hom_alt).
    If multiple VCF files are supplied, plots will be colored per file/label.
    Only *one* sample per VCF file can be plotted.

    Your VCF file *MUST* contain an AD column in the FORMAT field.
    Your VCF file *MUST* have contig names and lengths placed in the header.
    Your VCF file *MUST* be indexed with tabix.

    VCF files preferably have the same contigs,
    i.e. they are produced with the same reference.
    If this is not the case, this script will select the vcf file with the largest number of contigs.

    You may exclude contigs by supplying a regex pattern to the -e parameter.
    This parameter may be repeated.
    """
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--vcf", type=str,
                        help="Input vcf file(s)",
                        required=True,
                        action="append")
    parser.add_argument("-l", "--label", type=str,
                        help="Labels to vcf file(s)",
                        required=True,
                        action="append")
    parser.add_argument("-s", "--sample", type=str,
                        help="Sample identifiers (1 per vcf). Uses first sample in vcf by default",
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
                       help="Make scatter plot of AFs per chromosome", action="store_true")
    group.add_argument("--histogram",
                       help="Make histogram of AFs per chromosome",
                       action="store_true")
    group.add_argument("--distance",
                       help="Create scatter plot of distances to expected AFs",
                       action="store_true")

    parser.add_argument("-e", "--exclude-pattern",
                        type=str,
                        action="append",
                        help="Regex pattern to exclude from contig list")

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

    if args.scatter:
        scatter_main(readers, args.label, samples, contigs, args.output, args.dpi)
    elif args.histogram:
        histogram_main(readers, args.label, samples, contigs, args.output, args.dpi, args.kde_only)
    elif args.distance:
        distance_main(readers, args.label, samples, contigs, args.output, args.dpi)


if __name__ == '__main__':
    main()
