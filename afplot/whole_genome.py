"""
afplot.whole_genome
~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

from __future__ import print_function
import sys
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import progressbar
import seaborn as sns

from .utils import NEW_VCF
from .variation import get_all_allele_freqs, \
    get_distance_to_exp, get_variant_type


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
                    maf.append([record.POS, freq,
                                get_variant_type(record, sample), dist])
                else:
                    maf.append([record.POS, freq, label, dist])
    return np.array(maf)


def build_dataframe(readers, labels, samples, contigs):
    assert len(readers) == len(labels) and len(readers) == len(samples)
    the_dict = OrderedDict()
    for r, l, s in zip(readers, labels, samples):
        sample_dict = OrderedDict()
        for chrom in contigs:
            message = "Processing chromosome {0} " \
                      "for sample {1}".format(chrom, s)
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
    Clean dataframe so that it removes categories
    where all values of column are 0
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


def histogram_main(readers, labels, samples, contigs,
                   png, dpi=300, kde_only=False):
    df = build_dataframe(readers, labels, samples, contigs)
    df = clean_df(df, contigs)
    g = sns.FacetGrid(df, col="chromosome", hue="label",
                      aspect=1, col_wrap=4, sharey=False)
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
