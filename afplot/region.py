"""
afplot.region
~~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

from __future__ import print_function
from os.path import join
from warnings import warn

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .utils import NEW_VCF, region_key
from .variation import get_all_allele_freqs, \
    get_distance_to_exp, get_variant_type


def build_df_for_region(reader, region, sample=None, label=None):
    if label is None:
        label = "dummy"  # this is a hack, but FacetGrid won't work with None
    if sample is None:
        sample = reader.samples[0]
    if NEW_VCF:
        iterator = reader.fetch(
            region.chr,
            int(region.start),
            int(region.end)
        )
    else:
        iterator = reader.fetch(
            region.chr,
            int(region.start)+1,
            int(region.end)
        )
    maf = []
    for record in iterator:
        ad = get_all_allele_freqs(record, sample)
        distances = get_distance_to_exp(record, sample)
        if len(ad) == 0:
            continue
        for freq, dist in zip(ad, distances):
                maf.append([record.POS, freq,
                            get_variant_type(record, sample), dist])
    arr = np.array(maf)
    if len(arr) == 0:
        return None
    df = pd.DataFrame(
        {"pos": [int(x) for x in arr[:, 0]],
         "af": [float(x) for x in arr[:, 1]],
         "label": arr[:, 2],
         "distance": [float(x) for x in arr[:, 3]],
         "chrom": [label]*len(arr[:, 0])
         }
    )
    return df


def plot_single_histogram(dataframe, output, dpi=300,
                          kde_only=False, label=None):
    g = sns.FacetGrid(dataframe, col="chrom", hue="label", col_wrap=2)
    if kde_only:
        g = (g.map(sns.distplot, "af", hist=False).add_legend().set_titles(""))
    else:
        g = (g.map(sns.distplot, "af").add_legend().set_titles(""))
    for x in g.axes:
        if x.get_ylim()[1] > 10:
            x.set_ylim(0, 10)
        x.set_xlim(-0.5, 1.5)
    if label is not None:
        plt.title(label)
    plt.savefig(output, dpi=dpi)
    for x in g.axes:
        plt.sca(x)
    plt.close(g.fig)


def plot_single_scatter(dataframe, output, category="af", dpi=300, label=None):
    f = sns.lmplot("pos", category, dataframe, col="chrom",
                   col_wrap=1.2, fit_reg=False,
                   hue="label", scatter_kws={"alpha": 0.3}, aspect=3)
    f.add_legend()
    f.set_titles("")
    for x in f.axes:
        x.set_ylim(0, 1.0)
    if label is not None:
        plt.title(label)
    plt.savefig(output, dpi=dpi)
    for x in f.axes:
        plt.sca(x)
    plt.close(f.fig)


def region_histogram_main(reader, output_dir, regions,
                          label, dpi=300, kde_only=False):
    for reg in regions:
        name = region_key(reg)
        opath = join(output_dir, "{0}.png".format(name))
        df = build_df_for_region(reader, reg, label=label)
        if df is None:
            warn("Region {0} is empty".format(name))
            continue
        plot_single_histogram(df, opath, dpi, kde_only, label=label)


def region_scatter_main(reader, output_dir, regions, label, dpi=300):
    for reg in regions:
        name = region_key(reg)
        opath = join(output_dir, "{0}.png".format(name))
        df = build_df_for_region(reader, reg, label=label)
        if df is None:
            warn("Region {0} is empty".format(name))
            continue
        plot_single_scatter(df, opath, "af", dpi=dpi, label=label)


def region_distance_main(reader, output_dir, regions, label, dpi=300):
    for reg in regions:
        name = region_key(reg)
        opath = join(output_dir, "{0}.png".format(name))
        df = build_df_for_region(reader, reg, label=label)
        if df is None:
            warn("Region {0} is empty".format(name))
            continue
        plot_single_scatter(df, opath, "distance", dpi=dpi, label=label)
