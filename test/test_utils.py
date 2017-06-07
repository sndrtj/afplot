"""
test_utils
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""
import re

from os.path import realpath, join, dirname
from unittest.mock import Mock

import vcf

from afplot.utils import region_key, Region, _is_vcf_version_at_least_0_6_8, \
    get_contigs, get_longest_contig_list, bed_reader

long = join(dirname(realpath(__file__)), "data/header_vcf/test.vcf")
short = join(dirname(realpath(__file__)), "data/header_vcf/test.autosomes.vcf")
mini_bed = join(dirname(realpath(__file__)), "data/mini.bed")


class TestUtils(object):

    def test_region_key(self):
        assert region_key(Region("chr1", 1, 100)) == "chr1_1-100"

    def test_is_vcf_0_6_8(self):
        mock_0_6_7 = Mock()
        mock_0_6_7.VERSION = "0.6.7"
        mock_0_6_8 = Mock()
        mock_0_6_8.VERSION = "0.6.8"
        mock_1 = Mock()
        mock_1.VERSION = "1.0.0"
        assert not _is_vcf_version_at_least_0_6_8(mock_0_6_7)
        assert _is_vcf_version_at_least_0_6_8(mock_0_6_8)
        assert _is_vcf_version_at_least_0_6_8(mock_1)

    def test_longest_contig(self):
        l_reader = vcf.Reader(filename=long)
        s_reader = vcf.Reader(filename=short)
        assert len(get_longest_contig_list([l_reader, s_reader])) == 84

    def test_get_contig(self):
        gl_re = re.compile("^.*gl.*$")
        l_reader = vcf.Reader(filename=long)
        s_reader = vcf.Reader(filename=short)
        assert sorted(get_contigs([l_reader, s_reader], [gl_re])) == sorted([
            "chrM", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
            "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
            "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
            "chr22", "chrX", "chrY"
        ])

    def test_bed_reader(self):
        no_margin_regions = [x for x in bed_reader(mini_bed)]
        assert no_margin_regions[0] == Region("chr1", 100, 200)
        assert no_margin_regions[1] == Region("chr1", 1000, 2000)
        assert no_margin_regions[2] == Region("chr1", 10000, 20000)

    def test_bed_reader_w_margins(self):
        margin_regions = [x for x in bed_reader(mini_bed, 500)]
        assert margin_regions[0] == Region("chr1", 0, 700)
        assert margin_regions[1] == Region("chr1", 500, 2500)
        assert margin_regions[2] == Region("chr1", 9500, 20500)
