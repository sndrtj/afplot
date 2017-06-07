"""
test_variation
~~~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

from os.path import realpath, join, dirname

import pytest
import vcf

from afplot.variation import get_variant_type, get_all_allele_freqs, get_distance_to_exp

mini_vcf = join(dirname(realpath(__file__)), "data/mini.vcf")


@pytest.fixture
def vcf_records():
    reader = vcf.Reader(filename=mini_vcf)
    return [x for x in reader]


class TestVariation(object):

    def test_get_variant_type(self, vcf_records):
        assert get_variant_type(vcf_records[0], "SAMPLE1") == "hom_alt"
        assert get_variant_type(vcf_records[1], "SAMPLE1") == "hom_ref"
        assert get_variant_type(vcf_records[2], "SAMPLE1") == "het"
        assert get_variant_type(vcf_records[3], "SAMPLE1") == "no_call"
        assert get_variant_type(vcf_records[4], "SAMPLE1") == "no_call"

    def test_get_allele_freqs(self, vcf_records):
        assert get_all_allele_freqs(vcf_records[0], "SAMPLE1") == [(1.0/101), (100.0/101)]
        assert get_all_allele_freqs(vcf_records[1], "SAMPLE1") == [1.0, 0.0]
        assert get_all_allele_freqs(vcf_records[2], "SAMPLE1") == [(50.0/85), (35.0/85)]
        assert get_all_allele_freqs(vcf_records[3], "SAMPLE1") == [0.5, 0.5]
        assert get_all_allele_freqs(vcf_records[4], "SAMPLE1") == []

    def test_get_distances(self, vcf_records):
        assert get_distance_to_exp(vcf_records[0], "SAMPLE1") == [(1.0/101), 1-(100.0/101)]
        assert get_distance_to_exp(vcf_records[1], "SAMPLE1") == [0, 0]
        assert get_distance_to_exp(vcf_records[2], "SAMPLE1") == [(50.0/85)-0.5, 0.5-(35.0/85)]
        assert get_distance_to_exp(vcf_records[3], "SAMPLE1") == [0, 0]
        assert get_distance_to_exp(vcf_records[4], "SAMPLE1") == []
