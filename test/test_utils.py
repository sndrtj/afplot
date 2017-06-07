"""
test_utils
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""
from unittest.mock import Mock

from afplot.utils import region_key, Region, _is_vcf_version_at_least_0_6_8, \
    get_contigs, get_longest_contig_list, bed_reader


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





