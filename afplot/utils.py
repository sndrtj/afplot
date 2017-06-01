"""
afplot.utils
~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""

from collections import namedtuple
import vcf

Region = namedtuple("Region", ["chr", "start", "end"])


def _is_vcf_version_at_least_0_6_8():
    """
    The behaviour of vcfReader.fetch changed significantly
    from version 0.6.8 onwards
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
