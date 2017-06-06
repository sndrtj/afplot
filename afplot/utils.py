"""
afplot.utils
~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""


import re
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


def get_contigs(readers, exclude_patterns):
    """
    Get list of contigs to be used
    from a list of VCF readers and patters to exclude
    :param readers: list of VCF readers
    :param exclude_patterns: Regex patterns to exclude
    :return: list of usable contig names
    """
    contigs = get_longest_contig_list(readers)
    for pattern in exclude_patterns:
        regex = re.compile(pattern)
        contigs = [x for x in contigs if not regex.match(x)]
    return contigs


def region_key(r):
    """
    Create writeable name for a region
    :param r: region
    :return: str
    """
    return "{0}_{1}-{2}".format(r.chr, r.start, r.end)


def bed_reader(path, margin=0):
    """
    Generator of Regions for a bed file
    :param path: path to bed file
    :param margin: Margin to append to region
    :return: generator of Region
    """
    with open(path) as handle:
        for line in handle:
            s = line.strip().split("\t")
            yield Region(
                chr=s[0],
                start=int(s[1])-margin,
                end=int(s[2])+margin
            )
