"""
afplot.variation
~~~~~~~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""


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
    i.e. the number of AD values EXACTLY matches the number
    of alleles (INCLUDING ref allele)
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
        # freq of ref allele should be 1.0,
        # freq of all other alleles should be 0.0
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
        # freq of affected alleles should be 0.5,
        # all other alleles should be 0.0
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
