#!/usr/bin/env python

import re
import gzip

translation = {"B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "K": "[GT]",
               "M": "[AC]", "N": ".", "R": "[AG]", "S": "[CG]",
               "V": "[ACG]", "W": "[AT]", "Y": "[CT]",
               "A": "A", "C": "C", "G": "G", "T": "T"}


def fastaReader(filestream):
    name, value = None, None
    for line in filestream:
        if line.startswith(">"):
            if value != None:
                yield name, value
            value = None
            name = line[1:].strip()
        else:
            value = line.strip()
    if value != None:
        yield name, value


def nucleotideruns(seq):
    longestrun = 0
    lbase, llen = 'N', 0
    for base in seq:
        if base == lbase:
            llen += 1
        else:
            if llen > longestrun:
                longestrun = llen
            llen = 1
            lbase = base
    return longestrun


def Site2RegEx(seq):
    reqEx = ""
    for elem in seq:
        if elem in translation:
            reqEx += translation[elem]
    return reqEx


def regions_filter(regions_file, repeatIndex, TSSIndex, CTCFIndex, max_repeats):
    # , names=["chrom", "start", "end", "id", "qfilter", "strand"], sep="\t", skiprows = 1)
    regions = gzip.open(regions_file, 'rt')
    failed_list = []
    fail_reasons = {"TSS": 0,
                    "repeats": 0,
                    "CTCF": 0}

    for region in regions:
        failed = False
        region_split = region.strip().split("\t")
        rchrom, rstart, rend, rid = region_split[0:4]
        rstart, rend = int(rstart), int(rend)
        # filter simple repeats
        for line in repeatIndex.fetch(rchrom, rstart, rend):
            fields = line.split("\t")
            tstart, tend = int(fields[1]), int(fields[2])
            if (min(tend, rend)-max(tstart, rstart))/float(rend-rstart) > max_repeats:
                failed = True
                fail_reasons["repeats"] += 1
                break

        # filter TSS
        if any(TSSIndex.fetch(rchrom, rstart, rend)):
            failed = True
            fail_reasons["TSS"] += 1

        # filter CTCF
        if any(CTCFIndex.fetch(rchrom, rstart, rend)):
            failed = True
            fail_reasons["CTCF"] += 1

        if failed:
            failed_list.append(rid)

    return failed_list, fail_reasons


def seqs_filter(seqs, max_hom):
    names = {}
    sites = []
    names["CCTGCA^GG"] = "SbfI"
    sites.append("CCTGCA^GG")
    names["G^AATTC"] = "EcoRI"
    sites.append("G^AATTC")
    restriction_sites = list(zip(map(lambda site: re.compile(Site2RegEx(
        site), re.IGNORECASE), sites), map(lambda site: site.find("^"), sites)))
    
    failed_list = []
    fail_reasons = {"restrictions": 0, "hompol": 0}

    if seqs[-3:] == ".gz":
        seqfile = gzip.open(seqs, 'rt')
    else:
        seqfile = open(seqs)

    fasta = fastaReader(seqfile)
    for cid, seq in fasta:
        if (max_hom == None) or (nucleotideruns(seq) <= max_hom):
            foundRestrictionSite = False
            for ind, (site, pos) in enumerate(restriction_sites):
                if any(site.finditer(seq)):
                    foundRestrictionSite = True
                    break
            if foundRestrictionSite:
                fail_reasons["restrictions"] += 1
                failed_list.append(cid)
        else:
            fail_reasons["hompol"] += 1
            failed_list.append(cid)

    seqfile.close()
    return failed_list, fail_reasons
