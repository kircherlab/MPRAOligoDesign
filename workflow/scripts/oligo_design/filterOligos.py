#!/usr/bin/env python

import re
import pysam
import os
import pandas as pd
from optparse import OptionParser
import gzip
import numpy as np

LEFT = "AGGACCGGATCAACT"
RIGHT = "CATTGCGTGAACCGA"

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
    regions = gzip.open(regions_file, 'rt') #, names=["chrom", "start", "end", "id", "qfilter", "strand"], sep="\t", skiprows = 1)
    failed_list = []
    failed = False
    fail_reasons = {"TSS": 0, 
                "repeats": 0,
                "CTCF": 0}

    for region in regions:
        rchrom, rstart, rend, rid, _, _ = region.strip().split("\t")
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
    restriction_sites = zip(map(lambda site: re.compile(Site2RegEx(
       site), re.IGNORECASE), sites), map(lambda site: site.find("^"), sites))

    failed_list = []
    fail_reasons = {"restrictions": 0, "hompol": 0}

    if seqs[-3:] == ".gz":
        seqfile = gzip.open(seqs, 'rt')
    else:
        seqfile = open(seqs)

    fasta = fastaReader(seqfile)
    for cid, seq in fasta:
        cseq = LEFT + seq + RIGHT[:5]
        failed = False
        if (max_hom == None) or (nucleotideruns(cseq) <= max_hom):
            foundRestrictionSite = False
            for ind, (site, pos) in enumerate(restriction_sites):
                if any(site.finditer(cseq)):
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


def write_output(failed, map_file, remove_unused, out_file):

    map = pd.read_csv(map_file, sep="\t")

    for rid in failed["regions"]:
        if rid in map["Region"].values:
            map = map.drop(map[map["Region"] == rid].index)

    for sid in failed["seqs"]:
        if sid in map["REF_ID"].values:
            map = map.drop(map[map["REF_ID"] == sid].index)
        if remove_unused and sid in map["ALT_ID"].values:
            print(map[map["ALT_ID"] == sid])
            map = map.drop(map[map["ALT_ID"] == sid].index)

    map.to_csv(out_file, compression='gzip', sep="\t", index=False)


if __name__ == "__main__":
    parser = OptionParser("%prog [options]")
    parser.add_option("-s", "--seqs", dest="seqs",
                      help="File containing the designed oligos (def '')", default="results/oligo_design/test_1/design.fa")
    parser.add_option("-r", "--regions", dest="regions",
                      help="File containing the regions (def '')", default="results/oligo_design/test_1/design.regions.bed.gz")
    parser.add_option("-m", "--seq-map", dest="seq_map", 
                      help="Map that maps sequences to regions and variants", default="results/oligo_design/test_1/variant_region_map.tsv.gz")
    parser.add_option("-z", "--max_homopolymer_length", dest="maxHomLength",
                      help="Maximum homopolymer length (def 10)", default=10, type="int")  # deactivated before
    parser.add_option("-f", "--repeat", dest="repeat",
                      help="Maximum fraction explained by a single simple repeat annotation (def 0.25)", default=0.25, type="float")
    parser.add_option("--simple-repeats", dest="simple_repeats_file",
                      help="Bedfile for simple repeats", default="reference/simpleRepeat.bed.gz")
    parser.add_option("--tss-positions", dest="tss_pos_file",
                      help="Bedfile for TSS positions",default="reference/TSS_pos.bed.gz")
    parser.add_option("--ctcf-motifs", dest="ctcf_motif_file",
                      help="Bedfile for ctcf motifs", default="reference/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz")
    parser.add_option("-o", "--output-map", dest="map_out",
                      help="Output file of filtered sequence map", default="results/oligo_design/test_1/filtered.variant_region_map.tsv.gz")
    parser.add_option("-i", "--output-failed-variants", dest="variants_failed_out",
                      help="Output file of failed variant ids", default="results/oligo_design/test_1/failed.var.ids.txt")
    parser.add_option("-u", "--remove-unused-regions", dest="remove_unused",
                      help="Whether sequences without variants should be kept (def True)", default="true")                      
    (options, args) = parser.parse_args()

    repeatIndex = pysam.Tabixfile(options.simple_repeats_file)
    TSSIndex = pysam.Tabixfile(options.tss_pos_file)
    CTCFIndex = pysam.Tabixfile(options.ctcf_motif_file)

    remove_unused = options.remove_unused.lower() == "true"
 
    selectedSeqs = {}
    seqFilterMap = {}
    failedSeqs = []

    nr_failed = 0
    fail_reasons = {}
    failed = {}


    if options.regions and os.path.exists(options.regions):
        failed["regions"], fail_reasons = regions_filter(options.regions, repeatIndex, TSSIndex, CTCFIndex, options.repeat)

    if options.seqs and os.path.exists(options.seqs):
        failed["seqs"], reasons = seqs_filter(options.seqs, options.maxHomLength)
        fail_reasons |= reasons

    write_output(failed, options.seq_map, remove_unused, options.map_out)

    with open(options.variants_failed_out, 'wt') as out:
        out.writelines(failedSeqs)

    # print("""Failed sequences: 
    # %d due to homopolymers 
    # %d due to simple repeats 
    # %d due to TSS site overlap 
    # %d due to CTCF site overlap 
    # %d due to EcoRI or SbfI restriction site overlap""" % (homFailed, repeatsFailed, TSSfailed, CTCFfailed, restrictionsFailed))
    # print("Total failed: %d" % (nrFailed))
    # print("Total passed sequences: %d" % (len(selectedSeqs)))
