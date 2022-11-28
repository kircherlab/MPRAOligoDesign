#!/usr/bin/env python

import re
import pysam
import os
from optparse import OptionParser
import gzip

LEFT = "AGGACCGGATCAACT"
MIDDLE = "CCTGCAGGGAATTC"
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


if __name__ == "__main__":
    parser = OptionParser("%prog [options]")
    parser.add_option("-l", "--oligos", dest="oligos",
                      help="File containing the designed oligos (def '')", default="data/mEB_IGVF_MPRA_array_Y1_all_55k_20221012.fa.gz")
    parser.add_option("-t", "--total", dest="total",
                      help="Total number of oligos to design (def 300000)", default=300000, type="int")
    parser.add_option("-z", "--max_homopolymer_length", dest="maxHomLength",
                      help="Maximum homopolymer length (def 10)", default=10, type="int")  # deactivated before
    parser.add_option("-f", "--repeat", dest="repeat",
                      help="Maximum fraction explained by a single simple repeat annotation (def 0.25)", default=0.25, type="float")
    parser.add_option("-o", "--outfile", dest="outfile",
                      help="Output file (def 'design.tsv')", default="design.tsv")
    (options, args) = parser.parse_args()

    repeatIndex = pysam.Tabixfile("reference/simpleRepeat.bed.gz")
    TSSIndex = pysam.Tabixfile("reference/TSS_pos.bed.gz")
    CTCFIndex = pysam.Tabixfile(
        "reference/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz")

    names = {}
    sites = []
    names["CCTGCA^GG"] = "SbfI"
    sites.append("CCTGCA^GG")
    names["G^AATTC"] = "EcoRI"
    sites.append("G^AATTC")

    restriction_sites = zip(map(lambda site: re.compile(Site2RegEx(
        site), re.IGNORECASE), sites), map(lambda site: site.find("^"), sites))
    selectedSeqs = {}

    TSSfailed = 0
    repeatsFailed = 0
    CTCFfailed = 0
    homFailed = 0
    restrictionsFailed = 0
    if os.path.exists(options.oligos):
        if options.oligos[-2:] == "gz":
            file = gzip.open(options.oligos, 'rt')  # , encoding='utf-8')
        else:
            file = open(options.oligos)
        for cid, seq in fastaReader(file):
            cseq = LEFT + seq + MIDDLE[:5]
            failed = False
            if cid[-1] == ')':
                region = cid.split("_")[-1].strip("()")
                rchrom = region.split(":")[0]  # .replace("chr", "")
                rstart, rend = map(int, region.split(":")[-1].split("-"))

                # filter simple repeats
                for line in repeatIndex.fetch(rchrom, rstart, rend):
                    fields = line.split("\t")
                    tstart, tend = int(fields[1]), int(fields[2])
                    if (min(tend, rend)-max(tstart, rstart))/float(rend-rstart) > options.repeat:
                        failed = True
                        repeatsFailed += 1

                # filter TSS
                if any(TSSIndex.fetch(rchrom, rstart, rend)):
                    failed = True
                    TSSfailed += 1

                # filter CTCF
                if any(CTCFIndex.fetch(rchrom, rstart, rend)):
                    failed = True
                    CTCFfailed += 1

            if failed:
                continue
            if (options.maxHomLength == None) or (nucleotideruns(cseq) <= options.maxHomLength):
                foundRestrictionSite = False
                for ind, (site, pos) in enumerate(restriction_sites):
                    for match in site.finditer(cseq):
                        foundRestrictionSite = True
                        break
                if not foundRestrictionSite:
                    if seq not in selectedSeqs:
                        selectedSeqs[seq] = cid
                    else:
                        if cid.split("_")[0] not in selectedSeqs[seq]:
                            selectedSeqs[seq] = cid + ":" + selectedSeqs[seq]
                else:
                    restrictionsFailed += 1
            else:
                homFailed += 1
        file.close()

    nrFailed = homFailed + repeatsFailed + TSSfailed + CTCFfailed + restrictionsFailed
    print("""Failed sequences: 
    %d due to homopolymers 
    %d due to simple repeats 
    %d due to TSS site overlap 
    %d due to CTCF site overlap 
    %d due to EcoRI or SbfI restriction site overlap""" % (homFailed, repeatsFailed, TSSfailed, CTCFfailed, restrictionsFailed))
    print("Regular failed: %d" % (nrFailed))
    print("Total sequences included so far: %d" % (len(selectedSeqs)))
