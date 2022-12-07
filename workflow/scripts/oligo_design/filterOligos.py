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

def write_output(regions, seqmap, regions_out, design_out, variants_out):

    selected_ids = []
    selected_vars = []
    with open(design_out, 'w') as out_fa:
        for base_id in seqmap:
            if seqmap[base_id]["ref"]:
                for id, seq in seqmap[base_id]["seqs"].items():
                    out_fa.write(">{}\n{}\n".format(id, seq))
                    if len(id.split("_ID")) == 2:
                        selected_vars.append(id.split("_")[-1])
                selected_ids.append(base_id)
    with gzip.open(regions_out, 'wt') as out_bed:
        filtered = regions.loc[regions.index.isin(selected_ids)]
        filtered = filtered.reset_index(level=0)
        filtered.to_csv(out_bed, compression='gzip', sep="\t", 
        header=False, index=False, columns=["chrom", "start", "end", "id", "qfilter", "strand"])
    with open(variants_out, 'w') as out_vars:
        out_vars.writelines("\n".join(selected_vars))
            


if __name__ == "__main__":
    parser = OptionParser("%prog [options]")
    parser.add_option("-s", "--seqs", dest="seqs",
                      help="File containing the designed oligos (def '')", default="results/oligo_design/test_1/design.fa")
    parser.add_option("-r", "--regions", dest="regions",
                      help="File containing the regions (def '')", default="results/oligo_design/test_1/design.regions.bed.gz")
    parser.add_option("-z", "--max_homopolymer_length", dest="maxHomLength",
                      help="Maximum homopolymer length (def 10)", default=10, type="int")  # deactivated before
    parser.add_option("-f", "--repeat", dest="repeat",
                      help="Maximum fraction explained by a single simple repeat annotation (def 0.25)", default=0.25, type="float")
    parser.add_option("-o", "--output-regions", dest="regions_out",
                      help="Output file of filtered regions")
    parser.add_option("-d", "--output-design", dest="design_out",
                      help="Output file of diltered sequences")
    parser.add_option("-v", "--output-variants", dest="variants_out",
                      help="Output file of filtered variant ids")  
    parser.add_option("-i", "--output-failed-variants", dest="variants_failed_out",
                      help="Output file of failed variant ids")  
    (options, args) = parser.parse_args()

    repeatIndex = pysam.Tabixfile("reference/simpleRepeat.bed.gz") # FIXME: hardcoded path
    TSSIndex = pysam.Tabixfile("reference/TSS_pos.bed.gz") # FIXME: hardcoded path
    CTCFIndex = pysam.Tabixfile(
        "reference/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz") # FIXME: hardcoded path

    names = {}
    sites = []
    names["CCTGCA^GG"] = "SbfI"
    sites.append("CCTGCA^GG")
    names["G^AATTC"] = "EcoRI"
    sites.append("G^AATTC")

    restriction_sites = zip(map(lambda site: re.compile(Site2RegEx(
        site), re.IGNORECASE), sites), map(lambda site: site.find("^"), sites))
    selectedSeqs = {}
    seqFilterMap = {}
    failedSeqs = []

    nrFailed = 0
    TSSfailed = 0
    repeatsFailed = 0
    CTCFfailed = 0
    homFailed = 0
    restrictionsFailed = 0
    if os.path.exists(options.seqs) and os.path.exists(options.regions):
        if options.seqs[-3:] == ".gz":
            seqfile = gzip.open(options.seqs, 'rt') 
        else:
            seqfile = open(options.seqs)
        regions = pd.read_csv(options.regions, names=["chrom", "start", "end", "id", "qfilter", "strand"], sep="\t")
        regions = regions.set_index("id")
        for full_id, seq in fastaReader(seqfile):
            is_ref = full_id[:3] == "REF" 
            cid = full_id[4:].split("_ID")[0]
            cseq = LEFT + seq + RIGHT[:5]
            failed = False
            if cid in regions.index:
                rchrom, rstart, rend,_,_ = regions.loc[cid]
                # filter simple repeats
                for line in repeatIndex.fetch(rchrom, rstart, rend):
                    fields = line.split("\t")
                    tstart, tend = int(fields[1]), int(fields[2])
                    if (min(tend, rend)-max(tstart, rstart))/float(rend-rstart) > options.repeat:
                        failed = True
                        repeatsFailed += 1
                        break

                # filter TSS
                if any(TSSIndex.fetch(rchrom, rstart, rend)):
                    failed = True
                    TSSfailed += 1

                # filter CTCF
                if any(CTCFIndex.fetch(rchrom, rstart, rend)):
                    failed = True
                    CTCFfailed += 1
            
                if not failed:
                    if (options.maxHomLength == None) or (nucleotideruns(cseq) <= options.maxHomLength):
                        foundRestrictionSite = False
                        for ind, (site, pos) in enumerate(restriction_sites):
                            if any(site.finditer(cseq)):
                                foundRestrictionSite = True
                                break
                        if not foundRestrictionSite:
                            if cid not in selectedSeqs:
                                selectedSeqs[cid] = {"ref": is_ref, "seqs": {full_id: seq}}
                            else:
                                selectedSeqs[cid]["ref"] = is_ref or selectedSeqs[cid]["ref"]
                                selectedSeqs[cid]["seqs"] |= {full_id: seq}
                        else:
                            failed =True
                            restrictionsFailed += 1
                    else:
                        homFailed += 1
                        failed = True
                
                if failed:
                    nrFailed += 1
                    failedSeqs.append(full_id)
            else:
                print(Warning("ID "+ cid + " does not have associated coordinates in " + options.regions))
        seqfile.close()

    write_output(regions, selectedSeqs, options.regions_out, options.design_out, options.variants_out)

    with open(options.variants_failed_out, 'wt') as out:
        out.writelines(failedSeqs)

    print("""Failed sequences: 
    %d due to homopolymers 
    %d due to simple repeats 
    %d due to TSS site overlap 
    %d due to CTCF site overlap 
    %d due to EcoRI or SbfI restriction site overlap""" % (homFailed, repeatsFailed, TSSfailed, CTCFfailed, restrictionsFailed))
    print("Total failed: %d" % (nrFailed))
    print("Total passed sequences: %d" %(len(selectedSeqs)))
