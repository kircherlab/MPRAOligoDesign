#!/usr/bin/env python

import re
import pysam
import os
import pandas as pd
from optparse import OptionParser
import gzip
from filter import seqs_filter, regions_filter

def write_output(failed, map_file, remove_unused, out_file):

    map = pd.read_csv(map_file, sep="\t")

    total = len(pd.concat([map["REF_ID"], map["ALT_ID"]]).unique())

    for rid in failed["regions"]:
        if rid in map["Region"].values:
            map = map.drop(map[map["Region"] == rid].index)

    for sid in failed["seqs"]:
        if sid in map["REF_ID"].values:
            map = map.drop(map[map["REF_ID"] == sid].index)
        if remove_unused and sid in map["ALT_ID"].values:
            map = map.drop(map[map["ALT_ID"] == sid].index)

    map.to_csv(out_file, compression='gzip', sep="\t", index=False)

    removed = total - len(pd.concat([map["REF_ID"], map["ALT_ID"]]).unique())

    return (total, removed)


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
                      help="Bedfile for TSS positions", default="reference/TSS_pos.bed.gz")
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
    failedSeqs = []

    fail_reasons = {}
    failed = {}

    if options.regions and os.path.exists(options.regions):
        failed["regions"], fail_reasons = regions_filter(options.regions, repeatIndex, TSSIndex, CTCFIndex, options.repeat)

    if options.seqs and os.path.exists(options.seqs):
        failed["seqs"], reasons = seqs_filter(options.seqs, options.maxHomLength)
        fail_reasons |= reasons

    total, removed = write_output(failed, options.seq_map, remove_unused, options.map_out)

    with open(options.variants_failed_out, 'wt') as out:
        out.writelines(failedSeqs)

    print("""Failed sequences: 
    %d due to homopolymers 
    %d due to simple repeats 
    %d due to TSS site overlap 
    %d due to CTCF site overlap 
    %d due to EcoRI or SbfI restriction site overlap""" % (fail_reasons.get("hompol", 0), fail_reasons.get("repeats", 0), fail_reasons.get("TSS", 0), fail_reasons.get("CTCF", 0), fail_reasons.get("restrictions", 0)))
    print("Total failed: %d" % (removed))
    print("Total passed sequences: %d" % (total-removed))
