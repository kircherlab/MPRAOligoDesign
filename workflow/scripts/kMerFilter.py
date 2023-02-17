#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *01.09.2016
"""

import sys
import os
from optparse import OptionParser

# cat ~/regulatory_tests/liverEnhancer/design/inserts.fa | ./shuffleFasta -n 100 | ./kMerFilter.py -i ~/regulatory_tests/liverEnhancer/design/inserts.fa -k 6


def read_fastq(filehandle):
    """ Reads fastq and fasta entries from filehandle
        supports multiline fasta and fastq
    """
    count = 0
    seq = ""
    id = ""
    qual = None
    is_fasta = False
    for line in filehandle:
        line = line.rstrip("\n\r")
        if ((count == 0) and (line.startswith("@") or line.startswith(">"))):  # Read identifier
            id = line[1:]
            count += 1
            is_fasta = line.startswith(">")

        elif count == 1:        # read sequence
            seq = line
            count += 1

        # multiple case: a) quality identifier (fastq), b) more sequence (fastq,fasta), c) next sequence identifier (fasta)
        elif count == 2:
            if line.startswith("+"):  # case a)
                id2 = line[1:]
                if ((len(id2) > 0) and (id2 != id)):  # optional sanity check
                    sys.stderr.write("[NOTE] sequences identifier does not match quality identifier: " + id + " " + id2 + "\n")
                count += 1

            elif is_fasta and (line.startswith(">") or line.startswith("@")):  # case c)
                yield id, seq, None
                id, seq, qual = None, None, None
                count = 1
                id = line[1:]
                is_fasta = line.startswith(">")

            else:  # case b)
                seq = seq + line

        elif count == 3:
            if qual == None:
                qual = line
            else:
                qual = qual + line

            if (len(qual) > len(seq)):  # Another sanity check
                sys.stderr.write("[NOTE] sequence and quality line differ in length\n")
            if (len(qual) >= len(seq)):
                count = 0
                yield id, seq, qual
                id, seq, qual = None, None, None
        else:
            sys.stderr.write("Unexpected line:" + str(line.strip()) + "\n")
            count = 0

    if id and seq:
        yield id, seq, qual



parser = OptionParser("%prog [options]")
parser.add_option("-i", "--infile", dest="infile", help="Input fasta filename (def input.fa)", default="input.fa")
parser.add_option("-k", "--kmers", dest="kmers", help="Kmer length to filter (def 6)", default=6, type="int")
parser.add_option("--inclMinOverlap", dest="inclMinOverlap",
                  help="Include minimum kmer overlap in output file (def Off)", default=False, action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn on verbose output (def Off)",
                  default=False, action="store_true")
(options, args) = parser.parse_args()

if not os.path.exists(options.infile):
    sys.stderr.write("Input Fasta file does not exist.\n")
    sys.exit()

seqDict = {}
infile = open(options.infile)
for cid, seq, _ in read_fastq(infile):
    seqDict[cid] = seq
infile.close()

outputDict = {}
lKmers = set()
lName = None
for cid, seq, _ in read_fastq(sys.stdin):
    cname = "_".join(cid.split("_")[:-1])
    if cname not in seqDict:
        sys.stderr.write("Sequence missing from input fasta.\n")
        sys.exit()
    else:
        cKmers = set()
        for i in range(0, len(seq)-options.kmers+1):
            cKmers.add(seq[i:i+options.kmers])
        if cname != lName:
            lName = cname
            lKmers = set()
            rseq = seqDict[cname]
            for i in range(0, len(rseq)-options.kmers+1):
                lKmers.add(rseq[i:i+options.kmers])
        res = len(cKmers.intersection(lKmers))
        if cname not in outputDict:
            outputDict[cname] = seq, res
        elif outputDict[cname][1] > res:
            outputDict[cname] = seq, res

for cid, (seq, kmerCount) in outputDict.items():
    if options.inclMinOverlap:
        sys.stdout.write(">%s %d\n%s\n" % (cid, kmerCount, seq))
    else:
        sys.stdout.write(">%s\n%s\n" % (cid, seq))
