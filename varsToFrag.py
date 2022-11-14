#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: martin.kircher@bih-charite.de
:Date: *03.08.2021
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import pysam

genome = pysam.Fastafile("reference/hg38.fa.gz")

parser = OptionParser("%prog [options]")
(options, args) = parser.parse_args()

rsize = 191
prefix = ""
for line in sys.stdin:
  if line.startswith('#'): continue
  fields = line.split()
  if len(fields) >= 5:
    chrom = fields[0].replace('chr','')
    pos = int(fields[1])
    rsid = fields[2]
    ref = fields[3]
    alt = fields[4]
  else:
    continue

  try:
    seq = genome.fetch(prefix+chrom,pos-1,pos+len(ref)-1).upper()
  except:
    try:
      seq = genome.fetch("chr"+chrom,pos-1,pos+len(ref)-1).upper()
      prefix = "chr"
    except:
      sys.stderr.write("Error retrieving sequence from reference: %s:%d-%d\n"%(prefix+chrom,pos,pos+len(ref)-1))
      sys.exit()
  
  if seq != ref:
    sys.stderr.write("Reference base does not match reported allele (%s:%d Ref:%s Gen:%s)!\n"%(ref,seq))
    sys.exit()

  left=rsize//2-len(ref)//2
  
  seq = genome.fetch(prefix+chrom,pos-left-1,pos+(rsize-left-1)).upper()
  sys.stdout.write(">%s:%d:%s>%s:Ref_%s[chr%s:%d-%d]\n%s\n"%(chrom,pos,ref,alt,rsid,chrom,pos-left,pos+(rsize-left-1),seq))

  left=rsize//2-len(alt)//2
  seq = genome.fetch(prefix+chrom,pos-left-1,pos+(rsize-left+len(ref)-len(alt)-1)).upper()
  seq = seq[:left]+alt+seq[left+len(ref):rsize+len(ref)-len(alt)]
  sys.stdout.write(">%s:%d:%s>%s:Alt_%s[chr%s:%d-%d]\n%s\n"%(chrom,pos,ref,alt,rsid,chrom,pos-left,pos+(rsize-left+len(ref)-len(alt)-1),seq))
  #chrom,pos,ref,alt,"%s[%s/%s]%s"%(seq[:rsize],seq[rsize:rsize+len(ref)],alt,seq[rsize+len(ref):])
