#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: martin.kircher@bih-charite.de
:Date: *03.08.2021
"""

import sys, os
import pysam
from optparse import OptionParser

parser = OptionParser("%prog [options]")
parser.add_option("-s","--SNVs", dest="SNVs", help="Path to SNV output file (def STDOUT)",default=None)
parser.add_option("-i","--InDels", dest="InDels", help="Path to InDel output file (def STDOUT)",default=None)
parser.add_option("-r","--reference", dest="reference", help="Reference genome to perform reference check and missing context (fasta index needs to be available) (def reference/hg38.fa.gz)",default='reference/hg38.fa.gz')
parser.add_option("--refCheck", dest="refCheck", help="Perform reference check (def off)",default=False,action="store_true")
parser.add_option("-k","--keep", dest="keep", help="Keep additional fields (def off)",default=False,action="store_true")
parser.add_option("-v","--verbose", dest="verbose", help="Report fixed and failed variants (def off)",default=False,action="store_true")
(options, args) = parser.parse_args()

reference = None
if not os.path.exists(options.reference) or not os.path.exists(options.reference+".fai"): 
  options.refCheck = False
  sys.stderr.write("Warning: Reference not available, can not perform reference check or allele fix.\n")
else:
  reference = pysam.Fastafile(options.reference)

def is_nucleotide(seq):
  for base in seq.upper():
    if base not in 'ACGT': return False
  return True

def sharedPrefix(s1,s2):
  minLength = min(len(s1),len(s2))
  shared = 0
  for ind in range(minLength-1):
    if s1[ind] == s2[ind]: shared+=1
    else: break
  if minLength == 1:
    return max(0,shared-1)
  else:
    return shared

def sharedSuffix(s1,s2):
  minLength = min(len(s1),len(s2))-1
  shared = 0
  for ind in range(minLength*-1,0)[::-1]:
    if s1[ind] == s2[ind]: shared+=1
    else: break
  return shared

###################
# VCF FIELDS
###################

fchr = 0
fpos = 1
fdbsnp = 2
fref_allele = 3
falt_allele = 4
fgeno_qual = 5
fflag = 6
finfo = 7
fformat = 8
fvalues = 9

outsnvs = sys.stdout
outindels = sys.stdout
if options.SNVs != None:
  outsnvs = open(options.SNVs,'w')
if options.InDels != None:
  outindels = open(options.InDels,'w')

countSwitch = 0
countHomozygote = 0
countFailed = 0
countInDelFormat = 0

for line in sys.stdin:
  if line.upper().startswith("#CHROM"):
    fields = line.split()
    if options.keep:
      outsnvs.write("\t".join(fields)+"\n")
    else:
      outsnvs.write("\t".join(fields[:falt_allele+1])+"\n")
    if options.SNVs != None or options.InDels != None:
      if options.keep:
        outindels.write("\t".join(fields)+"\n")
      else:
        outindels.write("\t".join(fields[:falt_allele+1])+"\n")
  if line.startswith("#"): continue
  fields = line.split()
  
  if len(fields) > falt_allele:
    altAlleles = fields[falt_allele].upper().split(',')
    if is_nucleotide(fields[fref_allele]): 
      ref = fields[fref_allele].upper()
    elif fields[fref_allele] == "-":
      if reference == None: continue
      pos = int(fields[fpos])
      try:
        ref = reference.fetch(fields[fchr],pos-2,pos-1).upper()
      except:
        try:
          ref = reference.fetch("chr"+fields[fchr],pos-2,pos-1).upper()
        except:
          ref = None
      if ref == None: continue
      fields[fpos] = "%d"%(pos-1)
      altAlleles = map(lambda x:ref+x,altAlleles)
      for alt in altAlleles:
        if options.verbose:
          sys.stderr.write("FIXED: %s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
        countInDelFormat+=1
    else:
      continue

    nprefix = False
    naltAlleles = []
    for alt in altAlleles:
      if is_nucleotide(alt): 
        naltAlleles.append(alt.upper())
      elif alt == "-":
        nprefix = True
        naltAlleles.append("")
    altAlleles = naltAlleles
    
    if nprefix:
      if reference == None: continue
      pos = int(fields[fpos])
      try:
        prefixSeq = reference.fetch(fields[fchr],pos-2,pos-1).upper()
      except:
        try:
          prefixSeq = reference.fetch("chr"+fields[fchr],pos-2,pos-1).upper()
        except:
          prefixSeq = None
      if prefixSeq == None: continue
      fields[fpos] = "%d"%(pos-1)
      ref = prefixSeq+ref
      altAlleles = map(lambda x:prefixSeq+x,altAlleles)
      for alt in altAlleles:
        if options.verbose:
          sys.stderr.write("FIXED: %s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
        countInDelFormat+=1
        
    if options.refCheck:
      pos = int(fields[fpos])
      try:
        href = reference.fetch(fields[fchr],pos-1,pos-1+len(ref)).upper()
      except:
        try:
          href = reference.fetch("chr"+fields[fchr],pos-1,pos-1+len(ref)).upper()
        except:
          href = None
      if (href != ref) and (href not in altAlleles) and (ref != alt): 
        for alt in altAlleles:
          if options.verbose:
            sys.stderr.write("FAILED: %s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
          countFailed+=1
        continue
      elif (href != ref) and (href in altAlleles):
        altAlleles.remove(href)
        altAlleles.append(ref)
        ref = href
        for alt in altAlleles:
          if options.verbose:
            sys.stderr.write("FIXED: %s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
          countSwitch+=1
      elif (href != ref) and (ref == alt):
        ref = href
        for alt in altAlleles:
          if options.verbose:
            sys.stderr.write("FIXED: %s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
          countHomozygote+=1
          
    for alt in altAlleles:
      if ref == alt: continue

      if len(alt) == len(ref) and len(ref) == 1:
        if options.keep:
          outsnvs.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt,"\t".join(fields[falt_allele+1:])))
        else:
          outsnvs.write("%s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],ref,alt))
      else:
        trimValue = sharedPrefix(ref,alt)
        if trimValue != 0:
          nref = ref[trimValue:]
          nalt = alt[trimValue:]
        else:
          nref,nalt = ref,alt
        trimValue2 = sharedSuffix(nref,nalt)
        if trimValue2 != 0:
          nref = nref[:-trimValue2]
          nalt = nalt[:-trimValue2]
        if len(nalt) == len(nref) and len(ref) == 1:
          if options.keep:
            outsnvs.write("%s\t%d\t%s\t%s\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,fields[fdbsnp],nref,nalt,"\t".join(fields[falt_allele+1:])))
          else:
            outsnvs.write("%s\t%d\t%s\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,fields[fdbsnp],nref,nalt))
        elif (trimValue == 0):
          if options.keep:
            outindels.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],nref,nalt,"\t".join(fields[falt_allele+1:])))
          else:
            outindels.write("%s\t%s\t%s\t%s\t%s\n"%(fields[fchr],fields[fpos],fields[fdbsnp],nref,nalt))
        else:
          if options.keep:
            outindels.write("%s\t%d\t%s\t%s\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,fields[fdbsnp],nref,nalt,"\t".join(fields[falt_allele+1:])))
          else:
            outindels.write("%s\t%d\t%s\t%s\t%s\n"%(fields[fchr],int(fields[fpos])+trimValue,fields[fdbsnp],nref,nalt))

if options.SNVs != None:
  outsnvs.close()
if options.InDels != None:
  outindels.close()

sys.stderr.write("Fixed InDel %d, Switched %d, Homozygote %d, Failed %d\n"%(countInDelFormat,countSwitch,countHomozygote,countFailed))
