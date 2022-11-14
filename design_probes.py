#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: martin.kircher@bih-charite.de
:Date: *10.08.2021
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import pysam
import re
import random
import gzip

translation = {"B":"[CGT]","D":"[AGT]","H":"[ACT]","K":"[GT]",
               "M":"[AC]","N":".","R":"[AG]","S":"[CG]",
               "V":"[ACG]","W":"[AT]","Y":"[CT]",
               "A":"A","C":"C","G":"G","T":"T"}
table = str.maketrans('TGCABDHKMRSVWYtgcabdhkmrsvwy','ACGTVHDMKYSBWRACGTVHDMKYSBWR') # COMPLEMENT DNA IGNORING CASE

def nucleotideruns(seq):
  longestrun = 0
  lbase,llen = 'N',0
  for base in seq:
    if base == lbase:
      llen+=1
    else:
      if llen > longestrun: longestrun = llen
      llen = 1
      lbase = base
  return longestrun

def Site2RegEx(seq):
  reqEx = ""
  for elem in seq:
    if elem in translation:
      reqEx+=translation[elem]
  return reqEx

def fastaReader(filestream):
  name, value = None,None
  for line in filestream:
    if line.startswith(">"):
      if value != None:
        yield name,value
      value = None
      name = line[1:].strip()
    else:
      value = line.strip()
  if value != None:
    yield name,value
  
parser = OptionParser("%prog [options]")
parser.add_option("-b", "--barcodes", dest="barcodes", help="List of possible barcodes (def 'data/15nt_comb5nt_dist2_ind1_filter-final.txt.gz')",default="data/15nt_comb5nt_dist2_ind1_filter-final.txt.gz")
parser.add_option("-r", "--regular", dest="regular", help="List of comma separated files without line pairing to be included in the design (def '')",default="")
parser.add_option("-s", "--snvs", dest="snvs", help="List of comma separated filenames of paired SNV fasta files to be included in the design (def '')",default="")
#parser.add_option("-r", "--regions", dest="regions", help="List of comma separated BED files to be tiled in the design (def '')",default="")
parser.add_option("-t", "--total", dest="total", help="Total number of oligos to design (def 300000)",default=300000,type="int")
parser.add_option("-z", "--max_homopolymer_length", dest="maxHomLength", help="Maximum homopolymer length (def 10)",default=10,type="int") ## deactivated before
parser.add_option("-n", "--replication", dest="replication", help="Number of features per designed sequence (def 50)",default=50,type="int")
parser.add_option("-f", "--repeat", dest="repeat", help="Maximum fraction explained by a single simple repeat annotation (def 0.25)",default=0.25,type="float") 
parser.add_option("-l", "--length", dest="oligoLength", help="Designed oligo length (def 191)",default=191,type="int") ## deactivated before
parser.add_option("-o", "--outfile", dest="outfile", help="Output file (def 'design.tsv')",default="design.tsv")
(options, args) = parser.parse_args()

repeatIndex = pysam.Tabixfile("reference/simpleRepeat.bed.gz")
#genome = pysam.Fastafile("reference/hg38.fa.gz")

left = "AGGACCGGATCAACT"
middle = "CCTGCAGGGAATTC"
right = "CATTGCGTGAACCGA"
oligoLength = options.oligoLength

names = {}
sites = []
names["CCTGCA^GG"] = "SbfI"
sites.append("CCTGCA^GG")
names["G^AATTC"] = "EcoRI"
sites.append("G^AATTC")

restriction_sites = zip(map(lambda site: re.compile(Site2RegEx(site),re.IGNORECASE),sites),map(lambda site: site.find("^"),sites))

bases = set(['A','C','G','T'])

oligos = {}

if not os.path.exists(options.barcodes):
  sys.stderr.write("Error: Valid barcode file required!\n")
  sys.exit()

  
infile = gzip.open(options.barcodes,"rt") if options.barcodes.endswith(".gz") else open(options.barcodes)
barcodes = []
for line in infile:
  if line.startswith("#"): continue
  barcode = line.split()[0]
  seq = middle[-5:]+barcode+right

  foundRestrictionSite = False
  for ind,(site,pos) in enumerate(restriction_sites):
    for match in site.finditer(seq):
      foundRestrictionSite = True
      break
      break
  
  if not foundRestrictionSite:
    barcodes.append(barcode)
infile.close()
random.shuffle(barcodes)

print("Read barcodes: %d"%len(barcodes))
print("Site replication: %d\n"%(options.replication))

selectedSeqs = {}
increasedReplicationSeqs = set() # for pos control sequences

## READ REGULAR SEQUENCES FROM FASTA
regFailed = 0
for filename in options.regular.split(','):
  if os.path.exists(filename):
    incRep = "posCtl" in filename
    infile = open(filename)
    for cid,seq in fastaReader(infile):
      cseq = left+seq+middle[:5]
      failed = False
      if "[" in cid:
        region = cid.split("_")[-1].strip("[]")
        rchrom = region.split(":")[0].replace("chr","")
        rstart,rend = map(int,region.split(":")[-1].split("-"))
        for line in repeatIndex.fetch(rchrom,rstart,rend):
          fields = line.split("\t")
          tstart,tend = int(fields[1]),int(fields[2])
          if (min(tend,rend)-max(tstart,rstart))/float(rend-rstart) > options.repeat: failed = True
      if failed: 
        regFailed += 1
        continue
      if (options.maxHomLength == None) or (nucleotideruns(cseq) <= options.maxHomLength):
        foundRestrictionSite = False
        for ind,(site,pos) in enumerate(restriction_sites):
          for match in site.finditer(cseq):
            foundRestrictionSite = True
            break
            break
        if not foundRestrictionSite:
          if incRep:
            increasedReplicationSeqs.add(seq)
          if seq not in selectedSeqs:
            selectedSeqs[seq] = cid
          else:
            if cid.split("_")[0] not in selectedSeqs[seq]: 
              selectedSeqs[seq] = cid+":"+selectedSeqs[seq]
        else:
          regFailed += 1
      else:
        regFailed += 1
    infile.close()
print("Regular failed: %d"%(regFailed))
print("Total sequences included so far: %d"%(len(selectedSeqs)))


seqsFailed = 0
## READ PAIRED SEQUENCES (i.e. two entries from fasta file)
fpass,fseq,fname = True, None, None
for filename in options.snvs.split(','):
  #incRep = False 
  incRep = "posCtl" in filename
  if os.path.exists(filename):
    infile = open(filename)
    seqCounter = 0
    seqPassCounter = 0
    for seqind,(cid,seq) in enumerate(fastaReader(infile)):
      seqCounter +=1
      if seqind % 2 == 0:
        fpass,fseq,fname = True, None, None

      failed = False
      if "[" in cid:
        region = cid.split("_")[-1].split("[")[1].strip("]")
        rchrom = "chr"+region.split(":")[0].replace("chr","")
        rstart,rend = map(int,region.split(":")[-1].split("-"))
        for line in tabIndex.fetch(rchrom,rstart,rend):
          fields = line.split("\t")
          tstart,tend = int(fields[1]),int(fields[2])
          if (min(tend,rend)-max(tstart,rstart))/float(rend-rstart) > options.repeat: failed = True
      if failed: 
        if seqind % 2 == 1: seqsFailed += 2
        else: fpass = False
        continue
          
      cseq = left+seq+middle[:5]
      if (options.maxHomLength == None) or (nucleotideruns(cseq) <= options.maxHomLength):
        foundRestrictionSite = False
        for ind,(site,pos) in enumerate(restriction_sites):
          for match in site.finditer(cseq):
            foundRestrictionSite = True
            break
            break
        if not foundRestrictionSite:
          if seqind % 2 == 1:
            if fpass:
              if incRep:
                increasedReplicationSeqs.add(seq)
                increasedReplicationSeqs.add(fseq)
                
              if fseq not in selectedSeqs:
                selectedSeqs[fseq]=fname
              else:
                if fname.split("_")[0] not in selectedSeqs[fseq]: 
                  selectedSeqs[fseq] = selectedSeqs[fseq].split("_")[0]+"-"+fname.split("_")[0]+"_"+"_".join(selectedSeqs[fseq].split("_")[1:])

              if seq not in selectedSeqs:
                selectedSeqs[seq]=cid
              else:
                if cid.split("_")[0] not in selectedSeqs[seq]: 
                  selectedSeqs[seq] = selectedSeqs[seq].split("_")[0]+"-"+cid.split("_")[0]+"_"+"_".join(selectedSeqs[seq].split("_")[1:])              
              seqPassCounter+=1
            else:
              seqsFailed += 2
          else:
            fseq,fname = seq, cid
        else:
          if seqind % 2 == 1: seqsFailed += 2
          else: fpass = False
      else:
        if seqind % 2 == 1: seqsFailed += 2
        else: fpass = False
    infile.close()
    print("File %s: %d sites (%d sequences) [%d sites PASS]"%(filename,seqCounter//2,seqCounter,seqPassCounter))

print("\nIncluded sites: %d (%d sequences) / Failed sites: %d\n"%(len(selectedSeqs)//2,len(selectedSeqs),seqsFailed//2))

remaining = (options.total-len(selectedSeqs)*options.replication-len(increasedReplicationSeqs)*(100-options.replication))
#print remaining,remaining//options.replication,remaining-(remaining//options.replication)*options.replication

#tbases = 0
#for filename in options.regions.split(','):
  #if os.path.exists(filename):
    #infile = open(filename)
    #for line in infile:
      #fields = line.split()
      #if len(fields) > 3:
        #chrom,start,end = fields[0],int(fields[1]),int(fields[2])
        #tbases+=end-start-oligoLength/2
#tiling = tbases//(remaining//options.replication)
#print "Inferred tiling:",tiling

#tileFailed,tilePassed = 0,0
#for filename in options.regions.split(','):
  #if os.path.exists(filename):
    #infile = open(filename)
    #for line in infile:
      #fields = line.split()
      #if len(fields) > 3:
        #chrom,start,end = fields[0],int(fields[1]),int(fields[2])
        #name = fields[-1]
        #seq = genome.fetch(chrom,start,end)
        #for pos in xrange(0,len(seq)-oligoLength+tiling,tiling):
          #subseq = seq[pos:pos+oligoLength]
          #rstart,rend = start+pos+1,start+pos+oligoLength
          #if len(subseq) != oligoLength: 
            #subseq = seq[-oligoLength:]
            #rstart,rend = end-oligoLength+1,end
          #seqname = "Tiles:%s:%s:%d-%d_[chr%s:%d-%d]"%(name,chrom,rstart,rend,chrom,rstart,rend)
          
          #failed = False
          #for line in tabIndex.fetch(rchrom,rstart-1,rend):
            #tfields = line.split("\t")
            #tstart,tend = int(tfields[1]),int(tfields[2])
            #if (min(tend,rend)-max(tstart,rstart))/float(rend-rstart) > options.repeat: failed = True
          #if failed: 
            #tileFailed += 1
            #continue
          
          #cseq = left+subseq+middle[:5]
          #if (options.maxHomLength == None) or (nucleotideruns(cseq) <= options.maxHomLength):
            #foundRestrictionSite = False
            #for ind,(site,pos) in enumerate(restriction_sites):
              #for match in site.finditer(cseq):
                #foundRestrictionSite = True
                #break
                #break
            #if not foundRestrictionSite:
              #tilePassed += 1
              #if subseq not in selectedSeqs:
                #selectedSeqs[subseq] = seqname
              #else: 
                #if seqname.split("_")[0].replace("Tiles:","") not in selectedSeqs[subseq]: 
                  #selectedSeqs[subseq] = selectedSeqs[subseq].split("_")[0]+"-"+seqname.split("_")[0].replace("Tiles:","")+"_"+"_".join(selectedSeqs[subseq].split("_")[1:])              
          #else:
            #tileFailed += 1
#print "Created tiles:",tilePassed,"/ Failed tiles:",tileFailed

replication = (options.total-len(increasedReplicationSeqs)*100)//(len(selectedSeqs)-len(increasedReplicationSeqs))
remaining = options.total-replication*(len(selectedSeqs)-len(increasedReplicationSeqs))
replication2 = 0 if len(increasedReplicationSeqs) == 0 else remaining//len(increasedReplicationSeqs)
print("Inferred replication (%d): %d"%(len(selectedSeqs)-len(increasedReplicationSeqs),replication))
print("Inferred increased replication (%d): %d"%(len(increasedReplicationSeqs),replication2))

for seq,name in selectedSeqs.items():
  if seq not in increasedReplicationSeqs:
    for i in range(replication):
      oligos["%s:%03d"%(name,i+1)]=left+seq+middle+barcodes.pop()+right
  else:
    for i in range(replication2):
      oligos["%s:%03d"%(name,i+1)]=left+seq+middle+barcodes.pop()+right

oremaining=options.total-len(oligos) 
print("Oligos remaining: %d\n"%(oremaining))
rseqs = list(selectedSeqs.keys())
random.shuffle(rseqs)
while options.total-len(oligos) > 0:
  seq = rseqs.pop()
  if seq not in increasedReplicationSeqs:
    oligos["%s:%03d"%(selectedSeqs[seq],replication+1)]=left+seq+middle+barcodes.pop()+right

if oremaining > 0:
  print("Remaining oligos filled with additional replicates.")

print("Final number of oligos: %d\n"%(len(oligos)))

outfile = open(options.outfile,'w')
#outfile.write('#name\tsequence\n') ## Twist does not like the header line
toolong = 0
nhist = defaultdict(int)
for key,value in sorted(oligos.items()):
  nlength = nucleotideruns(value)
  nhist[nlength]+=1
  if nlength >= 8:
    toolong += 1
  outfile.write("%s\t%s\n"%(key,value))
outfile.close()

print("Oligos with long homopolymers (>=8): %d %d %.2f%%"%(toolong,len(oligos),toolong/float(len(oligos))*100))
for key,value in sorted(nhist.items()):
  print("%d\t%d"%(key,value))
