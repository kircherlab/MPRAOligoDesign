
# Create reference folder

mkdir reference

## Reference genome

```bash
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz; bgzip hg38.fa
samtools faidx hg38.fa.gz
cd ..
```

## Simple repeat annotation

```bash
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
zcat simpleRepeat.txt.gz | cut -f 2- | sort -k1,1 -k2,2n -k3,3n | bgzip -c > simpleRepeat.bed.gz
tabix -p bed simpleRepeat.bed.gz
cd ..
```

## TSS positions

```bash
cd reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
zcat gencode.v42.annotation.gtf.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($3 == "transcript") { if ($7 == "+") { print $1,$4-1,$4,".",0,"+" } else { print $1,$5-1,$5,".",0,"-" } } }' | sort -k1,1 -k2,2n -u | bgzip -c > TSS_pos.bed.gz
tabix -p bed simpleRepeat.tsv.gz
cd ..
```

# CTCF sites (ChIP peak + motif)

- Previously identified CTCF sites are in the data folder. Other definitions should be considered.

data/CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz


# Applying available scripts

## Checking VCF files for minimal variant representation and correct reference allele

```bash
for i in *.vcf.gz; do zcat $i | python VCFrefFix.py --refCheck | bgzip -c > ${i/.vcf.gz/.fix.vcf.gz}; done
```

## Generating centered Ref/Alt sequences from VCF file

```bash
for i in *.vcf.fix.gz; do zcat $i | python varsToFrag.py | bgzip -c > ${i/.vcf.fix.gz/snvs.fa.gz}; done
```

## Design from pre-defined sequences with including barcodes in synthesis
- Note their is some deactivated tiling code based on BED files in this source code file
- Misses filter for TSS or CTCF, but probably easily added using similar code as for simple repeats
- Recently, we would no longer design with barcodes, but add barcodes through amplification; here only the outer sequences are added to the design 

```bash
./design.py -s input1.snvs.fa.gz,input1.snvs.fa.gz -r control_regions.fa.gz -o new_design.tsv 
```

