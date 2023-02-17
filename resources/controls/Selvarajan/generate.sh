#!/bin/bash
(
    cat Selvarajan.HepG2.STARR-seq.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
        sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-10,$3+10,$9"|"$10,".","+"}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > Selvarajan.HepG2.STARR-seq.shifted.controls.bed.gz

(cat << EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=rsid,Number=1,Type=String,Description="RS id">
#CHROM	POS	ID REF	ALT	QUAL	FILTER	INFO
EOF
cat Selvarajan.HepG2.STARR-seq.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
    sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$4,$9,$6,$7,".","PASS","rsid="$9}' | \
    sort -k 1,1 -k2,2n | uniq;
) | \
bgzip -c > Selvarajan.HepG2.STARR-seq.shifted.controls.vcf.gz
