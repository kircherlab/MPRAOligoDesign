#!/bin/bash
(
    cat Kircher.controls.fasta | egrep ">" | sed 's/>//g' | \
        sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-10,$3+10,$9"|"$10,".","+"}' | \
        sed 's/:/|/g' ;
) | sort -k 1,1 -k2,2n | uniq | \
bedtools merge -s -c 4,5,6 -o distinct -d -1 -i - | \
bgzip -c > Kircher.controls.bed.gz

(cat << EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM	POS	ID REF	ALT	QUAL	FILTER	INFO
EOF
cat Kircher.controls.fasta | egrep ">" | sed 's/>//g' | \
    sed 's/_/\t/g' | sed 's/[.:]/_/g' | awk -v "OFS=\t" '{print $1,$4,$9,$6,$7,".","PASS","."}' | \
    sort -k 1,1 -k2,2n | uniq;
) | \
bgzip -c > Kircher.controls.vcf.gz
