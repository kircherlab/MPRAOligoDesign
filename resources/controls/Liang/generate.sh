#!/bin/bash
(
    cat Liang.neuronal.caQTL.eQTL.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
        sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-10,$3+10,$9"|"$10,".","+"}';
) | sort -k 1,1 -k2,2n | uniq | \
bedtools merge -s -c 4,5,6 -o distinct -d -1 -i - | \
bgzip -c > Liang.neuronal.caQTL.eQTL.shifted.controls.bed.gz

(cat << EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=rsid,Number=1,Type=String,Description="RS id">
#CHROM	POS	ID REF	ALT	QUAL	FILTER	INFO
EOF
cat Liang.neuronal.caQTL.eQTL.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
    sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$4,$9,$6,$7,".","PASS","rsid="$9}' | \
    sort -k 1,1 -k2,2n | uniq;
) | \
bgzip -c > Liang.neuronal.caQTL.eQTL.shifted.controls.vcf.gz
