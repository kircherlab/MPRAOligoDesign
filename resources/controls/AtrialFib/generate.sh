#!/bin/bash
(
    cat STARR-seq.atrial.fib.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
        sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-10,$3+10,$9"|"$10"|"$11,".","+"}';
) | sort -k 1,1 -k2,2n | uniq | \
bedtools merge -s -c 4,5,6 -o distinct -d -1 -i - | \
bgzip -c > STARR-seq.atrial.fib.shifted.controls.bed.gz

(cat << EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=target,Number=1,Type=String,Description="Relative gene">
#CHROM	POS	ID REF	ALT	QUAL	FILTER	INFO
EOF
cat STARR-seq.atrial.fib.shifted.controls.fasta | egrep ">" | sed 's/>//g' | \
    sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$4,$9,$6,$7,".","PASS","target="$10}' | \
    sort -k 1,1 -k2,2n | uniq;
) | \
bgzip -c > STARR-seq.atrial.fib.shifted.controls.vcf.gz
