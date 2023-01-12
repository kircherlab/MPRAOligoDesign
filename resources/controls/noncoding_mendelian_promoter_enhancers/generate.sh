#!/bin/bash

(cat << EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=gene,Number=1,Type=String,Description="Relative gene">
##INFO=<ID=PMID,Number=1,Type=Integer,Description="Pubmed ID of publication">
##INFO=<ID=Enhancer,Number=1,Type=Flag,Description="Enhancer variant">
##INFO=<ID=Promoter,Number=1,Type=Flag,Description="Promoter variant">
#CHROM	POS	ID REF	ALT	QUAL	FILTER	INFO
EOF
cat NonCoding.Mendelian.Enhancer.Promoter.Mutation.control.sequences.fasta | egrep ">" | sed 's/>//g' | sed 's/_/\t/g' | awk -v "OFS=\t" '{print $1,$4,$1":"$4$6">"$7"|"$9,$6,$7,".","PASS","gene="$9";PMID="$10";"$11}' | sort -k 1,1 -k2,2n | uniq; ) | \
bgzip -c > noncoding_mendelian_promoters_enhancers.vcf.gz
