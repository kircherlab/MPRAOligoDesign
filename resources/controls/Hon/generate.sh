#!/bin/bash
(
    cat Hon.congenital.heart.disease.enhancers.fasta | egrep ">" | sed 's/>//g' | \
        sed 's/_/\t/' |sed 's/_/\t/' |sed 's/_/\t/' | awk -v "OFS=\t" '{print $1,$2-1-10,$3+10,$4"|"$1":"$2"-"$3,".","+"}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > Hon.congenital.heart.disease.enhancers.bed.gz

