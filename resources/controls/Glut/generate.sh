#!/bin/bash
(
    cat Glutamatergic_neurons_active.fasta | egrep ">" | sed 's/>//g' | \
        egrep -v "^Glut_shuffle35" | \
        sed 's/[-]/\t/' | sed 's/[:_]/\t/g' | awk -v "OFS=\t" '{print $2,$3-1,$4,$1"|"$2":"$3"-"$4"|"$5"|"$6,".",$5}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > Glutamatergic_neurons_active.bed.gz

