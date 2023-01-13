#!/bin/bash
(
    cat GABAergic_neurons_active.fasta | egrep ">" | sed 's/>//g' | \
        egrep -v "^GABA_shuffle" | \
        sed 's/[-]/\t/' | sed 's/[:_]/\t/g' | awk -v "OFS=\t" '{print $2,$3-1,$4,$1"|"$2":"$3"-"$4"|"$5"|"$6,".",$5}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > GABAergic_neurons_active.bed.gz

