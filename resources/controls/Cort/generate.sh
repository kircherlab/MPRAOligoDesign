#!/bin/bash
(
    cat Primary_fetal_cortical_cells-active-hg38.fasta | egrep ">" | sed 's/>//g' | \
        egrep -v "^negative_da" | egrep -v "^negative_vaccarino" | egrep -v "^positive_white_lab" | \
        sed 's/[-:]/\t/g' | sed 's/_/\t/' | sed 's/_/\t/' | awk -v "OFS=\t" '{print $2,$3,$4,$1"|"$2":"$3+1"-"$4"|"$5"|"$6,".","+"}';
    cat Primary_fetal_cortical_cells-active-hg38.fasta | egrep ">" | sed 's/>//g' | \
        grep "positive_white_lab" | \
        sed 's/[-:]/\t/g' | sed 's/_/-/' | sed 's/_/-/' | sed 's/_/\t/g' | awk -v "OFS=\t" '{print $2,$3,$4,$1"|"$2":"$3+1"-"$4"|"$5,".","+"}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > primary_fetal_cortical_cells-active-hg38.bed.gz

