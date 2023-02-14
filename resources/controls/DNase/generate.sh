#!/bin/bash
(
    cat DNase_NegControls.BloodT-cell.fa | egrep ">" | sed 's/>//g' | \
        sed 's/[-]/\t/' | sed 's/[:_]/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-50,$3+49,$1":"$2-50"-"$3+49"_"$4"_"$5"_"$6"_"$7,".","."}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > DNase_NegControls.BloodT-cell.bed.gz

(
    cat DNase_NegControls.Brain.fa | egrep ">" | sed 's/>//g' | \
        sed 's/[-]/\t/' | sed 's/[:_]/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-50,$3+49,$1":"$2-50"-"$3+49"_"$4"_"$5"_"$6"_"$7,".","."}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > DNase_NegControls.Brain.bed.gz


(
    cat DNase_PosControls.fa | egrep ">" | sed 's/>//g' | \
        sed 's/[-]/\t/' | sed 's/[:_]/\t/g' | awk -v "OFS=\t" '{print $1,$2-1-50,$3+49,$1":"$2-50"-"$3+49"_"$4"_"$5"_"$6,".","."}';
) | sort -k 1,1 -k2,2n | uniq | \
bgzip -c > DNase_PosControls.bed.gz

