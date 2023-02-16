reference=/fast/projects/cubit/current/static_data/reference/hg38/ucsc/hg38.fa
for i in DNase_NegControls.BloodT-cell.bed.gz DNase_NegControls.Brain.bed.gz DNase_PosControls.bed.gz; do
    name=`basename $i .bed.gz`
    bedtools getfasta -nameOnly -fi $reference -bed <(zcat $i | sed 's/ /_/') | \
    awk '{if ($0 ~ /^>/){print $0} else { print toupper($0)}}' > $name.tmp
    
    fasta-dinucleotide-shuffle -f $name.tmp -s 42 -c 200 | \
    python ../../../workflow/scripts/kMerFilter.py -i $name.tmp -k 6 > $name.diShuffled.fa
    rm $name.tmp
done;