#!/bin/bash
cat original_regions.bed | sort -k 1,1 -k2,2n | \
bedtools slop -s -l 6 -r 0 -g /fast/projects/cubit/current/static_data/reference/GRCh38/hs38/hs38.fa.genome -i - | 
sort -k 1,1 -k2,2n | uniq | \
bgzip -c > Engreitz.housekeeping.promoters.controls.bed.gz

