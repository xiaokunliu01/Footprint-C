#!/bin/bash

# 1. takes file: .pair file and hg38_chromsize.txt file
#    applies tool: genomeCoverageBed and bedGraphToBigWig
#    produces output: .bw file

# 2. takes file: .bw file and .motif file
#    applies tool: computeMatrix and plotProfile in deeptools
#    produces output: average profile plot


# enrichment
mkdir enrichment/
cd enrichment/

awk '{print($1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6)}' ../fragment_contact_pairs/$1'.'pair > $1'.'bed
a=$(cat $1'.'bed | wc -l)
sort -k1,1 -k2,2n -k3,3n -S 9G $1'.'bed | genomeCoverageBed -bg -i - -g /ssd/genome/hg38_chromsize.txt | awk -v a=$a '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/a*1000000}' | awk '{$4/=1;print}' OFS='\t' > $1'.'bdg
bedGraphToBigWig $1'.'bdg /ssd/genome/hg38_chromsize.txt $1'.'bw

rm $1'.'bed
rm $1'.'bdg

computeMatrix reference-point --referencePoint center -S $1'.'bw -R motif/K562_CTCF_RAD21_DNase.motif -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o $1'_'CTCFmotif.gz
plotProfile -m $1'_'CTCFmotif.gz -o $1'_'CTCFmotif_profile.png --outFileNameData $1'_'CTCFmotif_profile.tab --plotHeight 6 --plotWidth 8 --perGroup --refPointLabel 'CTCF motif' --yMin 0 -z '' --samplesLabel $1

rm $1'_'CTCFmotif.gz
