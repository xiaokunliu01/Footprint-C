#!/bin/bash

# 1. takes file: fragment contact pairs file and motif file
#    applies tools: closestBed
#    produces output: v-plot data file



# plot V-plot
mkdir vplot/
cd vplot/

awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)+1"\t"NR"\t"$3-$2"\t"$10"\n"$4"\t"int(($5+$6)/2)"\t"int(($5+$6)/2)+1"\t"NR"_2""\t"$6-$5"\t"$11}' ../fragment_contact_pairs/$1'.'pair | sort -k1,1 -k2,2n -S 9G > $1'.'bed

closestBed -a $1'.'bed -b motif/K562_CTCF_RAD21_DNase.motif -d -t first > $1'_'closest_CTCFmotif.txt

awk -f /home/lxk/private/optionData/vplot.awk $1'_'closest_CTCFmotif.txt > FragmentLengthVsDistance.csv

rm $1'.'bed
rm $1'_'closest_CTCFmotif.txt
