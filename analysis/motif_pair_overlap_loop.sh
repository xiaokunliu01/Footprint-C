#!/bin/bash

cp /mnt/disk2/1/share/DNase-C/loop/K562_DNaseC_1200M_5kb.bedpe ./
sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n -S 9G K562_DNaseC_1200M_5kb.bedpe | awk '{print $0"\t"NR}' > K562_DNaseC_1200M_5kb_sort.bedpe
awk '{print $1"\t"$2"\t"$3"\t"NR}' K562_DNaseC_1200M_5kb_sort.bedpe > K562_DNaseC_1200M_5kb_1.bed
awk '{print $4"\t"$5"\t"$6"\t"NR}' K562_DNaseC_1200M_5kb_sort.bedpe > K562_DNaseC_1200M_5kb_2.bed

awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$11"\t"$4"\t"$12"\t"$6"\t"$14"\t"$17}' /home/lxk/private/DNase-C/DNaseC_paper/heatmap_CTCF_CTCF/CTCF/CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n -S 9G | awk '{print $0"\t"NR}' > CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt
/home/lxk/share/software/pgltools-3.0.1/sh/pgltools intersect -a CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt -b ../K562_DNaseC_1200M_5kb_sort.bedpe -wo > motif_pair_intersect_loop.txt


for i in KLF1 MAZ SP3
do
mkdir $i
cd $i
awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$11"\t"$4"\t"$12"\t"$6"\t"$14"\t"$17}' /home/lxk/private/DNase-C/DNaseC_paper/heatmap_CTCF_CTCF/$i/$i'_'$i'_'motif_pair_cis_dis320bp_groupby.txt | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n -S 9G | awk '{print $0"\t"NR}' > $i'_'$i'_'motif_pair_cis_dis320bp_groupby.txt
/home/lxk/share/software/pgltools-3.0.1/sh/pgltools intersect -a $i'_'$i'_'motif_pair_cis_dis320bp_groupby.txt -b ../K562_DNaseC_1200M_5kb_sort.bedpe -wo > motif_pair_intersect_loop.txt
for (( j=0; j<12; j=j+1 ))
do
  awk -v a=$j '{if($11>a)n=n+1}END{print a"\t"n}' $i'_'$i'_'motif_pair_cis_dis320bp_groupby.txt >> pair.tab
  awk -v a=$j '{if($19>a)n=n+1}END{print a"\t"n}' motif_pair_intersect_loop.txt >> intersect.tab
done
cd ..
done


