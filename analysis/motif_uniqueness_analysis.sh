#!/bin/bash

# calculate frequency of unique motif annotation of fragments from different Hi-C datasets
# 1. takes file: .pair file
#    applies tool: cooltools
#    produces output: PC1 value file


awk '{print $1"\t"int($2/2+$3/2)"\t"int($2/2+$3/2)+1"\t"NR"\t"$5"\t"$6"\t"$4}' /mnt/disk4/public/RefBed/motif/homer.KnownMotifs.hg38.191020.bed > homer_hg38.motif
cut -f1-6 homer_hg38.motif > homer_hg38_1.motif
bigWigAverageOverBed /mnt/disk4/public/K562/DNase-seq/K562_DNase.bw homer_hg38.motif a.tab -sampleAroundCenter=30
sort -k1,1n a.tab | paste homer_hg38.motif - | awk '{printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\t%s\n",$1,$2,$3,$4,$11*100000,$6}' > K562_hg38_motif.bed
sort -k1,1n a.tab | paste homer_hg38.motif - | awk '{printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\t%s\t%s\t%.0f\n",$1,$2,$3,$4,$12*100000,$6,$7}' > K562_hg38_motif.bed


# FootprintC
closestBed -a FootprintC_K562_FA_UMI_sort_name.singlepair -b K562_hg38_motif_sort.bed -d -t all > FootprintC_K562_closest_motif.txt
sort -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6n -k7,7n -u --buffer-size=10% -T ./ --parallel=80 FootprintC_K562_closest_motif.txt > FootprintC_K562_closest_motif_u.txt
awk '{if($13==0) print $0}' FootprintC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 6 -o min,max > FootprintC_K562_closest_motif_u_0_groupby.txt
awk '{if($13==0 && $6-$2>10 && $3-$7>10) print $0}' FootprintC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 6 -o min,max
awk '{if($13==0 && $6-$2>10 && $3-$7>10) print $0}' FootprintC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 6 -o min,max | awk '{if($6-$5<=20) print $0}' | wc -l


# MicroC
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/ && $8!="no_adapterno_adapterno_adapter") print $0}' /home/lxk/private/Micro-C/Public/K562/fragment_Length/K562_MicroC_R1r3_rmdup_sample.pair > K562_MicroC_R1r3.pair
shuf -n 822793 K562_MicroC_R1r3.pair | awk '{if(int($2/2+$3/2)-75>0 && int($5/2+$6/2)-75>0) print $1"\t"int($2/2+$3/2)-75"\t"int($2/2+$3/2)+75"\n"$4"\t"int($5/2+$6/2)-75"\t"int($5/2+$6/2)+75}' | sort -k1,1 -k2,3n -S 9G | awk '{print $0"\t"NR}' > K562_MicroC_R1r3.singlepair
closestBed -a K562_MicroC_R1r3.singlepair -b ../K562_hg38_motif_sort.bed -d -t all > MicroC_K562_closest_motif.txt
sort -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6n -k7,7n -u --buffer-size=10% -T ./ --parallel=80 MicroC_K562_closest_motif.txt > MicroC_K562_closest_motif_u.txt
awk '{if($13==0) print $0}' MicroC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 6 -o min,max > MicroC_K562_closest_motif_u_0_groupby.txt


# in situ Hi-C
awk '{if($2!~/chr[CLMT]/ && $5!~/chr[CLMT]/) print $0}' ../../decay/insituHiC_K562_R2_val.allValidPairs > insituHiC_K562_R2.allValidPairs
shuf -n 822793 insituHiC_K562_R2.allValidPairs > insituHiC_K562_R2_shuf.allValidPairs
awk '{print $9"\n"$10}' insituHiC_K562_R2_shuf.allValidPairs | sort -k1,1n -S 9G | awk '{print $0"\t"NR}' | join -1 1 -2 4 - hg38_Mbol.bed | awk '{print $3"\t"$4"\t"$5"\t"$2"\t"$1}' | sort -k1,1 -k2,3n -S 9G > insituHiC_K562.singlepair
closestBed -a insituHiC_K562.singlepair -b ../K562_hg38_motif_sort.bed -d -t all > insituHiC_K562_closest_motif.txt
sort -k1,1 -k2,2n -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n -u --buffer-size=10% -T ./ --parallel=80 insituHiC_K562_closest_motif.txt > insituHiC_K562_closest_motif_u.txt
awk '{if($14==0) print $0}' insituHiC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 7 -o min,max > insituHiC_K562_closest_motif_u_0_groupby.txt


# Hi-TrAC
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' /home/lxk/private/K562/FootprintC/public/HiTrAC/Hi-TrAC_K562_R1r1.bedpe > Hi-TrAC_K562_R1r1.pair
shuf -n 822793 Hi-TrAC_K562_R1r1.pair > Hi-TrAC_K562_R1r1_shuf.pair
awk '{if($9=="+" && $10=="+"){print $1"\t"$2"\t"$2+1"\n"$4"\t"$5"\t"$5+1}if($9=="+" && $10=="-"){print $1"\t"$2"\t"$2+1"\n"$4"\t"$6-1"\t"$6}if($9=="-" && $10=="+"){print $1"\t"$3-1"\t"$3"\n"$4"\t"$5"\t"$5+1}if($9=="-" && $10=="-"){print $1"\t"$3-1"\t"$3"\n"$4"\t"$6-1"\t"$6}}' Hi-TrAC_K562_R1r1_shuf.pair | sort -k1,1 -k2,3n -S 9G > Hi-TrAC_K562_R1r1_shuf.singlepair
closestBed -a Hi-TrAC_K562_R1r1_shuf.singlepair -b hg38_MluCI_NlaIII_sort.bed -d -t all | awk '{if($10==0)print $4"\t"$5"\t"$6"\t"NR"\t"$7}' | sort -k1,1 -k2,3n -S 9G > Hi-TrAC_K562.singlepair
closestBed -a Hi-TrAC_K562.singlepair -b ../K562_hg38_motif_sort.bed -d -t all > Hi-TrAC_K562_closest_motif.txt
sort -k1,1 -k2,2n -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n -u --buffer-size=10% -T ./ --parallel=80 Hi-TrAC_K562_closest_motif.txt > Hi-TrAC_K562_closest_motif_u.txt
awk '{if($14==0) print $0}' Hi-TrAC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 7 -o min,max > Hi-TrAC_K562_closest_motif_u_0_groupby.txt


#BL-Hi-C
awk '{if($2!~/chr[CLMT]/ && $5!~/chr[CLMT]/) print $0}' /mnt/disk1/6/lxk/private/BL-Hi-C/K562/HiCPro_Analysis_BLHiC_K562_R1r1/hicpro_results/hic_results/data/BLHiC_K562_R1r1_val_trimLk3/BLHiC_K562_R1r1_val_trimLk3.allValidPairs > BLHiC_K562_R1r1_nochrCLMT.allValidPairs
shuf -n 822793 BLHiC_K562_R1r1_nochrCLMT.allValidPairs > BLHiC_K562_R1r1_nochrCLMT_shuf.allValidPairs
awk '{print $2"\t"$3"\t"$3+1"\n"$5"\t"$6"\t"$6+1}' BLHiC_K562_R1r1_nochrCLMT_shuf.allValidPairs | sort -k1,1 -k2,3n -S 9G >  BLHiC_K562_R1r1_nochrCLMT_shuf.singlepair
digest_genome.py -r GG^CC -o hg38_HaeIII.bed /ssd/index/bwa/hg38.fa
sort -k1,1 -k2,3n -S 9G hg38_HaeIII.bed > hg38_HaeIII_sort.bed
closestBed -a BLHiC_K562_R1r1_nochrCLMT_shuf.singlepair -b hg38_HaeIII_sort.bed -d -t all | awk '{if($10==0)print $4"\t"$5"\t"$6"\t"NR"\t"$7}' | sort -k1,1 -k2,3n -S 9G > BLHiC_K562.singlepair
closestBed -a BLHiC_K562.singlepair -b ../K562_hg38_motif_sort.bed -d -t all > BLHiC_K562_closest_motif.txt
sort -k1,1 -k2,2n -k3,3n -k4,4n -k6,6 -k7,7n -k8,8n -u --buffer-size=10% -T ./ --parallel=80 BLHiC_K562_closest_motif.txt > BLHiC_K562_closest_motif_u.txt
awk '{if($14==0) print $0}' BLHiC_K562_closest_motif_u.txt | groupBy -g 1,2,3,4 -c 7 -o min,max > BLHiC_K562_closest_motif_u_0_groupby.txt
