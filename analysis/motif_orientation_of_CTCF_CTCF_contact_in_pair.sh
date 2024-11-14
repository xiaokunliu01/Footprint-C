#!/bin/bash

# calculate frequency of motif orientation of CTCF-CTCF contacts from different Hi-C datasets
# 1. takes file: .pair file and .motif file
#    applies tool: closestBed and stat_noextent.awk script
#    produces output: .stat file of motif orientation


# FootprintC
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/ && $8!="no_adapterno_adapterno_adapter") print $0}' /mnt/disk1/6/lxk/private/FootprintC/220725_FootprintC_K562_FA_EGS,FA_0percentSDS_4,8ulDNase_30m_noHI_rtaq_S,L-link_liga16C4h_exo1h_gDNA50-150/fragment_Length_R1r4/FootprintC_K562_FA_UMI_rmdup.pair > FootprintC_K562_FA_UMI_rmdup_wlink_nochrCLMT.pair
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$10}' FootprintC_K562_FA_UMI_rmdup_wlink_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 1_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"NR"\t"$11}' FootprintC_K562_FA_UMI_rmdup_wlink_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 2_sort.downpair
closestBed -a 1_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF1+.motif -d -t first > 1_sort_uppair_closest_CTCFmotif.txt
closestBed -a 2_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF1+.motif -d -t first > 2_sort_downpair_closest_CTCFmotif.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 1_sort_uppair_closest_CTCFmotif.txt > 1_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 2_sort_downpair_closest_CTCFmotif.txt > 2_sort_downpair_closest_CTCFmotif_sort.txt
paste 1_sort_uppair_closest_CTCFmotif_sort.txt 2_sort_downpair_closest_CTCFmotif_sort.txt > pair_closet_CTCFmotif.txt
awk -f /home/lxk/private/optionData/script/FootprintC/stat_noextend.awk pair_closet_CTCFmotif.txt > pair_closet_CTCFmotif.stat


# Hi-TrAC
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$9}' ../../motif_unique/HiTrAC/Hi-TrAC_K562_R1r1.pair | sort -k1,1 -k2,2n -S 9G > 1_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"NR"\t"$10}' ../../motif_unique/HiTrAC/Hi-TrAC_K562_R1r1.pair | sort -k1,1 -k2,2n -S 9G > 2_sort.downpair
closestBed -a 1_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 1_sort_uppair_closest_CTCFmotif.txt
closestBed -a 2_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 2_sort_downpair_closest_CTCFmotif.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 1_sort_uppair_closest_CTCFmotif.txt > 1_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 2_sort_downpair_closest_CTCFmotif.txt > 2_sort_downpair_closest_CTCFmotif_sort.txt
paste 1_sort_uppair_closest_CTCFmotif_sort.txt 2_sort_downpair_closest_CTCFmotif_sort.txt > HiTrAC_pair_closet_CTCFmotif.txt
awk -f /home/lxk/private/optionData/script/FootprintC/stat_noextend.awk HiTrAC_pair_closet_CTCFmotif.txt > HiTrAC_pair_closet_CTCFmotif.stat

# Micro-C
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' /mnt/disk1/6/lxk/private/Micro-C/Public/K562/enrichment/K562_MicroC_R1r3_rmdup.pair > K562_MicroC_R1r3_rmdup_nochrCLMT.pair
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$9}' K562_MicroC_R1r3_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 1_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"NR"\t"$10}' K562_MicroC_R1r3_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 2_sort.downpair
closestBed -a 1_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 1_sort_uppair_closest_CTCFmotif.txt
closestBed -a 2_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 2_sort_downpair_closest_CTCFmotif.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 1_sort_uppair_closest_CTCFmotif.txt > 1_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 2_sort_downpair_closest_CTCFmotif.txt > 2_sort_downpair_closest_CTCFmotif_sort.txt
paste 1_sort_uppair_closest_CTCFmotif_sort.txt 2_sort_downpair_closest_CTCFmotif_sort.txt > pair_closet_CTCFmotif.txt
awk -f /home/lxk/private/optionData/script/FootprintC/stat_noextend.awk pair_closet_CTCFmotif.txt > pair_closet_CTCFmotif.stat


#in situ Hi-C
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' /mnt/disk1/6/lxk/private/in-situ-Hi-C/K562/fragment_Length_R1/insituHiC_K562_R2_UMI_rmdup.pair > insituHiC_K562_R2_UMI_rmdup_nochrCLMT.pair
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$9}' insituHiC_K562_R2_UMI_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 1_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"NR"\t"$10}' insituHiC_K562_R2_UMI_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 2_sort.downpair
closestBed -a 1_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 1_sort_uppair_closest_CTCFmotif.txt
closestBed -a 2_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > 2_sort_downpair_closest_CTCFmotif.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 1_sort_uppair_closest_CTCFmotif.txt > 1_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 2_sort_downpair_closest_CTCFmotif.txt > 2_sort_downpair_closest_CTCFmotif_sort.txt
paste 1_sort_uppair_closest_CTCFmotif_sort.txt 2_sort_downpair_closest_CTCFmotif_sort.txt > pair_closet_CTCFmotif.txt
awk -f /home/lxk/private/optionData/script/FootprintC/stat_noextend.awk pair_closet_CTCFmotif.txt > pair_closet_CTCFmotif.stat


# BL-Hi-C
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' /mnt/disk1/6/lxk/private/BL-Hi-C/K562/fragment_Length_R1r1/BLHiC_K562_R1r1_UMI_rmdup.pair > BLHiC_K562_R1r1_UMI_rmdup_nochrCLMT.pair  
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$9}' BLHiC_K562_R1r1_UMI_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 1_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"NR"\t"$10}' BLHiC_K562_R1r1_UMI_rmdup_nochrCLMT.pair | sort -k1,1 -k2,2n -S 9G > 2_sort.downpair
closestBed -a 1_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF1+.motif -d -t first > 1_sort_uppair_closest_CTCFmotif.txt
closestBed -a 2_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF1+.motif -d -t first > 2_sort_downpair_closest_CTCFmotif.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 1_sort_uppair_closest_CTCFmotif.txt > 1_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4n --buffer-size=10% -T ./ --parallel=80 2_sort_downpair_closest_CTCFmotif.txt > 2_sort_downpair_closest_CTCFmotif_sort.txt
paste 1_sort_uppair_closest_CTCFmotif_sort.txt 2_sort_downpair_closest_CTCFmotif_sort.txt > pair_closet_CTCFmotif.txt
awk -f /home/lxk/private/optionData/script/FootprintC/stat_noextend.awk pair_closet_CTCFmotif.txt > pair_closet_CTCFmotif.stat
