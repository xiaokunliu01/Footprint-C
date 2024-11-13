#!/bin/bash

# 1. takes file: fragment contact pairs file and CTCF motifs file, 
     applies tools: closestBed, 
     produces output: v-plot data file.



# plot V-plot
mkdir vplot/
cd vplot/

awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)+1"\t"NR"\t"$3-$2"\t"$10"\n"$4"\t"int(($5+$6)/2)"\t"int(($5+$6)/2)+1"\t"NR"_2""\t"$6-$5"\t"$11}' ../fragment_contact_pairs/$1'.'pair | sort -k1,1 -k2,2n -S 9G > $1'.'bed

closestBed -a $1'.'bed -b motif/K562_CTCF_RAD21_DNase.motif -d -t first > $1'_'closest_CTCFmotif.txt

awk -f /home/lxk/private/optionData/vplot.awk $1'_'closest_CTCFmotif.txt > FragmentLengthVsDistance.csv

rm $1'.'bed
rm $1'_'closest_CTCFmotif.txt


# cd ..

# # Dimerization

# mkdir Dimerization_R1r1/
# cd Dimerization_R1r1/

# mkdir K562_noRNase/ 
# cd K562_noRNase/

# awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$10}' ../../fragment_Length_R1r1/DNaseC_K562_noRNase_UMI_rmdup.pair | sort -k1,1 -k2,2n -S 9G > DNaseC_K562_noRNase_UMI_rmdup_sort.uppair
# awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$11}' ../../fragment_Length_R1r1/DNaseC_K562_noRNase_UMI_rmdup.pair | sort -k1,1 -k2,2n -S 9G > DNaseC_K562_noRNase_UMI_rmdup_sort.downpair

# closestBed -a DNaseC_K562_noRNase_UMI_rmdup_sort.uppair -b /home/lxk/private/optionData/function_motif/K562_CTCF_RAD21_DNase.motif -d -t first > DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif.txt
# closestBed -a DNaseC_K562_noRNase_UMI_rmdup_sort.downpair -b /home/lxk/private/optionData/function_motif/K562_CTCF_RAD21_DNase.motif -d -t first > DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif.txt

# rm DNaseC_K562_noRNase_UMI_rmdup_sort.uppair
# rm DNaseC_K562_noRNase_UMI_rmdup_sort.downpair

# sort -k4,4 -S 9G DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif.txt > DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif_sort.txt
# sort -k4,4 -S 9G DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif.txt > DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif_sort.txt

# rm DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif.txt
# rm DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif.txt

# paste DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif_sort.txt DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif_sort.txt | sort -k1,1 -k2,2n -k3,3n -S 9G > DNaseC_K562_noRNase_UMI_rmdup_pair_closet_CTCFmotif.txt

# rm DNaseC_K562_noRNase_UMI_rmdup_sort_uppair_closest_CTCFmotif_sort.txt
# rm DNaseC_K562_noRNase_UMI_rmdup_sort_downpair_closest_CTCFmotif_sort.txt

# awk -f /home/lxk/private/optionData/script/DNase-C/stat_noextend.awk DNaseC_K562_noRNase_UMI_rmdup_pair_closet_CTCFmotif.txt > DNaseC_K562_noRNase_UMI_rmdup_pair_closet_CTCFmotif.stat

# cd ../../

# # asChIP-seq_Analysis

# mkdir asChIP-seq_Analysis_R1r1/
# cd asChIP-seq_Analysis_R1r1/

# awk '{print($1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6)}' ../fragment_Length_R1r1/DNaseC_K562_noRNase_UMI_rmdup.pair > DNaseC_K562_noRNase_UMI_rmdup.bed
# a=$(cat DNaseC_K562_noRNase_UMI_rmdup.bed | wc -l)
# sort -k1,1 -k2,2n -k3,3n -S 9G DNaseC_K562_noRNase_UMI_rmdup.bed | genomeCoverageBed -bg -i - -g /ssd/genome/hg38_chromsize.txt | awk -v a=$a '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/a*1000000}' | awk '{$4/=1;print}' OFS='\t' > DNaseC_K562_noRNase_UMI_rmdup.bdg
# bedGraphToBigWig DNaseC_K562_noRNase_UMI_rmdup.bdg /ssd/genome/hg38_chromsize.txt DNaseC_K562_noRNase_UMI_rmdup.bw

# rm DNaseC_K562_noRNase_UMI_rmdup.bed
# rm DNaseC_K562_noRNase_UMI_rmdup.bdg


# # Mm
# awk '{print($1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6)}' ../fragment_Length_R1r1/K562_noRNase_mm_intra.pair > K562_noRNase_mm.bed
# a=$(cat K562_noRNase_mm.bed | wc -l)
# sort -k1,1 -k2,2n -k3,3n -S 9G K562_noRNase_mm.bed | genomeCoverageBed -bg -i - -g /ssd/genome/mm10_chromsize.txt | awk -v a=$a '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/a*1000000}' | awk '{$4/=1;print}' OFS='\t' > K562_noRNase_mm.bdg
# bedGraphToBigWig K562_noRNase_mm.bdg /ssd/genome/mm10_chromsize.txt K562_noRNase_mm.bw

# rm K562_noRNase_mm.bed
# rm K562_noRNase_mm.bdg


# mkdir plot_Average_Profile
# cd plot_Average_Profile

# computeMatrix reference-point --referencePoint center -S ../DNaseC_K562_noRNase_UMI_rmdup.bw -R /home/lxk/private/optionData/function_motif/K562_CTCF_RAD21_DNase.motif -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o DNaseC_K562_noRNase_CTCFmotif.gz
# plotProfile -m DNaseC_K562_noRNase_CTCFmotif.gz -o DNaseC_K562_noRNase_CTCFmotif_profile.png --outFileNameData DNaseC_K562_noRNase_CTCFmotif_profile.tab --plotHeight 6 --plotWidth 8 --perGroup --refPointLabel 'CTCF motif' --yMin 0 -z 'no RNase'
