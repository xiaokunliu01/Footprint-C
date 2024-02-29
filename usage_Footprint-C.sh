#!/bin/bash

mkdir fastqc_R1r1/
mkdir trimAdapt_R1r1/
mkdir trimLk_R1r1/


#DNaseC_K562_WT
fastqc DNaseC_K562_WT_R1r1_1.fq.gz -o fastqc_R1r1/
fastqc DNaseC_K562_WT_R1r1_2.fq.gz -o fastqc_R1r1/

/mnt/disk2-0/lxk/private/software/miniconda3/bin/trim_galore -j 7 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired DNaseC_K562_WT_R1r1_*.fq.gz --gzip -o trimAdapt_R1r1/ --path_to_cutadapt /mnt/disk2-0/lxk/private/software/miniconda3/bin/cutadapt

# merge & S-linker number
mkdir merge_R1r1/
cd merge_R1r1/
/xch_lab/xch_lab/glx/share/script/flash ../trimAdapt_R1r1/DNaseC_K562_WT_R1r1_1_val_1.fq.gz ../trimAdapt_R1r1/DNaseC_K562_WT_R1r1_2_val_2.fq.gz -z -M 150 -t 10 -o DNaseC_K562_WT_R1r1 2>&1 | tee DNaseC_K562_WT_R1r1.log
perl ../DNaseC_K562_WT_merge_S-lnknum.pl
cd ..


/mnt/disk2-0/lxk/private/software/miniconda3/bin/cutadapt -a file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker.fa -A file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk1_1.fq.gz -p trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk1_2.fq.gz trimAdapt_R1r1/DNaseC_K562_WT_R1r1_1_val_1.fq.gz trimAdapt_R1r1/DNaseC_K562_WT_R1r1_2_val_2.fq.gz > trimLk_R1r1/WT_report1.txt
/mnt/disk2-0/lxk/private/software/miniconda3/bin/cutadapt -a file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker.fa -A file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk2_1.fq.gz -p trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk2_2.fq.gz trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk1_1.fq.gz trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk1_2.fq.gz > trimLk_R1r1/WT_report2.txt
/mnt/disk2-0/lxk/private/software/miniconda3/bin/cutadapt -a file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker_noAT.fa -A file:/xch_lab/xch_lab/lxk/private/optionData/MicroC_S-linker_noAT.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_1.fq.gz -p trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_2.fq.gz trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk2_1.fq.gz trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk2_2.fq.gz > trimLk_R1r1/WT_report3.txt

fastqc trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_1.fq.gz -o fastqc_R1r1/
fastqc trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_2.fq.gz -o fastqc_R1r1/

mkdir HiCPro_Analysis_WT_R1r1/
mkdir HiCPro_Analysis_WT_R1r1/data
mkdir HiCPro_Analysis_WT_R1r1/data/DNaseC_K562_WT_R1r1_val_trimLk3

cp trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_1.fq.gz HiCPro_Analysis_WT_R1r1/data/DNaseC_K562_WT_R1r1_val_trimLk3/
cp trimLk_R1r1/DNaseC_K562_WT_R1r1_val_trimLk3_2.fq.gz HiCPro_Analysis_WT_R1r1/data/DNaseC_K562_WT_R1r1_val_trimLk3/

cp /xch_lab/xch_lab/lxk/private/optionData/DNaseC_config_hicpro_XY_mm.txt HiCPro_Analysis_WT_R1r1/config_hicpro_WT.txt

cd HiCPro_Analysis_WT_R1r1/

hic-pro -c config_hicpro_WT.txt -i data -o hicpro_results -s mapping -s quality_checks
hic-pro -c config_hicpro_WT.txt -i hicpro_results/bowtie_results/bwt2 -o hicpro_results -s proc_hic -s merge_persample -s quality_checks

cd ..


rm fastqc_R1r1/*.zip



#rmdup

mkdir rmdup_R1r1/
cd rmdup_R1r1/

#no UMI
bamToBed -bedpe -i ../HiCPro_Analysis_WT_R1r1/hicpro_results/bowtie_results/bwt2/DNaseC_K562_WT_R1r1_val_trimLk3/DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs.bam | sort -k1,1 -k2,3n -k4,4 -k5,6n -S 9G | awk '{if($1==$4 && $2>$5){print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9}else{print $0}}' | sort -k1,6 -k9,10 -u -S 9G | awk -f /xch_lab/xch_lab/lxk/private/optionData/HiCpro_rmdup.awk | awk '{if($2==$5 && $3>$6){print $1 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "\t" $3 "\t" $4 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}else{print $0}}' | sort -k2,2V -k3,3n -k5,5V -k6,6n -S 9G > DNaseC_K562_WT.allValidPairs


#UMI
bamToBed -bedpe -i ../HiCPro_Analysis_WT_R1r1/hicpro_results/bowtie_results/bwt2/DNaseC_K562_WT_R1r1_val_trimLk3/DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs.bam | sort -k1,1 -k2,3n -k4,4 -k5,6n -S 9G | awk '{if($1==$4 && $2>$5){print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9}else{print $0}}' | awk -F';' '{print $1 "\t" $2$3$4$5$6$7$8$9}' | sort -k1,6 -k8,8 -k10,11 -u -S 9G | awk -f /xch_lab/xch_lab/lxk/private/optionData/HiCpro_rmdup_UMI.awk | awk '{if($2==$5 && $3>$6){print $1 "\t" $5 "\t" $6 "\t" $7 "\t" $2 "\t" $3 "\t" $4 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}else{print $0}}' | sort -k2,2V -k3,3n -k5,5V -k6,6n -S 9G > DNaseC_K562_WT_UMI.allValidPairs


cd ..

# optical duplication

mkdir optical_duplication_R1r1/
cd optical_duplication_R1r1/

samtools sort -@ 60 -m 1G ../HiCPro_Analysis_WT_R1r1/hicpro_results/bowtie_results/bwt2/DNaseC_K562_WT_R1r1_val_trimLk3/DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs.bam > DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort.bam


# DIS 2500
java -jar /xch_lab/xch_lab/lxk/public/software/picard.jar MarkDuplicates -REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT -I DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort.bam -O DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort_markdup.bam -M DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort.metrics --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --TAGGING_POLICY All > DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs.log 2>&1

# DIS 10000
java -jar /xch_lab/xch_lab/lxk/public/software/picard.jar MarkDuplicates -REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT -I DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort.bam -O DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort_markdup_DIS10k.bam -M DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_sort_DIS10k.metrics --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --TAGGING_POLICY All > DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs_DIS10k.log 2>&1

cd ..


# fragment_Length

mkdir fragment_Length_R1r1/
cd fragment_Length_R1r1/

bamToBed -bedpe -i ../HiCPro_Analysis_WT_R1r1/hicpro_results/bowtie_results/bwt2/DNaseC_K562_WT_R1r1_val_trimLk3/DNaseC_K562_WT_R1r1_val__hg38XY+mm10XY.bwt2pairs.bam | sort -k1,1 -k2,3n -k4,4 -k5,6n -S 9G | awk '{if($1==$4 && $2>$5){print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9}else{print $0}}' | awk -F';' '{print $1 "\t" $2$3$4$5$6$7$8$9}' | sort -k1,6 -k8,8 -k10,11 -u -S 9G > DNaseC_K562_WT_UMI_rmdup.pair
awk '{if($1==$4 && $3-$2>=$6-$5){print $3-$2 "\t" $6-$5} if($1==$4 && $3-$2<$6-$5){print $6-$5 "\t" $3-$2}}' DNaseC_K562_WT_UMI_rmdup.pair > DNaseC_K562_WT_UMI_rmdup_cis.length
awk '{if($1!=$4 && $3-$2>=$6-$5){print $3-$2 "\t" $6-$5} if($1!=$4 && $3-$2<$6-$5){print $6-$5 "\t" $3-$2}}' DNaseC_K562_WT_UMI_rmdup.pair > DNaseC_K562_WT_UMI_rmdup_trans.length

#V-plot
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' DNaseC_K562_WT_UMI_rmdup.pair | awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)+1"\t"NR"\t"$3-$2"\t"$10"\n"$4"\t"int(($5+$6)/2)"\t"int(($5+$6)/2)+1"\t"NR"_2""\t"$6-$5"\t"$11}' | sort -k1,1 -k2,2n -S 9G > DNaseC_K562_WT_UMI_sort.singlepair
closestBed -a DNaseC_K562_WT_UMI_sort.singlepair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d > DNaseC_K562_WT_UMI_sort_singlepair_closest_CTCFmotif.txt
awk -f /xch_lab/xch_lab/lxk/private/optionData/FragmentLengthVsDistance_stat.awk DNaseC_K562_WT_UMI_sort_singlepair_closest_CTCFmotif.txt > FragmentLengthVsDistance_K562_WT.csv


cd ..

# Dimerization

mkdir Dimerization_R1r1/
cd Dimerization_R1r1/

mkdir WT/ 
cd WT/

awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' ../../fragment_Length_R1r1/DNaseC_K562_WT_UMI_rmdup.pair | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$10}' | sort -k1,1 -k2,2n -S 9G > DNaseC_K562_WT_UMI_rmdup_sort.uppair
awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' ../../fragment_Length_R1r1/DNaseC_K562_WT_UMI_rmdup.pair | awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$11}' | sort -k1,1 -k2,2n -S 9G > DNaseC_K562_WT_UMI_rmdup_sort.downpair

closestBed -a DNaseC_K562_WT_UMI_rmdup_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > DNaseC_K562_WT_UMI_rmdup_sort_uppair_closest_CTCFmotif.txt
closestBed -a DNaseC_K562_WT_UMI_rmdup_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > DNaseC_K562_WT_UMI_rmdup_sort_downpair_closest_CTCFmotif.txt

sort -k4,4 -S 9G DNaseC_K562_WT_UMI_rmdup_sort_uppair_closest_CTCFmotif.txt > DNaseC_K562_WT_UMI_rmdup_sort_uppair_closest_CTCFmotif_sort.txt
sort -k4,4 -S 9G DNaseC_K562_WT_UMI_rmdup_sort_downpair_closest_CTCFmotif.txt > DNaseC_K562_WT_UMI_rmdup_sort_downpair_closest_CTCFmotif_sort.txt
paste DNaseC_K562_WT_UMI_rmdup_sort_uppair_closest_CTCFmotif_sort.txt DNaseC_K562_WT_UMI_rmdup_sort_downpair_closest_CTCFmotif_sort.txt | sort -k1,1 -k2,2n -k3,3n -S 9G > DNaseC_K562_WT_UMI_rmdup_pair_closet_CTCFmotif.txt

# awk '{if(($14!=-1 && $14<=10) && ($28<=10 && $28!=-1) && ($11=="-" && $25=="+") && ($5=="+" && $19=="-"))print $0}' DNaseC_K562_WT_UMI_rmdup_pair_closet_CTCFmotif.txt |  wc -l
awk -f /mnt/disk2-0/lxk/private/optionData/script/DNase-C/stat.awk DNaseC_K562_WT_UMI_rmdup_pair_closet_CTCFmotif.txt > DNaseC_K562_WT_UMI_rmdup_pair_closet_CTCFmotif.stat

cd ../../

# asChIP-seq_Analysis

mkdir asChIP-seq_Analysis_R1r1/
cd asChIP-seq_Analysis_R1r1/

awk '{print($1 "\t" $2 "\t" $3)}' ../fragment_Length_R1r1/DNaseC_K562_WT_UMI_rmdup.pair | awk '$1!~/chr[CLMT]/' > DNaseC_K562_WT_UMI_rmdup.bed
awk '{print($4 "\t" $5 "\t" $6)}' ../fragment_Length_R1r1/DNaseC_K562_WT_UMI_rmdup.pair | awk '$1!~/chr[CLMT]/' >> DNaseC_K562_WT_UMI_rmdup.bed
a=$(cat DNaseC_K562_WT_UMI_rmdup.bed | wc -l)
cat DNaseC_K562_WT_UMI_rmdup.bed | sort -k1,1 -k2,2n -k3,3n -S 9G | genomeCoverageBed -bg -i - -g /ssd/genome/hg38_chromsize.txt | awk -v a=$a '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/a*1000000}' | awk '{$4/=1;print}' OFS='\t' > DNaseC_K562_WT_UMI_rmdup.bdg
bedGraphToBigWig DNaseC_K562_WT_UMI_rmdup.bdg /ssd/genome/hg38_chromsize.txt DNaseC_K562_WT_UMI_rmdup.bw


mkdir plot_Average_Profile
cd plot_Average_Profile

computeMatrix reference-point --referencePoint center -S /mnt/disk4/public/K562/DNase-seq/K562_DNase.bw ../DNaseC_K562_WT_UMI_rmdup.bw -R /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o DNaseC_K562_WT_CTCFmotif.gz
plotProfile -m DNaseC_K562_WT_CTCFmotif.gz -o DNaseC_K562_WT_CTCFmotif_profile.png --outFileNameData DNaseC_K562_WT_CTCFmotif_profile.tab --plotHeight 20 --plotWidth 30 --perGroup

# computeMatrix reference-point --referencePoint center -S /mnt/disk4/public/K562/DNase-seq/K562_DNase.bw ../DNaseC_K562_WT_UMI_rmdup.bw -R /mnt/disk4/public/RefBed/refGene_hg38/hg38_TSS.bed -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o DNaseC_K562_WT_TSS.gz
# plotProfile -m DNaseC_K562_WT_TSS.gz -o DNaseC_K562_WT_TSS_profile.png --outFileNameData DNaseC_K562_WT_TSS_profile.tab --plotHeight 20 --plotWidth 30 --perGroup

# computeMatrix reference-point --referencePoint center -S /mnt/disk4/public/K562/DNase-seq/K562_DNase.bw ../DNaseC_K562_WT_UMI_rmdup.bw -R /mnt/disk4/public/RefBed/refGene_hg38/hg38_TTS.bed -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o DNaseC_K562_WT_TTS.gz
# plotProfile -m DNaseC_K562_WT_TTS.gz -o DNaseC_K562_WT_TTS_profile.png --outFileNameData DNaseC_K562_WT_TTS_profile.tab --plotHeight 20 --plotWidth 30 --perGroup

computeMatrix reference-point --referencePoint center -S /mnt/disk4/public/K562/DNase-seq/K562_DNase.bw ../DNaseC_K562_WT_UMI_rmdup.bw -R /mnt/disk4/public/K562/DNase-seq/K562_DNase_dis.summit -a 500 -b 501 -bs 1 --missingDataAsZero -p max -o DNaseC_K562_WT_DHS.gz
plotProfile -m DNaseC_K562_WT_DHS.gz -o DNaseC_K562_WT_DHS_profile.png --outFileNameData DNaseC_K562_WT_DHS_profile.tab --plotHeight 20 --plotWidth 30 --perGroup
