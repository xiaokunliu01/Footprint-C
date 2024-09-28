#!/bin/bash


mkdir fastqc/
mkdir trimAdapt/
mkdir trimLk/

# fastqc
fastqc $1'_'1.fq.gz -o fastqc/
fastqc $1'_'2.fq.gz -o fastqc/

# trim adapter
/home/lxk/private/software/miniconda3/bin/trim_galore -j 7 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired $1'_'*.fq.gz --gzip -o trimAdapt/ --path_to_cutadapt /home/lxk/private/software/miniconda3/bin/cutadapt

# trim S-linker
/home/lxk/private/software/miniconda3/bin/cutadapt -a file:Footprint-C_S-linker_umi.fa -A file:Footprint-C_S-linker_umi.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk/$1'_'val_trimLk1_1.fq.gz -p trimLk/$1'_'val_trimLk1_2.fq.gz trimAdapt/$1'_'1_val_1.fq.gz trimAdapt/$1'_'2_val_2.fq.gz > trimLk/$1'_'report1.txt
/home/lxk/private/software/miniconda3/bin/cutadapt -a file:Footprint-C_S-linker_umi.fa -A file:Footprint-C_S-linker_umi.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk/$1'_'val_trimLk2_1.fq.gz -p trimLk/$1'_'val_trimLk2_2.fq.gz trimLk/$1'_'val_trimLk1_1.fq.gz trimLk/$1'_'val_trimLk1_2.fq.gz > trimLk/$1'_'report2.txt
/home/lxk/private/software/miniconda3/bin/cutadapt -a file:Footprint-C_S-linker_noAT_umi.fa -A file:Footprint-C_S-linker_noAT_umi.fa --rename='{id};{r1.adapter_name} {comment}' -j 10 --minimum-length=10 -o trimLk/$1'_'val_trimLk3_1.fq.gz -p trimLk/$1'_'val_trimLk3_2.fq.gz trimLk/$1'_'val_trimLk2_1.fq.gz trimLk/$1'_'val_trimLk2_2.fq.gz > trimLk/$1'_'report3.txt

# fastqc
fastqc trimLk/$1'_'trimLk3_1.fq.gz -o fastqc/
fastqc trimLk/$1'_'trimLk3_2.fq.gz -o fastqc/

# HiC-Pro
mkdir HiCPro_Analysis/
mkdir HiCPro_Analysis/data
mkdir HiCPro_Analysis/data/$1'_'val_trimLk3

cp trimLk/$1'_'val_trimLk3_1.fq.gz HiCPro_Analysis/data/$1'_'val_trimLk3/
cp trimLk/$1'_'val_trimLk3_2.fq.gz HiCPro_Analysis/data/$1'_'val_trimLk3/

cp /home/lxk/private/optionData/DNaseC_config_hicpro_XY.txt HiCPro_Analysis/config_hicpro.txt

cd HiCPro_Analysis/

hic-pro -c config_hicpro.txt -i data -o hicpro_results -s mapping -s quality_checks
hic-pro -c config_hicpro.txt -i hicpro_results/bowtie_results/bwt2 -o hicpro_results -s proc_hic -s merge_persample -s quality_checks

cd ..

rm fastqc/*.zip

# generate fragment contact pairs
mkdir fragment_contact_pairs/
cd fragment_contact_pairs/

bamToBed -bedpe -i ../HiCPro_Analysis/hicpro_results/bowtie_results/bwt2/$1'_'val_trimLk3/$1'_'val__hg38XY.bwt2pairs.bam | sort -k1,1 -k2,3n -k4,4 -k5,6n -S 9G | awk '{if($1==$4 && $2>$5){print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9}else{print $0}}' | awk -F';' '{print $1 "\t" $2$3$4$5$6$7$8$9}' | sort -k1,6 -k8,8 -k10,11 -u -S 9G | awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/ && $8!="no_adapterno_adapterno_adapter") print $0}' > $1'.'pair

cd ..


# fragment_Length

mkdir fragment_Length/
cd fragment_Length/

awk '{print $3-$2"\t"$6-$5}' ../fragment_contact_pairs/$1'.'pair > $1'.'length
