#!/bin/bash

# calculate different TF-TF contact
# 1. takes file: .pair file and .motif file
#    applies tool: closestBed
#    produces output: different TF-TF contact and interaction distance


for i in CDC5L KLF11 KLF16 KMT2A MTF2 RBAK SP3 TFDP1 VEZF1 ZFP69B
do
mkdir $i
cd $i
closestBed -a /home/lxk/private/K562/FootprintC/Dimerization/all_pair/FootprintC_K562_FA_sort.uppair -b ../motif/K562_$i'.'motif -d -t first > uppair_closest_motif.txt
closestBed -a /home/lxk/private/K562/FootprintC/Dimerization/all_pair/FootprintC_K562_FA_sort.downpair -b ../motif/K562_$i'.'motif -d -t first > downpair_closest_motif.txt
sort -k4,4 --buffer-size=30% -T ./ --parallel=80 uppair_closest_motif.txt > uppair_closest_motif_sort.txt
rm uppair_closest_motif.txt
sort -k4,4 --buffer-size=30% -T ./ --parallel=80 downpair_closest_motif.txt > downpair_closest_motif_sort.txt
rm downpair_closest_motif.txt
 paste uppair_closest_motif_sort.txt downpair_closest_motif_sort.txt > FootprintC_K562_FA_pair_closet_motif.txt
rm uppair_closest_motif_sort.txt
rm downpair_closest_motif_sort.txt
awk -f /mnt/disk2-0/lxk/private/optionData/script/DNase-C/stat_noextend_wlink.awk FootprintC_K562_FA_pair_closet_motif.txt > FootprintC_K562_FA_pair_closet_motif.stat
cd ..
done

for i in CDC5L KLF11 KLF16 KMT2A MTF2 RBAK SP3 TFDP1 VEZF1 ZFP69B
do
cd $i
cut -f6-14 FootprintC_K562_FA_pair_closet_motif.txt > up.bed
cut -f20-28 FootprintC_K562_FA_pair_closet_motif.txt > down.bed
cd ..
done
#mkdir CTCF
#cd CTCF
#cut -f1-14 /home/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/CTCF/FootprintC_K562_FA_pair_closet_motif.txt > up.bed
#cut -f15-28 /home/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/CTCF/FootprintC_K562_FA_pair_closet_motif.txt > down.bed
#cd ..

#for i in E2F4 EGR1 KLF1 MAZ SP1 ZNF143
#do
#mkdir $i
#cd $i
#cut -f6-14 /home/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/$i/FootprintC_K562_FA_pair_closet_motif.txt > up.bed
#cut -f20-28 /home/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/$i/FootprintC_K562_FA_pair_closet_motif.txt > down.bed
#cd ..
#done


for i in CTCF CDC5L E2F4 EGR1 KLF1 KLF11 KLF16 KMT2A MAZ MTF2 RBAK SP1 SP3 TFDP1 VEZF1 ZFP69B ZNF143
do 
echo $i/up.bed >> tmp1.txt
echo $i/down.bed >> tmp2.txt
done

awk '{printf("%s ", $1)}' tmp1.txt
# CTCF/up.bed CDC5L/up.bed E2F4/up.bed EGR1/up.bed KLF1/up.bed KLF11/up.bed KLF16/up.bed KMT2A/up.bed MAZ/up.bed MTF2/up.bed RBAK/up.bed SP1/up.bed SP3/up.bed TFDP1/up.bed VEZF1/up.bed ZFP69B/up.bed ZNF143/up.bed
awk '{printf("%s ", $1)}' tmp2.txt
# CTCF/down.bed CDC5L/down.bed E2F4/down.bed EGR1/down.bed KLF1/down.bed KLF11/down.bed KLF16/down.bed KMT2A/down.bed MAZ/down.bed MTF2/down.bed RBAK/down.bed SP1/down.bed SP3/down.bed TFDP1/down.bed VEZF1/down.bed ZFP69B/down.bed ZNF143/down.bed
rm tmp1.txt tmp2.txt

paste CTCF/up.bed CDC5L/up.bed E2F4/up.bed EGR1/up.bed KLF1/up.bed KLF11/up.bed KLF16/up.bed KMT2A/up.bed MAZ/up.bed MTF2/up.bed RBAK/up.bed SP1/up.bed SP3/up.bed TFDP1/up.bed VEZF1/up.bed ZFP69B/up.bed ZNF143/up.bed CTCF/down.bed CDC5L/down.bed E2F4/down.bed EGR1/down.bed KLF1/down.bed KLF11/down.bed KLF16/down.bed KMT2A/down.bed MAZ/down.bed MTF2/down.bed RBAK/down.bed SP1/down.bed SP3/down.bed TFDP1/down.bed VEZF1/down.bed ZFP69B/down.bed ZNF143/down.bed > total.pair

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' /mnt/disk1/6/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/CTCF/FootprintC_K562_FA_pair_closet_motif.txt > CTCF_CTCF_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_CTCF_motif_pair_length60bp.txt > CTCF_CTCF_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_CTCF_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_CTCF_motif_pair_length60bp_cis_dis320bp_groupby.txt

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' /mnt/disk1/6/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/KLF1/FootprintC_K562_FA_pair_closet_motif.txt > KLF1_KLF1_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' KLF1_KLF1_motif_pair_length60bp.txt > KLF1_KLF1_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G KLF1_KLF1_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF1_KLF1_motif_pair_length60bp_cis_dis320bp_groupby.txt

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' /mnt/disk1/6/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/MAZ/FootprintC_K562_FA_pair_closet_motif.txt > MAZ_MAZ_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' MAZ_MAZ_motif_pair_length60bp.txt > MAZ_MAZ_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G MAZ_MAZ_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > MAZ_MAZ_motif_pair_length60bp_cis_dis320bp_groupby.txt

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' /mnt/disk1/6/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/ZNF143/FootprintC_K562_FA_pair_closet_motif.txt > ZNF143_ZNF143_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' ZNF143_ZNF143_motif_pair_length60bp.txt > ZNF143_ZNF143_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G ZNF143_ZNF143_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > ZNF143_ZNF143_motif_pair_length60bp_cis_dis320bp_groupby.txt

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' /mnt/disk1/6/lxk/private/K562/FootprintC/Dimerization/all_pair/230203/SP1/FootprintC_K562_FA_pair_closet_motif.txt > SP1_SP1_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' SP1_SP1_motif_pair_length60bp.txt > SP1_SP1_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G SP1_SP1_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > SP1_SP1_motif_pair_length60bp_cis_dis320bp_groupby.txt

awk '{if(($3-$2<=60) && ($17-$16<=60) && ($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' FootprintC_K562_FA_pair_closet_motif.txt > KLF11_KLF11_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' KLF11_KLF11_motif_pair_length60bp.txt > KLF11_KLF11_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G KLF11_KLF11_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF11_KLF11_motif_pair_length60bp_cis_dis320bp_groupby.txt


#CTCF-X

# # CDC5L
# awk '{if(($3-$2<=60) && ($161-$160<=60) && ($14==0) && ($181==0)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$173"\t"$174"\t"$175"\t"$176"\t"$177"\t"$178"\t"$179"\t"$180"\t"$5$163}' ../total.pair > CTCF_CDC5L_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_CDC5L_motif_pair_length60bp.txt > CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp_groupby.txt

# # E2F4
# awk '{if(($3-$2<=60) && ($161-$160<=60) && ($14==0) && ($190==0)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$182"\t"$183"\t"$184"\t"$185"\t"$186"\t"$187"\t"$188"\t"$189"\t"$5$163}' ../total.pair > CTCF_E2F4_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_E2F4_motif_pair_length60bp.txt > CTCF_E2F4_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_E2F4_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_E2F4_motif_pair_length60bp_cis_dis320bp_groupby.txt

# #EGR1
# awk '{if(($3-$2<=60) && ($161-$160<=60) && ($14==0) && ($199==0)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$191"\t"$192"\t"$193"\t"$194"\t"$195"\t"$196"\t"$197"\t"$198"\t"$5$163}' ../total.pair > CTCF_EGR1_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_EGR1_motif_pair_length60bp.txt > CTCF_EGR1_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_EGR1_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_EGR1_motif_pair_length60bp_cis_dis320bp_groupby.txt

# # KLF1
# awk '{if(($3-$2<=60) && ($161-$160<=60) && ($14==0) && ($208==0)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$200"\t"$201"\t"$202"\t"$203"\t"$204"\t"$205"\t"$206"\t"$207"\t"$5$163}' ../total.pair > CTCF_KLF1_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_KLF1_motif_pair_length60bp.txt > CTCF_KLF1_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_KLF1_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_KLF1_motif_pair_length60bp_cis_dis320bp_groupby.txt


for i in 181 190 199 208 217 226 244 253 262 271 280 289 298 307 316
# CDC5L E2F4 EGR1 KLF1 KLF11 KLF16 MAZ MTF2 RBAK SP1 SP3 TFDP1 VEZF1 ZFP69B ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($14==0) && ($a==0)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$5$163}' ../total.pair > CTCF_$i'_'motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_$i'_'motif_pair_length60bp.txt > CTCF_$i'_'motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G CTCF_$i'_'motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_$i'_'motif_pair_length60bp_cis_dis320bp_groupby.txt
done

# KLF1-X
for i in 199 217 226 244 280 289 316
# EGR1 KLF11 KLF16 MAZ SP3 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($50==0) && ($a==0)) print $42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$5$163}' ../total.pair > KLF1_$i'_'motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' KLF1_$i'_'motif_pair_length60bp.txt > KLF1_$i'_'motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G KLF1_$i'_'motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF1_$i'_'motif_pair_length60bp_cis_dis320bp_groupby.txt
done

# MAZ-X
for i in 199 217 226 280 289 316
# EGR1 KLF11 KLF16 SP3 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($86==0) && ($a==0)) print $78"\t"$79"\t"$80"\t"$81"\t"$82"\t"$83"\t"$84"\t"$85"\t"$(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$5$163}' ../total.pair > MAZ_$i'_'motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' MAZ_$i'_'motif_pair_length60bp.txt > MAZ_$i'_'motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G MAZ_$i'_'motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > MAZ_$i'_'motif_pair_length60bp_cis_dis320bp_groupby.txt
done

# SP3-X
for i in 199 217 226 289 316
# EGR1 KLF11 KLF16 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($122==0) && ($a==0)) print $114"\t"$115"\t"$116"\t"$117"\t"$118"\t"$119"\t"$120"\t"$121"\t"$(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$5$163}' ../total.pair > SP3_$i'_'motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' SP3_$i'_'motif_pair_length60bp.txt > SP3_$i'_'motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G SP3_$i'_'motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > SP3_$i'_'motif_pair_length60bp_cis_dis320bp_groupby.txt
done


# X-CTCF

# awk '{if(($3-$2<=60) && ($161-$160<=60) && ($23==0) && ($172==0)) print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$164"\t"$165"\t"$166"\t"$167"\t"$168"\t"$169"\t"$170"\t"$171"\t"$5$163}' ../total.pair > CDC5L_CTCF_motif_pair_length60bp.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CDC5L_CTCF_motif_pair_length60bp.txt > CDC5L_CTCF_motif_pair_length60bp_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CDC5L_CTCF_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CDC5L_CTCF_motif_pair_length60bp_cis_dis320bp_groupby.txt


for i in 23 32 41 50 59 68 86 95 104 113 122 131 140 149 158
# CDC5L E2F4 EGR1 KLF1 KLF11 KLF16 MAZ MTF2 RBAK SP1 SP3 TFDP1 VEZF1 ZFP69B ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($a==0) && ($172==0)) print $(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$164"\t"$165"\t"$166"\t"$167"\t"$168"\t"$169"\t"$170"\t"$171"\t"$5$163}' ../total.pair > $i'_'CTCF_motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' $i'_'CTCF_motif_pair_length60bp.txt > $i'_'CTCF_motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G $i'_'CTCF_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > $i'_'CTCF_motif_pair_length60bp_cis_dis320bp_groupby.txt
done

# X-KLF1

for i in 41 59 68 86 122 131 158
# EGR1 KLF11 KLF16 MAZ SP3 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($a==0) && ($208==0)) print $(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$200"\t"$201"\t"$202"\t"$203"\t"$204"\t"$205"\t"$206"\t"$207"\t"$5$163}' ../total.pair > $i'_'KLF1_motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' $i'_'KLF1_motif_pair_length60bp.txt > $i'_'KLF1_motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G $i'_'KLF1_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > $i'_'KLF1_motif_pair_length60bp_cis_dis320bp_groupby.txt
done


# X-MAZ

for i in 41 59 68 122 131 158
# EGR1 KLF11 KLF16 SP3 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($a==0) && ($244==0)) print $(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$236"\t"$237"\t"$238"\t"$239"\t"$240"\t"$241"\t"$242"\t"$243"\t"$5$163}' ../total.pair > $i'_'MAZ_motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' $i'_'MAZ_motif_pair_length60bp.txt > $i'_'MAZ_motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G $i'_'MAZ_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > $i'_'MAZ_motif_pair_length60bp_cis_dis320bp_groupby.txt
done

# X-SP3

for i in 41 59 68 131 158
# EGR1 KLF11 KLF16 TFDP1 ZNF143
do
awk -v a=$i '{if(($3-$2<=60) && ($161-$160<=60) && ($a==0) && ($280==0)) print $(a-8)"\t"$(a-7)"\t"$(a-6)"\t"$(a-5)"\t"$(a-4)"\t"$(a-3)"\t"$(a-2)"\t"$(a-1)"\t"$272"\t"$273"\t"$274"\t"$275"\t"$276"\t"$277"\t"$278"\t"$279"\t"$5$163}' ../total.pair > $i'_'SP3_motif_pair_length60bp.txt
awk '{if($1==$9 && $10-$2>320) print $0}' $i'_'SP3_motif_pair_length60bp.txt > $i'_'SP3_motif_pair_length60bp_cis_dis320bp.txt
sort -k4,4n -k12,12n -S 9G $i'_'SP3_motif_pair_length60bp_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > $i'_'SP3_motif_pair_length60bp_cis_dis320bp_groupby.txt
done



cat CDC5L_CTCF_motif_pair_length60bp_cis_dis320bp.txt CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp.txt > CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp_merge.txt; awk '{if($6=="+" && $14=="-"){i=i+1}if($6=="+" && $14=="+"){j=j+1}if($6=="-" && $14=="+"){m=m+1}if($6=="-" && $14=="-"){n=n+1}}END{print i+j+m+n"\ttotal\n"i"\tCon\n"j"\tTandem F\n"m"\tDiv\n"n"\tTandem R"}' CTCF_CDC5L_motif_pair_length60bp_cis_dis320bp_merge.txt


# distance_density
rm homo.tab
for i in CTCF KLF1 MAZ SP3 ZNF143; do awk -v a=$i '{print $10-$2"\t"a"-"a}' ../../$i/$i'_'$i'_'motif_pair_length60bp_cis_dis320bp.txt >> homo.tab; done



# CTCF-X
for i in KLF1 MAZ SP3 ZNF143
do
awk -v a=$i '{print $10-$2"\tCTCF-"a}' ../../CTCF/CTCF_$i/CTCF_$i'_'motif_pair_length60bp_cis_dis320bp_merge.txt >> C-X.tab
done

# X-X
awk '{print $10-$2"\tKLF1-MAZ"}' ../../KLF1/KLF1_MAZ/KLF1_MAZ_motif_pair_length60bp_cis_dis320bp_merge.txt > X-X.tab;awk '{print $10-$2"\tKLF1-SP3"}' ../../KLF1/KLF1_SP3/KLF1_SP3_motif_pair_length60bp_cis_dis320bp_merge.txt >> X-X.tab;awk '{print $10-$2"\tKLF1-ZNF143"}' ../../KLF1/KLF1_ZNF143/KLF1_ZNF143_motif_pair_length60bp_cis_dis320bp_merge.txt >> X-X.tab;awk '{print $10-$2"\tMAZ-SP3"}' ../../MAZ/MAZ_SP3/MAZ_SP3_motif_pair_length60bp_cis_dis320bp_merge.txt >> X-X.tab;awk '{print $10-$2"\tMAZ-ZNF143"}' ../../MAZ/MAZ_ZNF143/MAZ_ZNF143_motif_pair_length60bp_cis_dis320bp_merge.txt >> X-X.tab;awk '{print $10-$2"\tSP3-ZNF143"}' ../../SP3/SP3_ZNF143/SP3_ZNF143_motif_pair_length60bp_cis_dis320bp_merge.txt >> X-X.tab

