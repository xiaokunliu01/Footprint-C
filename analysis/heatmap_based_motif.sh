#!/bin/bash

awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$10}' /mnt/disk1/6/lxk/private/DNase-C/Total_FootprintC_K562/FootprintC_K562_total_UMI_wlink_nochrCLMT.pairs | sort -k1,1 -k2,2n --buffer-size=40% -T ./ --parallel=80 > FootprintC_K562_sort.uppair
awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$11}' /mnt/disk1/6/lxk/private/DNase-C/Total_FootprintC_K562/FootprintC_K562_total_UMI_wlink_nochrCLMT.pairs | sort -k1,1 -k2,2n --buffer-size=40% -T ./ --parallel=80 > FootprintC_K562_sort.downpair


closestBed -a FootprintC_K562_sort.uppair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > uppair_closest_motif.txt
closestBed -a FootprintC_K562_sort.downpair -b /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif -d -t first > downpair_closest_motif.txt
sort -k4,4 --buffer-size=40% -T ./ --parallel=80 uppair_closest_motif.txt > uppair_closest_motif_sort.txt
rm uppair_closest_motif.txt
sort -k4,4 --buffer-size=40% -T ./ --parallel=80 downpair_closest_motif.txt > downpair_closest_motif_sort.txt
rm downpair_closest_motif.txt
paste uppair_closest_motif_sort.txt downpair_closest_motif_sort.txt > FootprintC_K562_FA_pair_closet_motif.txt
rm uppair_closest_motif_sort.txt
rm downpair_closest_motif_sort.txt
awk -f /home/lxk/private/optionData/script/DNase-C/stat_noextend.awk FootprintC_K562_FA_pair_closet_motif.txt > FootprintC_K562_FA_pair_closet_motif.stat


awk '{if(($14==0) && ($28==0) && ($23!=$9)) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$5$19}' FootprintC_K562_FA_pair_closet_motif.txt > CTCF_CTCF_motif_pair.txt

awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_CTCF_motif_pair.txt > CTCF_CTCF_motif_pair_cis_dis320bp.txt; sort -k4,4 -k12,12 -S 9G CTCF_CTCF_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt


for i in KLF1 KLF16 MAZ KLF11 SP3
do
mkdir $i
cd $i
closestBed -a ../FootprintC_K562_sort.uppair -b /mnt/disk1/6/lxk/private/DNase-C/FootprintC_paper/dimer/motif/K562_$i'.'motif -d -t first > uppair_closest_motif.txt
closestBed -a ../FootprintC_K562_sort.downpair -b /mnt/disk1/6/lxk/private/DNase-C/FootprintC_paper/dimer/motif/K562_$i'.'motif -d -t first > downpair_closest_motif.txt
sort -k4,4 --buffer-size=25% -T ./ --parallel=80 uppair_closest_motif.txt > uppair_closest_motif_sort.txt
rm uppair_closest_motif.txt
sort -k4,4 --buffer-size=25% -T ./ --parallel=80 downpair_closest_motif.txt > downpair_closest_motif_sort.txt
rm downpair_closest_motif.txt
paste uppair_closest_motif_sort.txt downpair_closest_motif_sort.txt > FootprintC_K562_FA_pair_closet_motif.txt
awk -f /home/lxk/private/optionData/script/DNase-C/stat_noextend.awk FootprintC_K562_FA_pair_closet_motif.txt > FootprintC_K562_FA_pair_closet_motif.stat
cd ..
done


cd CTCF
cut -f1-14 FootprintC_K562_FA_pair_closet_motif.txt > up.bed
cut -f15-28 FootprintC_K562_FA_pair_closet_motif.txt > down.bed
cd ..

for i in KLF1 KLF16 MAZ KLF11 SP3
do
cd $i
cut -f6-14 FootprintC_K562_FA_pair_closet_motif.txt > up.bed
cut -f20-28 FootprintC_K562_FA_pair_closet_motif.txt > down.bed
cd ..
done

paste CTCF/up.bed KLF1/up.bed KLF11/up.bed KLF16/up.bed MAZ/up.bed SP3/up.bed CTCF/down.bed KLF1/down.bed KLF11/down.bed KLF16/down.bed MAZ/down.bed SP3/down.bed > total.pair


# CTCF-CTCF, KLF1, MAZ
awk '{if($14==0 && $23!=0 && $50!=0 && $73==0 && $82!=0 && $109!=0) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$65"\t"$66"\t"$67"\t"$68"\t"$69"\t"$70"\t"$71"\t"$72"\t"$5$64}' ../total.pair > CTCF_CTCF_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_CTCF_motif_pair.txt > CTCF_CTCF_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_CTCF_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt

awk '{if($14==0 && $23!=0 && $50!=0 && $73!=0 && $82==0 && $109!=0) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$74"\t"$75"\t"$76"\t"$77"\t"$78"\t"$79"\t"$80"\t"$81"\t"$5$64}' ../total.pair > CTCF_KLF1_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_KLF1_motif_pair.txt > CTCF_KLF1_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_KLF1_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_KLF1_motif_pair_cis_dis320bp_groupby.txt

awk '{if($14==0 && $23!=0 && $50!=0 && $73!=0 && $82!=0 && $109==0) print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$101"\t"$102"\t"$103"\t"$104"\t"$105"\t"$106"\t"$107"\t"$108"\t"$5$64}' ../total.pair > CTCF_MAZ_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' CTCF_MAZ_motif_pair.txt > CTCF_MAZ_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G CTCF_MAZ_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > CTCF_MAZ_motif_pair_cis_dis320bp_groupby.txt

# KLF1, MAZ-CTCF
awk '{if($14!=0 && $23==0 && $50!=0 && $73==0 && $82!=0 && $109!=0) print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$65"\t"$66"\t"$67"\t"$68"\t"$69"\t"$70"\t"$71"\t"$72"\t"$5$64}' ../total.pair > KLF1_CTCF_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' KLF1_CTCF_motif_pair.txt > KLF1_CTCF_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G KLF1_CTCF_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF1_CTCF_motif_pair_cis_dis320bp_groupby.txt

awk '{if($14!=0 && $23!=0 && $50==0 && $73==0 && $82!=0 && $109!=0) print $42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$65"\t"$66"\t"$67"\t"$68"\t"$69"\t"$70"\t"$71"\t"$72"\t"$5$64}' ../total.pair > MAZ_CTCF_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' MAZ_CTCF_motif_pair.txt > MAZ_CTCF_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G MAZ_CTCF_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > MAZ_CTCF_motif_pair_cis_dis320bp_groupby.txt

# KLF1-KLF1, MAZ
awk '{if($14!=0 && $23==0 && $50!=0 && $73!=0 && $82==0 && $109!=0) print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$74"\t"$75"\t"$76"\t"$77"\t"$78"\t"$79"\t"$80"\t"$81"\t"$5$64}' ../total.pair > KLF1_KLF1_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' KLF1_KLF1_motif_pair.txt > KLF1_KLF1_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G KLF1_KLF1_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF1_KLF1_motif_pair_cis_dis320bp_groupby.txt

awk '{if($14!=0 && $23==0 && $50!=0 && $73!=0 && $82!=0 && $109==0) print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$101"\t"$102"\t"$103"\t"$104"\t"$105"\t"$106"\t"$107"\t"$108"\t"$5$64}' ../total.pair > KLF1_MAZ_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' KLF1_MAZ_motif_pair.txt > KLF1_MAZ_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G KLF1_MAZ_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > KLF1_MAZ_motif_pair_cis_dis320bp_groupby.txt

# MAZ-KLF1
awk '{if($14!=0 && $23!=0 && $50==0 && $73!=0 && $82==0 && $109!=0) print $42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$74"\t"$75"\t"$76"\t"$77"\t"$78"\t"$79"\t"$80"\t"$81"\t"$5$64}' ../total.pair > MAZ_KLF1_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' MAZ_KLF1_motif_pair.txt > MAZ_KLF1_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G MAZ_KLF1_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > MAZ_KLF1_motif_pair_cis_dis320bp_groupby.txt

# MAZ-MAZ
awk '{if($14!=0 && $23!=0 && $50==0 && $73!=0 && $82!=0 && $109==0) print $42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$101"\t"$102"\t"$103"\t"$104"\t"$105"\t"$106"\t"$107"\t"$108"\t"$5$64}' ../total.pair > MAZ_MAZ_motif_pair.txt; awk '{if($1==$9 && $10-$2>320) print $0}' MAZ_MAZ_motif_pair.txt > MAZ_MAZ_motif_pair_cis_dis320bp.txt; sort -k4,4n -k12,12n -S 9G MAZ_MAZ_motif_pair_cis_dis320bp.txt | groupBy -g 1-16 -c 17,17 -o count,freqdesc > MAZ_MAZ_motif_pair_cis_dis320bp_groupby.txt





awk '{if($1=="chr22" && $9=="chr22"){if($6=="+" && $14=="-"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tFR"} if($6=="+" && $14=="+"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tFF"} if($6=="-" && $14=="+"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tRF"} if($6=="-" && $14=="-"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tRR"}}}' ../CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' > CTCF_motif_pair_chr22.txt

awk '{if($1=="chr22") print $0}' /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif | sed 's/v//g' > CTCF_chr22.motif

awk '{if($4>=154998 && $4<=155038) printf("%s ", $4)}' CTCF_chr22.motif

for i in 154998 155001 155007 155010 155012 155013 155015 155017 155018 155019 155020 155021 155024 155026 155028 155029 155030 155031 155032 155035 155036 155037; do for j in 154998 155001 155007 155010 155012 155013 155015 155017 155018 155019 155020 155021 155024 155026 155028 155029 155030 155031 155032 155035 155036 155037; do if [ $i == $j ]; then echo -e "$i;$j\t$i\t$j\tNA" >> index.txt; else echo -e "$i;$j\t$i\t$j\t0" >> index.txt; fi; done; done


join -a1 index.txt CTCF_motif_pair_chr22.txt |  awk '{if(NF==4){print $2"\t"$3"\t"$4"\tNA"}if(NF==8){print $5"\t"$6"\t"$7"\t"$8}}' > tmp_total.txt

awk '{if($4=="FR"){print $1"\t"$2"\t"$3} if($4!="FR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' tmp_total.txt > tmp_FR.txt

awk '{if($4=="FF"){print $1"\t"$2"\t"$3} if($4!="FF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' tmp_total.txt > tmp_FF.txt

awk '{if($4=="RR"){print $1"\t"$2"\t"$3} if($4!="RR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' tmp_total.txt > tmp_RR.txt

awk '{print $2"\t"$1"\t"$3}' tmp_FR.txt > tmp_FR_1.txt; cat tmp_FR.txt tmp_FR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > tmp_FR_2.txt





awk '{if($1==$9){if($6=="+" && $14=="-"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tFR"} if($6=="+" && $14=="+"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tFF"} if($6=="-" && $14=="+"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tRF"} if($6=="-" && $14=="-"){print $4";"$12"\t"$4"\t"$12"\t"$17"\tRR"}}}' CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' > CTCF_CTCF_motif_pair_index.txt

awk '{print $4}' /mnt/disk4/public/RefBed/CTCF/K562_CTCF_RAD21_DNase.motif | sed 's/v//g' > CTCF_motif_index.txt


for i in 12022 14306 19504 20426 20630 56574 85372 97180 114361 126651 127815 130259 152501 155395 155532 189790 201790 204691 243636 249676
do
mkdir $i
cd $i
a=$(awk -v a=$i '{if($1>=(a-20) && $1<=(a+20)) printf("%s ", $1)}' ../../CTCF_motif_index.txt)
for m in $a; do for n in $a; do echo -e "$m;$n\t$m\t$n\t0" >> index.txt; done; done
join -a1 index.txt ../../CTCF_CTCF_motif_pair_index_sort.txt |  awk '{if(NF==4){print $2"\t"$3"\t"$4}if(NF==8){print $5"\t"$6"\t"$7"\t"$8}}' > total.txt
awk '{if($1!=$2){print $1"\t"$2"\t"$3}if($1==$2){print $1"\t"$2"\tNA"}}' total.txt > total_1.txt
awk '{print $2"\t"$1"\t"$3}' total_1.txt > total_2.txt; cat total_1.txt total_2.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > total_3.txt
awk '{if($4=="FR"){print $1"\t"$2"\t"$3} if($4!="FR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FR.txt;awk '{if($4=="FF"){print $1"\t"$2"\t"$3} if($4!="FF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FF.txt;awk '{if($4=="RF"){print $1"\t"$2"\t"$3} if($4!="RF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RF.txt;awk '{if($4=="RR"){print $1"\t"$2"\t"$3} if($4!="RR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RR.txt
awk '{print $2"\t"$1"\t"$3}' FR.txt > FR_1.txt; cat FR.txt FR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FR_2.txt;awk '{print $2"\t"$1"\t"$3}' FF.txt > FF_1.txt; cat FF.txt FF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FF_2.txt;awk '{print $2"\t"$1"\t"$3}' RF.txt > RF_1.txt; cat RF.txt RF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RF_2.txt;awk '{print $2"\t"$1"\t"$3}' RR.txt > RR_1.txt; cat RR.txt RR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RR_2.txt
Rscript ../heatmap.R $i total_3; Rscript ../heatmap.R $i FR_2; Rscript ../heatmap.R $i FF_2; Rscript ../heatmap.R $i RF_2; Rscript ../heatmap.R $i RR_2
cd ..
done



for i in 4245 6885 7584 19000 22931 24852 39032 55782 90952 114239 124776 153643 155239 165342 207896 210699 224189 239410 246214 249605 251960
do
mkdir $i
cd $i
a=$(awk -v a=$i '{if($1>=(a-20) && $1<=(a+20)) printf("%s ", $1)}' ../../CTCF_motif_index.txt)
for m in $a; do for n in $a; do echo -e "$m;$n\t$m\t$n\t0" >> index.txt; done; done
join -a1 index.txt ../../CTCF_CTCF_motif_pair_index_sort.txt | awk '{if(NF==4){print $2"\t"$3"\t"$4}if(NF==8){print $5"\t"$6"\t"$7"\t"$8}}' > total.txt
awk '{if($1!=$2){print $1"\t"$2"\t"$3}if($1==$2){print $1"\t"$2"\tNA"}}' total.txt > total_1.txt
awk '{print $2"\t"$1"\t"$3}' total_1.txt > total_2.txt; cat total_1.txt total_2.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > total_3.txt
awk '{if($4=="FR"){print $1"\t"$2"\t"$3} if($4!="FR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FR.txt;awk '{if($4=="FF"){print $1"\t"$2"\t"$3} if($4!="FF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FF.txt;awk '{if($4=="RF"){print $1"\t"$2"\t"$3} if($4!="RF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RF.txt;awk '{if($4=="RR"){print $1"\t"$2"\t"$3} if($4!="RR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RR.txt
awk '{print $2"\t"$1"\t"$3}' FR.txt > FR_1.txt; cat FR.txt FR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FR_2.txt;awk '{print $2"\t"$1"\t"$3}' FF.txt > FF_1.txt; cat FF.txt FF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FF_2.txt;awk '{print $2"\t"$1"\t"$3}' RF.txt > RF_1.txt; cat RF.txt RF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RF_2.txt;awk '{print $2"\t"$1"\t"$3}' RR.txt > RR_1.txt; cat RR.txt RR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RR_2.txt
Rscript ../heatmap.R $i total_3; Rscript ../heatmap.R $i FR_2; Rscript ../heatmap.R $i FF_2; Rscript ../heatmap.R $i RF_2; Rscript ../heatmap.R $i RR_2
cd ..
done


# 3w 
# CTCF-MAZ-CTCF
awk '{print $0"\tM"NR}' ../../../dimer/motif/K562_CTCF.motif > K562_CTCF_rename.motif
awk '{print $0"\tM"NR}' ../../../dimer/motif/K562_MAZ.motif > K562_MAZ_rename.motif
cat K562_CTCF_rename.motif K562_MAZ_rename.motif | sort -k1,1 -k2,2n -k3,3n -S 9G | awk '{print $0"\t"NR}' > K562_CTCF_MAZ_rename.motif

awk '{print $4"\t"$9"\t"$10}' K562_CTCF_MAZ_rename.motif | grep C | sort -k1,1 -S 9G > CTCF_motif.index
awk '{print $4"\t"$9"\t"$10}' K562_CTCF_MAZ_rename.motif | grep M | sort -k1,1 -S 9G > MAZ_motif.index

awk '{print $4"\t"$12"\t"$0}' ../../CTCF/CTCF_KLF1_MAZ/CTCF_CTCF_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' | sort -k1,1 -S 9G | join - CTCF_motif.index | sort -k2,2 -S 9G | join -1 2 -2 1 - CTCF_motif.index | awk '{if($3==$11){if($8=="+" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFR"} if($8=="+" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFF"} if($8=="-" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRF"} if($8=="-" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRR"}}}' > CTCF-CTCF.count

awk '{print $4"\t"$12"\t"$0}' ../../CTCF/CTCF_KLF1_MAZ/CTCF_MAZ_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' | sort -k1,1 -S 9G | join - CTCF_motif.index | sort -k2,2 -S 9G | join -1 2 -2 1 - MAZ_motif.index | awk '{if($3==$11){if($8=="+" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFR"} if($8=="+" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFF"} if($8=="-" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRF"} if($8=="-" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRR"}}}' > CTCF-MAZ.count

awk '{print $4"\t"$12"\t"$0}' ../../CTCF/CTCF_KLF1_MAZ/MAZ_CTCF_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' | sort -k1,1 -S 9G | join - MAZ_motif.index | sort -k2,2 -S 9G | join -1 2 -2 1 - CTCF_motif.index | awk '{if($3==$11){if($8=="+" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFR"} if($8=="+" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFF"} if($8=="-" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRF"} if($8=="-" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRR"}}}' > MAZ-CTCF.count

awk '{print $4"\t"$12"\t"$0}' ../../MAZ/CTCF_KLF1_MAZ/MAZ_MAZ_motif_pair_cis_dis320bp_groupby.txt | sed 's/v//g' | sort -k1,1 -S 9G | join - MAZ_motif.index | sort -k2,2 -S 9G | join -1 2 -2 1 - MAZ_motif.index | awk '{if($3==$11){if($8=="+" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFR"} if($8=="+" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tFF"} if($8=="-" && $16=="+"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRF"} if($8=="-" && $16=="-"){print $22";"$24"\t"$22"\t"$24"\t"$19"\tRR"}}}' > MAZ-MAZ.count

cat CTCF-CTCF.count CTCF-MAZ.count MAZ-CTCF.count MAZ-MAZ.count | sort -k1,1 -S 9G > merge_CTCF-MAZ.count

for i in 121275
do
mkdir $i
cd $i
a=$(awk -v a=$i '{if($10>=(a-3) && $10<=(a+3)) printf("%s ", $10)}' ../K562_CTCF_MAZ_rename.motif)
for m in $a; do for n in $a; do echo -e "$m;$n\t$m\t$n\t0" >> index.txt; done; done
join -a1 index.txt ../merge_CTCF-MAZ.count | awk '{if(NF==4){print $2"\t"$3"\t"$4}if(NF==8){print $5"\t"$6"\t"$7"\t"$8}}' > total.txt
awk '{if($1!=$2){print $1"\t"$2"\t"$3}if($1==$2){print $1"\t"$2"\tNA"}}' total.txt > total_1.txt
awk '{print $2"\t"$1"\t"$3}' total_1.txt > total_2.txt; cat total_1.txt total_2.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > total_3.txt
awk '{if($4=="FR"){print $1"\t"$2"\t"$3} if($4!="FR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FR.txt;awk '{if($4=="FF"){print $1"\t"$2"\t"$3} if($4!="FF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > FF.txt;awk '{if($4=="RF"){print $1"\t"$2"\t"$3} if($4!="RF"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RF.txt;awk '{if($4=="RR"){print $1"\t"$2"\t"$3} if($4!="RR"){if($1!=$2){print $1"\t"$2"\t"0}if($1==$2){print $1"\t"$2"\tNA"}}}' total.txt > RR.txt
awk '{print $2"\t"$1"\t"$3}' FR.txt > FR_1.txt; cat FR.txt FR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FR_2.txt;awk '{print $2"\t"$1"\t"$3}' FF.txt > FF_1.txt; cat FF.txt FF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > FF_2.txt;awk '{print $2"\t"$1"\t"$3}' RF.txt > RF_1.txt; cat RF.txt RF_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RF_2.txt;awk '{print $2"\t"$1"\t"$3}' RR.txt > RR_1.txt; cat RR.txt RR_1.txt | sort -k1,1n -k2,2n -k3,3nr | sort -k1,1 -k2,2 -u > RR_2.txt
Rscript ../heatmap.R $i total_3; Rscript ../heatmap.R $i FR_2; Rscript ../heatmap.R $i FF_2; Rscript ../heatmap.R $i RF_2; Rscript ../heatmap.R $i RR_2
cd ..
done
