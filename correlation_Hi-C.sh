/opt/HiC-Pro_3.0.0/bin/utils/hicpro2higlass.sh -i /mnt/disk1/6/lxk/private/DNase-C/220725_DNaseC_K562_FA_EGS,FA_0percentSDS_4,8ulDNase_30m_noHI_rtaq_S,L-link_liga16C4h_exo1h_gDNA50-150/HiCPro_Analysis_FA_R1r4/hicpro_results/hic_results/data/DNaseC_K562_FA_R1r4_val_trimLk3/DNaseC_K562_FA_R1r4_val_trimLk3.allValidPairs -r 100000 -c /ssd/genome/hg38_chromsize_min.txt -n -p 80

/opt/HiC-Pro_3.0.0/bin/utils/hicpro2higlass.sh -i /mnt/disk1/6/lxk/private/DNase-C/221122_DNaseC_K562_FA_0percentSDS_8,16ulDNase_30m_noHI_S_link_gDNA_110-150/HiCPro_Analysis_8ulDNase_R1r3/hicpro_results/hic_results/data/DNaseC_K562_8ulDNase_R1r3_val_trimLk3/DNaseC_K562_8ulDNase_R1r3_val_trimLk3.allValidPairs -r 100000 -c /ssd/genome/hg38_chromsize_min.txt -n -p 80

/opt/HiC-Pro_3.0.0/bin/utils/hicpro2higlass.sh -i /mnt/disk1/6/lxk/private/in-situ-Hi-C/K562/fragment_Length/K562_MboI_R2.allValidPairs -r 100000 -c /ssd/genome/hg38_chromsize_min.txt -n -p 80

/opt/HiC-Pro_3.0.0/bin/utils/hicpro2higlass.sh -i /mnt/disk4/public/K562/Micro-C/K562_MicroC_R2r4.allValidPairs -r 100000 -c /ssd/genome/hg38_chromsize_min.txt -n -p 80


hicrep DNaseC_K562_8ulDNase_R1r3_val_trimLk3.mcool DNaseC_K562_FA_R1r4_val_trimLk3.mcool outputSCC_8ulVSFA.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000

hicrep DNaseC_K562_8ulDNase_R1r3_val_trimLk3.mcool K562_MboI_R2.mcool outputSCC_8ulVSinsituHiC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000

hicrep DNaseC_K562_8ulDNase_R1r3_val_trimLk3.mcool K562_MicroC_R2r4.mcool outputSCC_8ulVSMicroC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000

hicrep K562_MboI_R2.mcool K562_MicroC_R2r4.mcool outputSCC_insituHiCVSMicroC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000
