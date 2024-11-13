#/bin/bash

# 1. takes file: .allValidPairs file and hg38_chromsize_min.txt file
#    applies tool: hicpro2higlass.sh 
#    produces output: .mcool file

# 2. takes file: .mcool file
#    applies tool: hicpro2higlass.sh 
#    produces output: the stratum-adjusted correlation coefficient file


# generate .mcool file
for i in FootprintC_R1 FootprintC_R2 MicroC HiC
do
/opt/HiC-Pro_3.0.0/bin/utils/hicpro2higlass.sh -i $i'.'allValidPairs -r 100000 -c /ssd/genome/hg38_chromsize_min.txt -n -p 80
done

# calculate the stratum-adjusted correlation coefficient using HiCRep
# Footprint-C R1 VS Footprint-C R2
hicrep FootprintC_R1.mcool FootprintC_R2.mcool outputSCC_FootprintC_R1_VS_FootprintC_R2.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000
# FootprintC R1 VS in situ Hi-C
hicrep FootprintC_R1.mcool HiC.mcool outputSCC_FootprintC_R1_VS_HiC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000
# Footprint-C R1 VS Micro-C
hicrep FootprintC_R1.mcool MicroC.mcool outputSCC_FootprintC_R1_VS_MicroC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000
# in situ HiC VS Micro-C
hicrep HiC.mcool MicroC.mcool outputSCC_HiC_VS_MicroC.txt --h 1 --dBPMax 100000 --excludeChr "chrY" --binSize 100000
