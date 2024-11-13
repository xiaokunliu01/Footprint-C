#/bin/bash

# 1. takes file: .allValidPairs file and hg38_chromsize.txt file
#    applies tool: hicpro2higlass.sh and cooler
#    produces output: .mcool file

# 2. takes file: .mcool file
#    applies tool: Mustache or Chromosight
#    produces output: loops file

# 3. takes file: .mcool file
#    applies tool: Stripenn or StripeCaller
#    produces output: stripes file


# generate .cool file
# shuf -n600000000 K562_$i.allValidPairs > K562_$i_600M.allValidPairs
for i in FootprintC MicroC HiC
do
/opt/HiC-Pro-3.0.0/bin/utils/hicpro2higlass.sh -i K562_$i'.'allValidPairs -r 5000 -c hg38_chromsize.txt -p 30
cooler balance K562_$i'.'cool
done

# call loop
# Mustache or Chromosight methods
for i in FootprintC MicroC HiC
do
mustache -f K562_$i'.'cool -r 5000 --pThreshold 0.1 -p 20 -o K562_$i'_'5kb.bedpe
chromosight detect -t12 --pattern=loops_small --min-dist=15000 --max-dist=2000000 K562_$i'_'600M.cool K562_$i_600M'_'loop.bed
done

# call stripes
# Stripenn or StripeCaller methods
for i in FootprintC MicroC HiC
do
mkdir $i
stripenn compute --cool K562_$i'.'cool --out $i/
call-stripes -O K562_$i'_'stripes.bed -p K562_$i'.'cool --nproc 40
done



