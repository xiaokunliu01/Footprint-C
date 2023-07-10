# awk 'BEGIN{srand()} {print rand()"\t"$0}' K562_DNaseC.allValidPairs | sort -nk 1 -S 50G -m | head -n600000000 | awk -F '\t' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' OFS='\t' > K562_DNaseC_600M.allValidPair
# /opt/HiC-Pro-3.0.0/bin/utils/hicpro2higlass.sh -i K562_DNaseC_600M.allValidPairs -r 5000 -c hg38_chromsize.txt -p 30
# cooler balance K562_DNaseC_600M.cool
# mustache -f K562_DNaseC_600M.cool -r 5000 --pThreshold 0.1 -p 20 -o K562_DNaseC_600M_5kb.bedpe
# stripenn compute --cool K562_MicroC_600M.cool --out K562_MicroC_600M/
# call-stripes -O K562_DNaseC_600M_stripes.bed -p K562_DNaseC_600M.cool --nproc 40
# chromosight detect -t12 --pattern=loops_small --min-dist=15000 --max-dist=2000000 K562_DNaseC_600M.cool K562_DNaseC_600M_loop.bed

shuf -n600000000 K562_DNaseC.allValidPairs > K562_DNaseC_600M.allValidPairs
/opt/HiC-Pro-3.0.0/bin/utils/hicpro2higlass.sh -i K562_DNaseC_600M.allValidPairs -r 5000 -c hg38_chromsize.txt -p 30
cooler balance K562_DNaseC_600M.cool
mustache -f K562_DNaseC_600M.cool -r 5000 --pThreshold 0.1 -p 20 -o K562_DNaseC_600M_5kb.bedpe
stripenn compute --cool K562_MicroC_600M.cool --out K562_MicroC_600M/
call-stripes -O K562_DNaseC_600M_stripes.bed -p K562_DNaseC_600M.cool --nproc 40
chromosight detect -t12 --pattern=loops_small --min-dist=15000 --max-dist=2000000 K562_DNaseC_600M.cool K562_DNaseC_600M_loop.bed
