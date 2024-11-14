# Footprint-C
We developed a protocol of Footprint-C, which yields high-resolution chromatin contact maps built upon intact and genuine footprints protected by transcription factor (TF) binding. When analyzed at one-dimensional level, the billions of chromatin contacts from Footprint-C enable genome-wide analysis at single footprint resolution, and reveal preferential modes of local TF co-occupancy. At pairwise contact level, Footprint-C exhibits higher efficiency in identifying chromatin structural features when compared with other Hi-C methods, segregates chromatin interactions emanating from adjacent TF footprints, and uncovers multiway interactions involving different TFs. Altogether, Footprint-C results suggest that rich regulatory modes of TF may underly both local residence and distal chromatin interactions, in terms of TF identity, valency, and conformational configuration.

This repository contains the code and supporting data used to make a pipeline for analysing Footprint-C data. In this work, We use the human hg38 genome as reference, and we obtained the public datasets of K562 and HEK293T cell lines from Gene Expression Omnibus (GEO) or ENDCODE database. All raw and processed data generated in this study can be obtained at Genome Sequence Archive (GSA) under accession number: HRA004768 (https://ngdc.cncb.ac.cn/gsa-human/s/6t0JFuwn).

## The following software prerequisites are required to run the pipeline
- FastQC (v0.11.8) 
- Trim Galore (v0.6.10)  
- Cutadapt (v4.2)  
- HiC-Pro (v3.0.0)  
- Samtools (v1.16.1)  
- bedtools (v2.30.0)  
- deepTools (v3.5.0)  
- Bowtie2 (v2.4.5)  
- Cooler (v0.9.1)  
- cooltools (v0.5.4)  
- coolpup.py (v1.1.0)  
- HiCPlotter (https://github.com/jiangfuqing/HiCPlotter) 
- MACS2 (v2.2.7.1)  
- HOMER (v4.11)  
- Mustache (v1.0.1)  
- Chromosight (v1.6.3)  
- Stripenn (v1.1.65.20)  
- StripeCaller (v0.1.0) 
- MEME (v5.4.1)  
- HiCRep (https://github.com/dejunlin/hicrep)  

## General Instructions
This repository includes the following:
1. Test Footprint-C datasets and TF motif files;
2. Footprint-C datasets preprocessing scripts;
3. Scripts related to individual analysis.

## Footprint-C Data preprocessing
Read pairs were first trimmed by the Illumine adapter and bridge linker sequence (AGCCCGGTNNACGCCCGT, both forward and reverse complementary) from both ends by Trim Galore (https://github.com/FelixKrueger/TrimGalore) and Cutadapt (https://github.com/marcelm/cutadapt/). Only read pairs with bridge linker sequence detected and with both mates ≥10 bp after trimming were kept. Valid Footprint-C fragment contact pairs were obtained from the HiC-Pro54 analysis pipeline. The detailed description and code can be found at https://github.com/nservant/HiC-Pro. In brief, a pair of trimmed fastq files were mapped to the human (hg38) or Drosophila (dm6) genome separately by Bowtie255 with ‘--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder’ mode. Aligned reads were paired by the read name. Pairs with multiple hits, low MAPQ (<=10), singleton, dangling end, and self-circle were removed. The fragment contact pairs were obtained from the output BAM files containing paired aligned reads by removing PCR duplicates and were used in downstream analyses.

Footprint-C datasets can be preprocessed using the ```Footprint-C_data_preprocessing/FootprintC_preprocessing.sh``` script.  
Place name FASTQ files as ```<example1>_1.fq.gz``` and ```<example1>_2.fq.gz```
```
bash FootprintC_preprocessing.sh <example1>
```
## Footprint-C Data preprocessing
