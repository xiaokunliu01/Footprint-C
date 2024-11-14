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
Read pairs were first trimmed by the Illumine adapter and bridge linker sequence (AGCCCGGTNNACGCCCGT, both forward and reverse complementary) from both ends by Trim Galore (https://github.com/FelixKrueger/TrimGalore) and Cutadapt (https://github.com/marcelm/cutadapt/). Only read pairs with bridge linker sequence detected and with both mates ≥10 bp after trimming were kept. Valid Footprint-C fragment contact pairs were obtained from the HiC-Pro analysis pipeline. The detailed description and code can be found at https://github.com/nservant/HiC-Pro. In brief, a pair of trimmed fastq files were mapped to the human (hg38) or Drosophila (dm6) genome separately by Bowtie2 with ```--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder``` mode. Aligned reads were paired by the read name. Pairs with multiple hits, low MAPQ (<=10), singleton, dangling end, and self-circle were removed. The fragment contact pairs were obtained from the output BAM files containing paired aligned reads by removing PCR duplicates and were used in downstream analyses.

Footprint-C datasets can be preprocessed using the ```Footprint-C_data_preprocessing/FootprintC_preprocessing.sh``` script.  
Place name FASTQ files as ```<example1>_1.fq.gz``` and ```<example1>_2.fq.gz```
```
bash FootprintC_preprocessing.sh <example1>
```

## Enrichment analysis
The fragment contact pairs of Footprint-C or other Hi-C datasets were first converted to single alignment BED files. Alignments of other Hi-C datasets were scaled to 60 bp, the median fragment length of Footprint-C. genomeCoverageBed in bedtools were used to compute normalized signal (reads per million for hg38). The average signal distribution relative to the occupied CTCF motifs and distal DNase I hypersensitive sites (dDHS) was computed across ±500 bp by computeMatrix in deepTools. The average signal of Read 1 or Read 2 alignments in HiCAR dataset around the occupied CTCF motifs and dDHS were calculated respectively. The list of functional CTCF motifs was generated as described previously. The sets of occupied CTCF motifs in K562 and GM12878 were screened by respective CTCF ChIP-seq datasets, and only those with an average RPM > 1 within a 30 bp window surrounding the motif center base were kept. The set of DHS in K562 or GM12878 cells was called from respective DNase-seq dataset using MACS2. DHS with q values > 20 and distance from TSS > 1 kb were defined dDHS.

The script is
```
analysis/enrichment_analysis.sh
``` 

## Loop calling
Loops were computed using Mustache or Chromosight. For Mustache, loops were called at 5-kb resolution using the options ```-r 5000 --pThreshold 0.1```. For Chromosight, loops were called by Chromosight detect function at 5-kb resolution. Loops smaller than 20 kb were removed. Common and specific loops between Footprint-C and Micro-C were characterized as described previously (PMID: 33239788). Loops with both anchors overlapping were called as common loops. Loops with one or neither anchor overlapping were called as specific loops. Anchors were extended ±20 kb, and pairtopair in bedtools were used to characterize common loops with parameters ```-type both -f 0.5```. Footprint-C or Micro-C specific loops were characterized with parameters ```-type notboth -slop 40000```.

## Stripe calling
Stripes were called from contact matrices using Stripenn or StripeCaller. The stripes were called by Stripenn at 5-kb resolution with the following settings: ```-m 0.95,0.96,0.97,0.98,0.99 -p 0.1```. Stripes were called by StripeCaller with the following settings: ```--local-num 2 --fold-enrichment 1.1 --min-seed-len 6```. Pile-up plots of stripes were obtained using coolpup.py, with the following settings: ```--local --rescale```. The loop domains in Supplementary Fig. 4k were annotated by horizontal and vertical stripes obtained by Stripenn, and were divided into four types of left-, right-, both-sided, and no stripe loops.

## Insulation score analysis
The insulation scores were calculated at 40-kb resolution using the cooltools package.

## Compartment analysis
The compartments were identified at 100-kb resolution using cooltools package. The eigenvector of the first principal component represents the compartment profile, with positive and negative values representing A and B compartments respectively.

## Motif uniqueness analysis
The fragment contact pairs of Footprint-C and other Hi-C datasets were first converted to single fragment BED files. The Footprint-C fragments less than or equal to 60 bp in length were kept for analysis. The fragments of in situ Hi-C were extended to the nearest GATC site (DpnII). The fragments of BL-Hi-C were extended to the nearest GGCC sites (HaeIII). The fragments of Hi-TrAC were extended to the nearest AATT or CATG sites (MluCI or NlaIII). The fragments of HiCAR were extended to the nearest GTAC sites (CviQI). The fragments of Micro-C were extended ±75 bp from the center base. Fragment coordinates were then intersected with all HOMER motif coordinates by intersectBed in bedtools. Finally, the proportions of fragments from Footprint-C, in situ Hi-C, BL-Hi-C, Hi-TrAC, HiCAR, or Micro-C datasets annotated with zero, one or multiple motifs were calculated.

## Analysis of motif orientation of CTCF-CTCF contacts
The closest CTCF motif to upstream and downstream read pair within contacts from Footprint-C or others Hi-C datasets was obtained by bedtools closestBed with parameters: -t first. Only the contacts with both alignments overlapping with center base of a CTCF motif were counted. The frequencies of CTCF contacts were calculated according to the orientations of the two motifs (+-: convergent, ++: tandem F, -+: divergent, --: tandem R).

## Construction of TF motif-based contact maps
The lists of functional CTCF or MAZ motifs were generated as described previously (PMID: 29590048). The motifs were sorted by genome coordinates, indexed, and used to construct a genome-wide 2D contact map. The Footprint-C fragments were annotated by the indexed motifs. The contact pairs with motifs annotated on both fragments were extracted and dumped to the respective bins in the 2D contact map. The counts in each bin were shaded by four different colors according to the orientations of the two motifs (+-, ++, -+, --). The motif-based contact maps were plotted using ggplot2, and were merged using Adobe Illustrator.


## Reproducibility analysis
The stratum-adjusted correlation coefficient using HiCRep73 was calculated for 100-kb resolution contact matrices by parameters settings of ```--h 1 --dBPMax 100000 --binSize 100000```. Comparisons were done between two biological replicates of Footprint-C libraries, or between Footprint-C and Micro-C or in situ Hi-C libraries.
