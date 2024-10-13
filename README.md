# Footprint-C
We developed a protocol of Footprint-C, which yields high-resolution chromatin contact maps built upon intact and genuine footprints protected by transcription factor (TF) binding. When analyzed at one-dimensional level, the billions of chromatin contacts from Footprint-C enable genome-wide analysis at single footprint resolution, and reveal preferential modes of local TF cooperative binding. At pairwise contact level, Footprint-C exhibits higher efficiency in identifying chromatin structural features when compared with other Hi-C methods, segregates chromatin interactions emanating from adjacent TF footprints, and uncovers multiway interactions involving different TFs. Altogether, Footprint-C reveals that a rich regulatory lexicon of TF cooperativity underlies both local residence and distal chromatin interactions, in terms of TF identity, valency, and conformational configuration.

This repository contains the code and supporting data used to make a pipeline for analysing Footprint-C data. In this work, We use the human hg38 genome as reference, and we obtained the public datasets of K562 and HEK293T cell lines from Gene Expression Omnibus (GEO) or ENDCODE database.

## The following software prerequisites are required to run the pipeline

FastQC  
Trim Galore (v0.6.10)  
Cutadapt (v4.2)  
HiC-Pro (v3.0.0)  
Samtools (v1.16.1)  
bedtools (v2.30.0)  
deepTools (v3.5.0)  
Bowtie2 (v2.4.5)  
Cooler (v0.9.1)  
cooltools (v0.5.4)  
coolpup.py (v1.1.0)  
HiCPlotter  
MACS2 (v2.2.7.1)  
HOMER (v4.11)  
Mustache (v1.0.1)  
Chromosight  
Stripenn (v1.1.65.20)  
StripeCaller  
MEME (v5.4.1)  
HiCRep (https://github.com/dejunlin/hicrep)  
