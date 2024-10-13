import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import cooler
import cooltools.lib.plotting

import cooltools
clr = cooler.Cooler('K562_MboI_600M.cool')
import bioframe
bins = clr.bins()[:]
hg38_genome = bioframe.load_fasta('/home/whh/private/DamID/hg38.fa');
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

view_df = pd.DataFrame({'chrom': clr.chromnames,
	'start': 0,
	'end': clr.chromsizes.values,
	'name': clr.chromnames}
	)


# obtain first 3 eigenvectors
cis_eigs = cooltools.eigs_cis(
	clr,
	gc_cov,
	view_df=view_df,
	n_eigs=3,
	)

# cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
eigenvector_track.to_csv('K562_MboI_600M.tsv',index=False,sep='\t')




import cooltools
clr = cooler.Cooler('K562_FootprintC_600M.cool')
import bioframe
bins = clr.bins()[:]
#hg38_genome = bioframe.load_fasta('/home/whh/private/linker_micro-c/20220420_ML/merge/hg38.fa');
#gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
#gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

view_df = pd.DataFrame({'chrom': clr.chromnames,
	'start': 0,
	'end': clr.chromsizes.values,
	'name': clr.chromnames}
	)


# obtain first 3 eigenvectors
cis_eigs = cooltools.eigs_cis(
	clr,
	gc_cov,
	view_df=view_df,
	n_eigs=3,
	)

# cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
eigenvector_track.to_csv('K562_FootprintC_600M.tsv',index=False,sep='\t')



import cooltools
clr = cooler.Cooler('K562_MicroC_600M.cool')
import bioframe
bins = clr.bins()[:]
#hg38_genome = bioframe.load_fasta('/home/whh/private/linker_micro-c/20220420_ML/merge/hg38.fa');
#gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
#gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

view_df = pd.DataFrame({'chrom': clr.chromnames,
	'start': 0,
	'end': clr.chromsizes.values,
	'name': clr.chromnames}
	)


# obtain first 3 eigenvectors
cis_eigs = cooltools.eigs_cis(
	clr,
	gc_cov,
	view_df=view_df,
	n_eigs=3,
	)

# cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
eigenvector_track.to_csv('K562_MicroC_600M.tsv',index=False,sep='\t')

