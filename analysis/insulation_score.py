import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cooler
import cooltools.lib.plotting
from cooltools import insulation
import cooltools


resolution = 40000
clr = cooler.Cooler('../loop/DNase-C/K562_DNaseC_60M.mcool::/resolutions/40000')
windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
insulation_table = insulation(clr, windows, verbose=True)
first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
#insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]
insulation_table.to_csv('DNase-C_insulation_table_600M.bed',index=False,sep='\t')


resolution = 40000
clr = cooler.Cooler('../loop_new/K562_MicroC_600M.mcool::/resolutions/40000')
windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
insulation_table = insulation(clr, windows, verbose=True)
first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
#insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]
insulation_table.to_csv('Micro-C_insulation_table_600M.bed',index=False,sep='\t')


resolution = 40000
clr = cooler.Cooler('../loop_new/K562_MboI_600M.mcool::/resolutions/40000')
windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
insulation_table = insulation(clr, windows, verbose=True)
first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
#insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]
insulation_table.to_csv('Hi-C_insulation_table_600M.bed',index=False,sep='\t')

