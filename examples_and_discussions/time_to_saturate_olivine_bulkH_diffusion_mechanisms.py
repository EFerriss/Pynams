# -*- coding: utf-8 -*-
"""
Created on Wed May 04 16:39:25 2016

@author: Ferriss

how long to saturate a block
"""
from HydrogenCpx import cpx_spectra
from pynams import pynams, styles
from pynams import diffusivity_library as dlib


#%%

# Demouchy et al. 2016 3GPa experiments
SD1 = pynams.Sample(twoA_list=[1060.], twoB_list=[2060.], twoC_list=[1740.])
SD2 = pynams.Sample(twoA_list=[1060.], twoB_list=[2080.], twoC_list=[1770.])
SD3 = pynams.Sample(twoA_list=[1060.], twoB_list=[2020.], twoC_list=[1730.])
SD5 = pynams.Sample(twoA_list=[1010.], twoB_list=[1090.], twoC_list=[2090.])
SD6 = pynams.Sample(twoA_list=[1040.], twoB_list=[1180.], twoC_list=[2040.])
SD7 = pynams.Sample(twoA_list=[1710.], twoB_list=[2040.], twoC_list=[1040.])

sample = SD6
celsius = 900.
hours = 10.

fig, ax = dlib.fast_vs_slow(celsius=celsius, sample=sample,
                            minutes=60.*hours)
print 'Finished'


#%%
reload(pynams)
reload(styles)
reload(cpx_spectra)


testspec = pynams.Spectrum()

idx = 0
profile = cpx_spectra.profile_K4_904C_154hr_C
spectrum = profile.spectra_list[idx]

#fig, ax = spectrum.plot_peakfit()
label = ''.join(('observed\n', '{:.2f}'.format(profile.positions_microns[idx]), 
                 ' $\mu$m || ', profile.direction))
fig, ax = spectrum.plot_peakfit_and_baseline(bottom=0, label_spectrum=label)
fig.set_size_inches(6, 6)
