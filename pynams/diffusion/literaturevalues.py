# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:49:59 2015

@author: Ferriss

Library of diffusivities, both data and fit lines, and 
some functions to play with them. The data in in a CSV in this folder.

The focus is on hydrogen diffusion in olivine and clinopyroxene.
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import pynams
from pynams import Sample
from pynams.diffusion.diffusivities import Diffusivities
import pandas as pd
import os

# pull data from spreadsheet
try:
    file = os.path.join(pynams.__path__[0], 'diffusion', 
                        'literaturevalues.csv')
except FileNotFoundError:
    print('Warning: Could not locate and import literature data')
    exit()
    
olivine = pd.read_csv(file)
olivine = olivine.dropna(how='all') 
olivine.fillna(0, inplace=True)
olivine["paper"] = olivine["Author"] + ' ' + olivine["Year"].map(str)
olivinegroups = olivine.groupby(['paper', 'mechanism', 
                                 'percentpv', 'orientation'])
    
## Uncomment to see group names
#for x, group in olivinegroups:
#    print(x)
#print()

SanCarlosOlivine = Sample(Fe3=0.003, Fe2=0.174, Mg=1.821, Ti=0., Al=0.)

# Mackwell and Kohlstedt 1990
SanCarlosOlivineKM = Sample(Mg=0.907, Fe=0.090)
MK90 = Diffusivities(description='Mackwell & Kohlstedt 1990',
                     sample=SanCarlosOlivineKM)
MK90.fill_in_data(df=olivinegroups, mech='bulk', percentpv=0)
MK90.solve_Ea_D0(printout=False)
MK90.activation_energy_kJmol[2] = MK90.activation_energy_kJmol[0] # from paper
MK90.D0_m2s[2] = 5E-6 # from paper

# Kohlstedt and Mackwell 1998
KM98_fast = Diffusivities(description='Kohlstedt & Mackwell 1998',
                          sample=SanCarlosOlivineKM)
KM98_fast.fill_in_data(df=olivinegroups, mech='bulk', percentpv=0)
KM98_fast.solve_Ea_D0(printout=False) # estimate D0
KM98_fast.activation_energy_kJmol = [145., 180., 110., 0] # from paper

KM98_slow = Diffusivities(description='Kohlstedt & Mackwell 1998',
                          sample=SanCarlosOlivineKM)
KM98_slow.fill_in_data(df=olivinegroups, mech='bulk', percentpv=100)
KM98_slow.solve_Ea_D0(printout=False)

# Demouchy and Mackwell 2003
DM03 = Diffusivities(description='Demouchy & Mackwell 2003',
                     sample=Sample(Fe2=0., Fe3=0., Al=0., Mg=2.0, Ti=0.))
DM03.fill_in_data(df=olivinegroups, mech='bulk', percentpv=100)

# Demouchy and Mackwell 2006
SanCarlosOlivineDM = Sample(Mg=0.904, Fe=0.092)
DM06_fast = Diffusivities(description = 'Demouchy & Mackwell 2006', 
                          sample = SanCarlosOlivineDM)
DM06_fast.fill_in_data(df=olivinegroups, mech='bulk', percentpv=0)

DM06_slow = Diffusivities(description = 'Demouchy & Mackwell 2006', 
                          sample = SanCarlosOlivineDM)
DM06_slow.fill_in_data(df=olivinegroups, mech='bulk', percentpv=100)

# pp and pv mechanisms in olivine
pp = Diffusivities(description='KM98 + DM06 "pp"')
pp.log10D = [d98+d06 for d98, d06 in \
             zip(KM98_fast.log10D, DM06_fast.log10D)]
pp.celsius = [d98+d06 for d98, d06 in \
             zip(KM98_fast.celsius, DM06_fast.celsius)]
pp.solve_Ea_D0(printout=False)

pv = Diffusivities(description='KM98 + DM03 + DM06 "pv"')
pv.log10D = [d98 + d03 + d06 for d98, d03, d06 in \
             zip(KM98_slow.log10D, DM03.log10D,DM06_slow.log10D)]
pv.celsius = [d98 + d03 + d06 for d98, d03, d06 in \
             zip(KM98_slow.celsius, DM03.celsius, DM06_slow.celsius)]
pv.solve_Ea_D0(printout=False) # quite different from DM06 results
pv.D0_m2s = [10**-4.5, 10**-4.5, 10**-1.4, 0] #DM06 pg 8
pv.activation_energy_kJmol = [204, 204, 258, 0] #DM06 pg 8


def mix_olivine_mechanisms(percent_slow, celsius):
    """
    Determines diffusivities assuming that some percentage
    of H diffuses by the fast 'proton-polaron' mechanism, and
    some diffuses by the slow 'proton-vacancy' mechanism.
    
    Required input:
        percentage of H that is moving by the slow mechanism
        temperature in Celsius
    
    Output: list of three log10 diffusivities, || a, b, and c in m2/s that 
    represent a linear mixture of the (unlogged) fast and slow mechanism
    diffusivities.
    """
    Ds_log = pv.whatIsD(celsius, printout=False)[0:3]
    Df_log = pp.whatIsD(celsius, printout=False)[0:3]
    Ds = [(10**logD) * percent_slow/100. for logD in Ds_log]
    Df = [(10**logD) * (100.-percent_slow)/100. for logD in Df_log]
    D = np.array(Ds) + np.array(Df)
    D = np.log10(D)
    return list(D)

# Padron-Navarta et al. 2014
# [Si]: 3613 peak in MgO-buffered Fo, 
# [Si]-[Ti]: 3613 peak in Ti-doped sample
# [Ti]: 3525 peak
# [Mg]: 3220 peak
pnav_Mg = Diffusivities(description = 'Padron-Navarta et al. 2014')
pnav_Mg.fill_in_data(df=olivinegroups, mech='[Mg]', percentpv=100)

pnav_Si = Diffusivities(description = 'Padron-Navarta et al. 2014')
pnav_Si.fill_in_data(df=olivinegroups, mech='[Si]', percentpv=100)

pnav_Ti = Diffusivities(description = 'Padron-Navarta et al. 2014')
pnav_Ti.fill_in_data(df=olivinegroups, mech='[Ti]', percentpv=100)

# Demouchy et al. 2016 hydration at 3GPa
D16 = Diffusivities(description='Demouchy et al. 2016',
                    sample=Sample(Fe=0.092, Mg=0.903)) # Ni=0.005, trace Cr

#%% The rest of this has data in an old format that I haven't had a 
### need to switch over into the CSV but I don't want to discard

##Wanamaker_Si = Diffusivities(
##               description = "San Carlos ol. $V''''_{Si}$\nWanamaker 1994",
##               celsius_u= [1100., 1200.], 
##               logDu = [-10.4, -9.99],
##               basestyle = {'color' : 'k', 'marker' : '*', 'alpha' : 0.5,
##                            'linestyle' : 'none', 'markersize' : 10})
##
##Wanamaker_Mg = Diffusivities(
##               description = "San Carlos ol. $V''_{Me}$\nWanamaker 1994",
##               celsius_u = [1100., 1200., 1300.],
#               logDu = [-11.20, -10.8, -10.3],
#               basestyle = {'color' : 'g', 'marker' : '*', 'alpha' : 0.5,
#                            'fillstyle' : 'none',
#                            'linestyle' : 'none', 'markersize' : 10})
#

##%% clinopyroxene
#H_CrDiopside_Ingrin95 = Diffusivities(
#                        description = ''.join(('H in Cr-rich diopside in air',
#                               '\nmainly 3645 cm$^{-1}$\nIngrin et al. 1995')),
#                        celsius_all = np.array([973, 1073, 1173, 1273]) - 273.,
#                        logDx = [np.log10(3E-14), np.log10(9E-14),
#                                 np.log10(5.5E-13), np.log10(1.9E-12)],
#                        logDz = [np.log10(5E-14), np.log10(13E-14), 
#                                np.log10(4.5E-13), np.log10(12E-13)],
#                        basestyle = {'color' : 'green', 'marker' : '^', 
#                                   'markersize' :  12, 
#                                   'linestyle' : 'none', 'alpha' : 1.})
#
#H_cpx_Wade2008 = Diffusivities(
#                 description = 'H in cpx phenocryst\nWade et al. 2008',
#                 celsius_all = [1100], logDu = [-13.],
#                 basestyle = {'color' : 'indigo', 'marker' : '$\clubsuit$', 
#                            'markersize' :  12,
#                            'markeredgewidth' : 0,
#                            'linestyle' : 'none', 'alpha' : 1.})
#
#H_diopside_Woods00 = Diffusivities(
#                     description = 'Jaipur bulk H\nWoods et al. 2000',
#                     logDx = np.log10(np.array([5e-13, 7e-12, 1.1e-11, 1.5e-11, 
#                                          2e-11, 1.8e-11, 3.5e-12, 3.5e-11])),
#                     logDy = np.log10(np.array([1.5e-12, 3e-12, 2.5e-12, 8e-13, 
#                                                2e-12])),
#                     logDz = np.log10(np.array([2e-12, 7e-12, 1.5e-11, 1.5e-11, 
#                                                6e-12, 3e-11])),
#                     celsius_x = [700, 750, 800, 800, 800, 850, 700, 850],
#                     celsius_y = [750, 800, 800, 800, 850],
#                     celsius_z = [700, 750, 800, 800, 800, 850],
#                     basestyle = {'color' : 'k', 'marker' : 's', 
#                                  'markerfacecolor' : 'turquoise',
#                                  'markersize' :  10, 
#                                  'linestyle' : 'none', 'alpha' : 1.,
#                                  'markeredgewidth' : 0.5})
#                            
#H_cpxBasanite_Xia00 = Diffusivities(
#                      description = ''.join(('H in basanite cpx\n3630 + 3500',
#                            ' cm$^{-1}$\nf$_{02}=10^{-14}$ Xia et al. 2000')),
#                      logDz = [np.log10(1.8E-12), np.log10(6.5E-12)],
#                      celsius_z = [850, 950],
#                      basestyle = {'color' : 'k', 'marker' : '>', 
#                                   'markerfacecolor' : 'k',
#                                   'markersize' :  8, 
#                                   'linestyle' : 'none', 'alpha' : 1.,
#                                   'markeredgewidth' : 1})
#
#H_diopside_noFe = Diffusivities(
#                  description = ''.join(('pure synth. diopside in air',
#                                         '\nSundvall et al. 2009')),
#                  logDx = [-12.3, -14.7], 
#                  celsius_x = np.array([1273, 1073]) - 273.,
#                  logDy = [-12.6, -13.9, -15.1],
#                  celsius_y = np.array([1273, 1173, 1073]) - 273.,
#                  basestyle = {'color' : 'b', 'marker' : '>', 
#                               'markerfacecolor' : 'b', 'markersize' :  8, 
#                               'linestyle' : 'none', 'alpha' : 1.,
#                               'markeredgewidth' : 1})
#
#H_diopside_Sundvall = Diffusivities(
#                      description = ''.join(('synth. Fe-bearing diopside\n',
#                                             'in air; Sundvall et al. 2009')),
#                      logDy = [-13.7, -15.3, -15.9],
#                      celsius_y = np.array([1000., 900., 800.]),
#                      basestyle = {'color' : 'y', 'marker' : 'v', 
#                                   'markerfacecolor' : 'y', 'markersize' :  13, 
#                                   'linestyle' : 'none', 'alpha' : 1.,
#                                   'markeredgewidth' : 1})
#
#Fe_diopside = Diffusivities()
#Fe_diopside.description = 'Fe in diopside\nAzough & Freer 2000'
#Fe_diopside.celsius_u = np.array([950., 1100.])
#Fe_diopside.activation_energy_kJmol = 161.5
#Fe_diopside.D0_m2s = 6.22E-15
##Fe_diopside.logDu = Fe_diopside.makelines_Ea_D0(Ea=161.5, log10D0=10.**(6.22E-15), 
##                                celsius_list=Fe_diopside.celsius_u)
#Fe_diopside.basestyle = {'color' : 'k', 'marker' : 'None', 
#                         'linestyle' : '-', 'alpha' : 1.}
#
#
## H in Kunlun diopside
## 1000 C, in whole-block paper, label: K5
##Dx_wholeblock = ufloat(-13.2, 0.2)
##Dy_wholeblock = ufloat(-13.4, 0.2)
##Dz_wholeblock = ufloat(-13.6, 0.3)
##Dy_slice_FTIR = ufloat(-13.1, 0.3)
##Dz_slice_FTIR = ufloat(-13.1, 0.2)
##Dy_slice_SIMS = ufloat(-13.3, 0.4)
##Dz_slice_SIMS = ufloat(-13.2, 0.4)
##Dy_ave = np.mean([Dy_slice_FTIR, Dy_slice_SIMS, Dy_wholeblock])
##Dz_ave = np.mean([Dz_slice_FTIR, Dz_slice_SIMS, Dz_wholeblock])
##
##
##Kunlun_bulkH = Diffusivities()
##Kunlun_bulkH.description = 'Kunlun diopside\nbulk H, QFM'
##Kunlun_bulkH.celsius_all = [1000.]
##Kunlun_bulkH.logDx = [-13.0]
##Kunlun_bulkH.logDy = [-13.5]
##Kunlun_bulkH.logDz = [-13.5]
##Kunlun_bulkH.basestyle = {'color' : 'black', 'marker' : 'D', 
##                          'markersize' :  markersizefloat, 
##                          'linestyle' : 'none', 'alpha' : 0.5,}
##Kunlun_bulkH.logDx_error = [0.1]
##Kunlun_bulkH.logDy_error = [0.1]
##Kunlun_bulkH.logDz_error = [0.1]
##
##Kunlun_peak3617 = Diffusivities()
##Kunlun_peak3617.description = 'Kunlun diopside\nPeak at 3617 cm$^{-1}$'
##Kunlun_peak3617.celsius_all = [1000.]
##Kunlun_peak3617.logDx = [-13.4]
##Kunlun_peak3617.logDy = [-12.3]
##Kunlun_peak3617.logDz = [-13.5]
##Kunlun_peak3617.logDx_error = [0.4]
##Kunlun_peak3617.logDy_error = [0.03]
##Kunlun_peak3617.logDz_error = [0.2]
##Kunlun_peak3617.basestyle = dict(Kunlun_bulkH.basestyle.items())
##Kunlun_peak3617.basestyle['color'] = 'red'
##
##Kunlun_peak3540 = Diffusivities()
##Kunlun_peak3540.description = 'Kunlun diopside\nPeak at 3540 cm$^{-1}$'
##Kunlun_peak3540.celsius_all = [1000.]
##Kunlun_peak3540.logDx = [-12.7]
##Kunlun_peak3540.logDy = [-12.2]
##Kunlun_peak3540.logDz = [-12.7]
##Kunlun_peak3540.logDx_error = [0.3]
##Kunlun_peak3540.logDy_error = [0.1]
##Kunlun_peak3540.logDz_error = [0.3]
##Kunlun_peak3540.basestyle = dict(Kunlun_bulkH.basestyle.items())
##Kunlun_peak3540.basestyle['color'] = 'blue'
#
## Jaipur diopside, Woods et al. 2001
#Jaipur_bulk = Diffusivities()
#Jaipur_bulk.description = 'Jaipur diopside bulk H\n(Woods et al. 2001)'
#Jaipur_bulk.celsius_x = [700, 750, 800, 800, 800, 850, 700, 850]
#Jaipur_bulk.celsius_y = [750, 800, 800, 800, 850]
#Jaipur_bulk.celsius_z = [700, 750, 800, 800, 800, 850]
#Jaipur_bulk.logDx = np.log10([5e-13, 7e-12, 1.1e-11, 1.5e-11, 2e-11, 1.8e-11, 
#                     3.5e-12, 3.5e-11])
#Jaipur_bulk.logDy = np.log10([1.5e-12, 3e-12, 2.5e-12, 8e-13, 2e-12])
#Jaipur_bulk.logDz = np.log10([2e-12, 7e-12, 1.5e-11, 1.5e-11, 6e-12, 3e-11])
#Jaipur_bulk.basestyle = {'color' : 'green', 'marker' : 's',
#                         'markersize' : markersizefloat, 
#                         'linestyle' : 'none'}
#
#O_diopside_self = Diffusivities()
#O_diopside_self.description = 'O self-diffusion in diopside\nRyerson & McKeegan 1994'
#O_diopside_self.logDz = [np.log10(1.96e-21), np.log10(2.04e-21), np.log10(1.48e-21),
#                         np.log10(2.29e-21), np.log10(1.19e-20), np.log10(2.47e-20),
#                         np.log10(3.63e-20), np.log10(1.9e-20), np.log10(3.17e-20),
#                         np.log10(1.75e-20), np.log10(2.17e-20), np.log10(3.74e-20),
#                         np.log10(8.94e-20), np.log10(9.64e-20), np.log10(1.22e-19)]
#O_diopside_self.celsius_z = [1104, 1104, 1104, 1105, 1150, 1200, 1200, 1200, 
#                             1200, 1202, 1202, 1202, 1251, 1251, 1251]
#O_diopside_self.basestyle = {'color' : 'red', 'marker' : '^', 
#                           'linestyle' : 'none',
#                           'fillstyle' : 'right', 'alpha' : 0.5}                             
#
#U_diopside = Diffusivities()
#U_diopside.description = 'U in diopside\nVan Orman et al. 1998'
#U_diopside.celsius_z = [1150, 1150, 1200, 1200, 1200, 1300, 1300]  
#U_diopside.logDz = [np.log10(4.54e-22), np.log10(9.44e-22), np.log10(2.90e-21),
#                    np.log10(3.95e-21), np.log10(2.35e-21), np.log10(2.35e-21),
#                    np.log10(2.18e-20)]
#U_diopside.basestyle = {'color' : 'orange', 'marker' : 'o', 
#                           'linestyle' : 'none', 'markersize' : 7, 
#                           'fillstyle' : 'right', 'alpha' : 0.5}                             
#
#Th_diopside = Diffusivities()
#Th_diopside.description = 'Th in diopside\nVan Orman et al. 1998'
#Th_diopside.celsius_z = [1150, 1150, 1200, 1200, 1200, 1300, 1300]  
#Th_diopside.logDz = [np.log10(1.59e-21), np.log10(1.2e-21), np.log10(5.26e-21),
#                     np.log10(3.80e-21), np.log10(3.22e-21), np.log10(2.12e-20),
#                     np.log10(2.75e-20)]
#Th_diopside.basestyle = {'color' : 'yellow', 'marker' : 'o', 
#                           'linestyle' : 'none', 'markersize' : 7, 
#                           'fillstyle' : 'right', 'alpha' : 0.5}                             
#
#Ce_diopside = Diffusivities()
#Ce_diopside.description = 'Ce in diopside\nVan Orman et al. 2001'
#Ce_diopside.celsius_z = [1300, 1275, 1250, 1250, 1225, 1200, 1200, 1200, 1200, 
#                         1175, 1150]
#Ce_diopside.logDz = [np.log10(31.9e-21), np.log10(25.0e-21), np.log10(11.5e-21),
#                     np.log10(6.83e-21), np.log10(5.83e-21), np.log10(2.53e-21),
#                     np.log10(4.45e-21), np.log10(4.01e-21), np.log10(4.53e-21),
#                     np.log10(0.68e-21), np.log10(0.62e-21)]
#Ce_diopside.basestyle = {'color' : 'lime', 'marker' : 'o', 
#                           'linestyle' : 'none', 'markersize' : 7, 
#                           'fillstyle' : 'right', 'alpha' : 0.5}                             
#
#Al_diopside = Diffusivities()
#Al_diopside.description = 'Al in diopside; Sautter\net al. EPSL 1988'
#Al_diopside.logDu = [-18.495]
##Al_diopside.logDu_error = [0.66]
#Al_diopside.celsius_u = [1180.]
#Al_diopside.basestyle = {'color' : 'g', 'marker' : 'x',
#                         'fillstyle' : 'full'}
#
#CaMg_diopside_2010 = Diffusivities()
#CaMg_diopside_2010.description = 'Ca-Mg interdiffusion in di.\nZhang et al. 2010'
#CaMg_diopside_2010.logDx = [-19.33, -19.37, -19.92, -19.55, -19.99, -19.97, 
#                            -20.39, -20.41, -20.82, -21.18]
#CaMg_diopside_2010.celsius_x = [1150.00, 1150.00, 1100.00, 1100.00, 1050.00, 
#                                1050.00, 1000.00, 1000.00, 950.00, 950.00]
#CaMg_diopside_2010.logDy = [-19.39, -19.38, -20.06, -20.07, -20.95, -21.08, 
#                            -20.92, -21.18, -21.26, -21.36]
#CaMg_diopside_2010.celsius_y = [1150.00, 1150.00, 1050.00, 1050.00, 1000.00, 
#                                1000.00, 1000.00, 1000.00, 950.00, 950.00]
#CaMg_diopside_2010.logDz = [-19.49, -19.44, -20.12, -20.24, -21.36, -21.38]
#CaMg_diopside_2010.celsius_z = [1150.00, 1150.00, 1050.00, 1050.00, 950.00, 950.00]
#                            
#CaMg_diopside_2010.basestyle = {'color' : 'orange', 'marker' : 's', 
#                                'linestyle' : 'none', 'markersize' : 7,
#                                'fillstyle' : 'full', 'alpha' : 0.8}
#
#CaMg_diopside_1983 = Diffusivities()
#CaMg_diopside_1983.description = 'Ca-Mg interdiffusion in cpx\nBrady & McCallister 1983'
#CaMg_diopside_1983.logDu = np.array([np.log10(1.4E-16), np.log10(6.9E-17), np.log10(5.6E-17),
#                                 np.log10(2.0E-17), np.log10(5.7E-18), np.log10(5.6E-16),
#                                 np.log10(2.8E-16), np.log10(1.7E-16), np.log10(6.9E-17),
#                                 np.log10(8.3E-16), np.log10(4.2E-16), np.log10(5.6E-16)]) - 2.
#CaMg_diopside_1983.celsius_all = [1100., 1100., 1100., 1100., 1100., 1150., 1150., 
#                             1150., 1150., 1200., 1200., 1250]
#CaMg_diopside_1983.basestyle = {'color' : 'orange', 'marker' : 's', 
#                           'linestyle' : 'none', 'markersize' : 9,
#                           'fillstyle' : 'none', 'alpha' : 0.5}
#                           
#FeMg_cpx_2013 = Diffusivities()
#FeMg_cpx_2013.description = 'Fe-Mg interdiffusion in cpx\nMueller et al. 2013'
#FeMg_cpx_2013.logDz = [-21, -20.8, -20.46, -18, -21.89, -21.05, -18.64, -19.75, 
#                       -18.3, -17.52, -19.52, -20.4, -20.46, -19.41, -18.7, 
#                       -20.4, -19.7, -18.7, -19.92, -19.52, -19.82, -19.92]
#FeMg_cpx_2013.celsius_z = [850, 900, 950, 1150, 800, 905, 1106, 1006, 1154, 
#                           1200, 1035, 924, 956, 1048, 1102, 945, 999, 1100, 
#                           1000, 1007, 1007, 1007]
#FeMg_cpx_2013.basestyle = {'color' : 'grey', 'marker' : 'h', 
#                           'linestyle' : 'none', 'markersize' : 7,
#                           'fillstyle' : 'right', 'alpha' : 0.8}
#                         
#FeMg_diopside = Diffusivities()
#FeMg_diopside.description = 'Fe-Mg interdiffusion\nDimanov & Wiedenbeck 2006'
#FeMg_diopside.logDz = np.array([-16.6649, -14.9891, -14.5297, -14.6098, -15.0794, 
#                       -15.6405, -15.9297, -15.8744, -16.139, -15.224, 
#                       -16.0353, -14.1226]) - 2. # original reported in log D (cm2/s)
#FeMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
#                           1050, 1150, 1100, 1100]
#FeMg_diopside.basestyle = {'color' : 'green', 'marker' : 'o', 
#                           'linestyle' : 'none',
#                           'fillstyle' : 'right', 'alpha' : 0.5}
#
#MnMg_diopside = Diffusivities()
#MnMg_diopside.description = 'Mn-Mg interdiffusion\nDimanov & Wiedenbeck 2006'
#MnMg_diopside.logDz = np.array([-16.6923, -15.0946, -14.8291, -14.8714, -15.223, 
#                       -15.711, -15.9744, -16.0584, -16.2665, -15.8388, 
#                       -16.2668, -14.3497]) - 2. # original reported in log D (cm2/s)
#MnMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
#                           1050, 1150, 1100, 1100]
#MnMg_diopside.basestyle = {'color' : 'yellow', 'marker' : 'o', 
#                           'linestyle' : 'none',
#                           'fillstyle' : 'right', 'alpha' : 0.5}
#
#FeMnMg_diopside = Diffusivities()
#FeMnMg_diopside.description = '(Fe,Mn)Mg interdiffusion in di.\nDimanov & Wiedenbeck 2006'
#FeMnMg_diopside.logDz = np.array([-16.6784, -15.0387, -14.6541, -14.7212, -15.1453, 
#                         -15.6743, -15.9515, -15.9567, -16.1981, -15.4307, 
#                         -16.1358, -14.2215]) - 2. # original reported in log D (cm2/s)
#FeMnMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
#                           1050, 1150, 1100, 1100]
#FeMnMg_diopside.basestyle = {'color' : 'g', 'marker' : 's', 
#                           'linestyle' : 'none',
#                           'fillstyle' : 'right', 'alpha' : 0.5}
#
#Ti_diopside = Diffusivities()
#Ti_diopside.description = 'Ti in diopside\nCherniak & Liang 2012'
#Ti_diopside.logDz = [-19.89, -19.90, -20.25, -20.43, -20.98, -20.92, -21.41, 
#                     -21.05, -21.43, -21.58, -22.22, -22.27, -22.24, -22.77]
#Ti_diopside.celsius_z = [1250, 1200, 1200, 1151, 1102, 1102, 1090, 1052, 999,
#                         1000, 952, 954, 905, 905]
#Ti_diopside.basestyle = {'color' : 'b', 'marker' : 'o', 
#                           'linestyle' : 'none', 'markersize' : 7, 
#                           'fillstyle' : 'right', 'alpha' : 0.5}
#He_cpx = Diffusivities()
#He_cpx.description = 'He in cpx\nTrull & Kurz 1993'
#He_cpx.logDu = np.array([np.log10(1.26E-10), np.log10(1.34E-10), 
#                                   np.log10(1.32E-10), np.log10(1.17E-10),
#                                   np.log10(7.24E-10), np.log10(9.29E-10), 
#                                   np.log10(8.04E-10), np.log10(6.52E-10),
#                                   np.log10(1.47E-8), np.log10(1.48E-8), 
#                                   np.log10(9.19E-9), np.log10(5.06E-9), 
#                                   np.log10(4.12E-9)]) - 2.
#He_cpx.celsius_u = [965., 965., 965., 965., 1070., 1070., 1070., 1070.,
#                             1170., 1170., 1170., 1170., 1170.]
#He_cpx.basestyle = {'color' : 'lawngreen', 'marker' : '$\spadesuit$', 
#                    'linestyle' : 'none', 
#                    'markersize' : 10, 'fillstyle' : 'full', 'alpha' : 1.}
#
#Li_cpx_interstitial = Diffusivities()
#Li_cpx_interstitial.description = 'Li in cpx as interstitial\nRichter et al. 2014'
#Li_cpx_interstitial.logDu = np.log10(np.array([9.2E-10, 1.6E-9, 1.6E-9, 
#                                     3.3E-8, 3.7E-8, 1.2E-8, 4.6E-10, 6.4E-12]))
#Li_cpx_interstitial.celsius_u = [900.]*len(Li_cpx_interstitial.logDu)
#Li_cpx_interstitial.basestyle = {'color' : 'darkorchid', 'marker' : '+', 
#                                 'linestyle' : 'none', 
#                                 'markersize' : 10, 'alpha' : 1.}
#                             
#Li_cpx_effective = Diffusivities()
#Li_cpx_effective.description = 'effective Li in cpx\nRichter et al. 2014'
#Li_cpx_effective.logDu = np.log10(np.array([1.8E-10, 2.6E-10, 2.5E-9, 2.1E-9]))
#Li_cpx_effective.celsius_u = [900.]*len(Li_cpx_effective.logDu)
#Li_cpx_effective.basestyle = {'color' : 'darkorchid', 'marker' : '+', 
#                                 'linestyle' : 'none', 'mew' : 3,
#                                 'markersize' : 10, 'alpha' : 1.}
#
#
#def fast_vs_slow(celsius, minutes, lengths_microns=[2000., 2000., 2000.], 
#                 sample=None, printout=False):
#    """ Takes: 
#        (1) temperature in celsius (required)
#        (2) time in minutes (required)
#        (3) sample block geometry as EITHER
#            (a) a sample with thickness information in 3 directions
#            (b) keyword lengths_microns (defaults to 2x2x2 mm)
#    Outputs 
#        (1) a plot showing fast and slow mechanism H diffusion progress
#            in olivine for the given temperature, time, and size based on 
#            Kohlstedt and Mackwell, 1998 Arrhenius law
#        (2) if printout=True, the diffusivities
#    """    
#    if sample is not None:
#        a, b, c = sample.thickness_microns[0:3]
#    else:
#        a, b, c = lengths_microns
#
#        
#    if printout is True:
#        print('slower proton-vacancy pv mechanism diffusivities')
#    D_slow = [KM98_slow.whatIsD(celsius, orient='x', printout=printout), 
#              KM98_slow.whatIsD(celsius, orient='y', printout=printout),          
#              KM98_slow.whatIsD(celsius, orient='z', printout=printout)]
#    
#    if printout is True:
#        print
#        print('faster proton-polaron pp mechanism diffusivites')
#    D_fast = [KM98_fast.whatIsD(celsius, orient='x', printout=printout), 
#              KM98_fast.whatIsD(celsius, orient='y', printout=printout),          
#              KM98_fast.whatIsD(celsius, orient='z', printout=printout)]
#    
#    fig, axes, v, x, y = diffusion3Dnpi(lengths_microns=[a, b, c],
#                                        log10Ds_m2s=D_slow, 
#                                        time_seconds=minutes*60., 
#                                        initial=1., final=0.,
#                                        plot3=True, centered=False)
#
#    v2, x2, y2 = diffusion3Dnpi(lengths_microns=[a, b, c], 
#                                styles3=[st.style_1]*3,
#                                log10Ds_m2s=D_fast, figaxis3=axes,
#                                time_seconds=minutes*60., 
#                                initial=1., final=0., plot3=True, 
#                                centered=False)
#                             
#    return fig, axes
#

#
##%%
#def print_info(Diff):
#    """ Takes Diffusivity object and prints out data nicely for copying """
#    orients = ['a', 'b', 'c', 'u']
#    for d3, temp, orient in zip(Diff.logD, Diff.celsius, orients):
#        for d, t in zip(d3, temp):
#            print(orient, t, d)
#            