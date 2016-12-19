6# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:49:59 2015

@author: Ferriss

Library of diffusivities, both data and fit lines
With focus on hydrogen diffusion in minerals

First set up Diffusivities as a class with function plotD (for all
orientations) take the figure axis handle and plot the data onto it.

Also has function add_to_legend to add that information to a list of 
legend handles. 

After the classes are set up, I have a list of diffusivities
with their own marker styles set up.

"""
import numpy as np
from uncertainties import ufloat
import diffusion
import pynams
import styles as st
from diffusion import Diffusivities

reload(diffusion)

GAS_CONSTANT = 0.00831 # kJ/mol K

# plotting detail defaults
markersizefloat = 8
style_unoriented = {'fillstyle' : 'none'}
style_Dx = {'fillstyle' : 'left'}
style_Dy = {'fillstyle' : 'bottom'}
style_Dz = {'fillstyle' : 'right'}
style_Dx_line = {'linestyle' : '--'}
style_Dy_line = {'linestyle' : '-.'}
style_Dz_line = {'linestyle' : ':'}
style_unoriented_line = {'linestyle' : '-'}

# Set up class and plotting commands
def make_line(celsius_list, logD_list, lowD=6.0, highD=10.0):
    """Take list of celsius and return 1e4/T in K"""
    if len(celsius_list) != len(logD_list):
        print 'List of temperatures not equal to list of diffusivities'
        return
    Tx = []
    for k in celsius_list:
        Tx.append(1.0e4 / (k+273.15))  
    p = np.polyfit(Tx, logD_list, 1)
    x = np.linspace(lowD, highD, 100) 
    y = np.polyval(p, x)
    return x, y


#%% olivine
Wanamaker_Si = Diffusivities(
               description = "San Carlos ol. $V''''_{Si}$\nWanamaker 1994",
               celsius_unoriented = [1100., 1200.], 
               logD_unoriented = [-10.4, -9.99],
               basestyle = {'color' : 'k', 'marker' : '*', 'alpha' : 0.5,
                            'linestyle' : 'none', 'markersize' : 10})

Wanamaker_Mg = Diffusivities(
               description = "San Carlos ol. $V''_{Me}$\nWanamaker 1994",
               celsius_unoriented = [1100., 1200., 1300.],
               logD_unoriented = [-11.20, -10.8, -10.3],
               basestyle = {'color' : 'g', 'marker' : '*', 'alpha' : 0.5,
                            'fillstyle' : 'none',
                            'linestyle' : 'none', 'markersize' : 10})
                            
DM03 = Diffusivities(
       description = 'bulk H forsterite\nDemouchy & Mackwell 2003',
       celsius_x = [1100, 1058, 1056, 1015], 
       celsius_y = [1109, 1099, 1060, 1056, 1001],
       celsius_z = [1109, 1060, 1057, 999, 949, 899, 898],
       logDx = [-12.2829, -12.6956, -12.6956, -12.9896],
       logDy = [-11.7728, -11.9740, -11.9688, -12.2108, -12.5100],
       logDz = [-11.4380, -11.5100, -11.6748, -11.9740, -12.1440, -12.5924, -12.6800],
       Fe2=0., Fe3=0., Al=0., Mg=2.0, Ti=0.,
       basestyle = {'color' : 'green', 'marker' : 'v', 'alpha' : 0.5,
                    'linestyle' : 'none', 'markersize' : markersizefloat})

# Padron-Navarta et al. 2014
# [Si]: 3613 peak in MgO-buffered Fo, [Si]-[Ti]: 3613 peak in Ti-doped sample
# [Ti]: 3525 peak, [Mg]: 3220 peak
pnav_Si = Diffusivities(description = 'forsterite, [Si]', 
          celsius_unoriented = [1000, 1100, 1200],
          logD_unoriented = [-15.66, -14.17, -13.02],
          basestyle={'color':'cornflowerblue', 'marker':'o', 'fillstyle':'full',
                     'markersize':10})

pnav_SiTi = Diffusivities(description = 'forsterite, [Si-Ti]',
            celsius_unoriented = [1000, 900, 800],
            logD_unoriented = [-13.35, -14.3, -15.29],)

pnav_Ti = Diffusivities(description = 'forsterite, [Ti]',
          celsius_unoriented = [1000, 900, 800],
          logD_unoriented = [-13.16, -14.1, -15.26],
          basestyle={'color':'orange', 'marker':'s', 'fillstyle':'full',
                     'markersize':10})
                        
pnav_Mg = Diffusivities(description = 'forsterite, [Mg]',
          celsius_unoriented = [1000, 900, 800],
          logD_unoriented = [-12.66, -13.7, -14.29],)

# Du Frane & Tybursky 2012 self-diffusion in olivine
DuFrane = Diffusivities(
          description = 'H self-diffusion\nDuFrane & Tybursky 2012',
          celsius_x = [750, 800, 900], logDx = [-12.3, -11.7, -11.2],
          celsius_z = [900], logDz = [-12], logDx_error = [0.3, 0.1, 0.1],
          logDz_error = [0.2],
          basestyle = {'color' : 'purple', 'marker' : 's', 'alpha' : 0.5,
                       'markersize' :  markersizefloat})

# Kohlstedt and Mackwell 1hr data "fast proton-polaron mechanism"
# Data from Frey and Prinz 1978
KM98_fast = Diffusivities(description = 'fast mech., KM98',
            celsius_x = [1000, 1000, 900, 900, 800],
            celsius_y = [1000, 1000, 900, 900, 800],
            celsius_z = [1000, 1000, 900, 900, 800],
            logDx = [-9.5634, -9.3928, -9.9876, -10.3970, -10.8071],
            logDy = [-11.0066, -10.6898, -11.3043, -11.9967, -12.2306],
            logDz = [-11.4068, -11.0066, -11.2551, -11.1582, -12.1140],
            Fe3=0.003, Fe2=0.174, Mg=1.821, Ti=0., Al=0.,
            basestyle = {'marker' : 'd', 'color' : 'grey', 'alpha' : 0.5,
                        'markersize' :  10})

# Demouchy&Mackwell 2006 - 1 hr for all "fast"
DM06_fast = Diffusivities(description = 'fast mech., DM06', 
            celsius_all = [900], logDx = [np.log10(4e-11)], 
            logDy = [np.log10(2e-12)], logDz = [np.log10(1e-12)],
            basestyle = {'marker' : 's', 'color' : 'navy', 'alpha' : 0.5,
                        'markersize' :  markersizefloat})

# Kohlstedt&Mackwell 1998, 8 hour data "slow"
KM98_slow = Diffusivities(description = 'slow mech., KM98', 
            celsius_all = [900, 1000], logDx = [-13.994, -12.783],
            logDy = [-14.17, -13.171], logDz = [-12.885, -11.776],
            basestyle = {'marker' : 'D', 'color' : 'green', 'alpha' : 0.5,
                      'markersize' :  10})

#% Demouchy & Mackwell 2006 - 20 hour and 5 hour - "slow"
DM06_slow = Diffusivities(description = 'slow mech., DM06', celsius_x = [900],
            celsius_y  = [900, 1000], celsius_z = [900], 
            logDx = [np.log10(5e-13)], 
            logDy = [np.log10(5e-14), np.log10(3e-13)], 
            logDz = [np.log10(1e-12)],
            basestyle = {'marker' : 's', 'color' : KM98_slow.basestyle['color'], 
                         'alpha' : 0.5, 'markersize' : markersizefloat+3})

# Ferriss et al. 2015- single mechanism in San Carlos olivine
Ferriss_single = Diffusivities(description = 'San Carlos ol.\nsingle mech.',
                 celsius_all = [900], logDx = [-11.3], logDy = [-13.6],
                 logDz = [-12.1], logDx_error = [0.1], logDy_error = [0.6],
                 logDz_error = [0.2], 
                 basestyle = {'color' : 'orangered', 'marker' : '*', 
                              'markersize' :  markersizefloat+5})

#%  Ferriss et al. 2015 - slow mechanism in simultaneous model 
Ferriss_slow = Diffusivities(
               description = 'Slow mechanism,\nsimultaneous model',
               celsius_all = [900], logDx = [-12.2], logDy = [-13.2], 
               logDz = [-12.3], logDx_error = [0.3], logDy_error = [0.5],
               logDz_error = [0.3],
               basestyle = {'marker' : 'D', 'color' : 'yellow', 'alpha' : 1,
                            'markersize' : markersizefloat})

# ---------------- natural olivine data -----------------------------
# See also OlivineNormalization.xls

# Loihi, Hawaii: Hauri 2002 Part 2
# Chemistry estimated from Garcia 1995 Table 2; Al and Ti not measured
Hauri02 = Diffusivities(description = 'Hauri 2002', 
                        celsius_unoriented = [1275], 
                        logD_unoriented = [np.log10(4e-9)], 
                        logDu_error = [], 
                        Fe3=0.004, Fe2=0.261, Mg=1.727, Al=0., Ti=0.,
                        basestyle = {'marker' : '+', 'color' : 'darkorange', 
                                     'markersize' : 12,
                                     'markeredgewidth' : '2'})

# Portnyagin 2008 minimum value, see new Miranov paper
# See Supplementary Materials. Al and Ti not measured in host olivines
Portnyagin = Diffusivities(
             description = 'Portnyagin et al. 2008\nMironov et al. 2015',
             celsius_unoriented = [1140.], logD_unoriented = [np.log10(5e-12)],
             Fe3=0.021, Fe2=0.276, Mg=1.7, Al=0., Ti=0.,
             basestyle = {'marker' : 's', 'color' : 'darkmagenta', 
                          'markerfacecolor' : 'white',
                          'markeredgecolor' : 'darkmagenta',
                          'markersize' : markersizefloat,
                          'markeredgewidth' : '1'})

# Chen et al. 2011
# Chemistry from Le Voyer 2014 Table 1 SOM1 olivine; Al and Ti not measured...
Chen11 = Diffusivities(description = 'Chen et al. 2011',
                       celsius_unoriented = [1533.-273.15, 1471.-273.15, 
                                             1437.-273.15, 1561.-273.15],
                       logD_unoriented  = [np.log10(2e-11), np.log10(2.5e-11), 
                                           np.log10(0.5e-11), np.log10(2.5e-11)],
                       Fe3=0.0, Fe2=0.192, Mg=1.784, Al=0., Ti=0.,
                       basestyle = {'marker' : 'x', 'color' : 'black',
                                     'markersize' : markersizefloat,
                                     'markeredgewidth' : '1'})
                    
# Gaetani et al. 2012 hydration data
# Moana Loa olivine with melt inclusion unoriented grain; 
# MID-1 data in Supplement, Ti and Al not measured...
# *** Si deficiency ***
Gaetani12 = Diffusivities(description = 'Gaetani et al., 2012',
                          celsius_unoriented = [1250.],
                          logD_unoriented = [np.log10(1.7e-11)],
                          Fe3=0.036, Fe2=0.206, Mg=1.769, Al=0., Ti=0.,
                          basestyle = {'marker' : 'h', 'color' : 'red',
                                       'markeredgecolor' : 'red', 
                                       'markersize' : 10,
                                       'markeredgewidth' : '1',
                                       'markerfacecolor' : 'white'})


generic = Diffusivities()
generic.basestyle = {'marker' : 's', 'color' : 'black', 'alpha' : 0.5,
                     'markersize' : markersizefloat, 'linestyle': 'none'}

#%% clinopyroxene
H_CrDiopside_Ingrin95 = Diffusivities(
                        description = ''.join(('H in Cr-rich diopside in air',
                               '\nmainly 3645 cm$^{-1}$\nIngrin et al. 1995')),
                        celsius_all = np.array([973, 1073, 1173, 1273]) - 273.,
                        logDx = [np.log10(3E-14), np.log10(9E-14),
                                 np.log10(5.5E-13), np.log10(1.9E-12)],
                        logDz = [np.log10(5E-14), np.log10(13E-14), 
                                np.log10(4.5E-13), np.log10(12E-13)],
                        basestyle = {'color' : 'green', 'marker' : '^', 
                                   'markersize' :  12, 
                                   'linestyle' : 'none', 'alpha' : 1.})

H_cpx_Wade2008 = Diffusivities(
                 description = 'H in cpx phenocryst\nWade et al. 2008',
                 celsius_all = [1100], logD_unoriented = [-13.],
                 basestyle = {'color' : 'indigo', 'marker' : '$\clubsuit$', 
                            'markersize' :  12,
                            'markeredgewidth' : 0,
                            'linestyle' : 'none', 'alpha' : 1.})

H_diopside_Woods00 = diffusion.Diffusivities(
                     description = 'Jaipur bulk H\nWoods et al. 2000',
                     logDx = np.log10(np.array([5e-13, 7e-12, 1.1e-11, 1.5e-11, 
                                          2e-11, 1.8e-11, 3.5e-12, 3.5e-11])),
                     logDy = np.log10(np.array([1.5e-12, 3e-12, 2.5e-12, 8e-13, 
                                                2e-12])),
                     logDz = np.log10(np.array([2e-12, 7e-12, 1.5e-11, 1.5e-11, 
                                                6e-12, 3e-11])),
                     celsius_x = [700, 750, 800, 800, 800, 850, 700, 850],
                     celsius_y = [750, 800, 800, 800, 850],
                     celsius_z = [700, 750, 800, 800, 800, 850],)
#                     basestyle = {'color' : 'k', 'marker' : 's', 
#                                  'markerfacecolor' : 'turquoise',
#                                  'markersize' :  10, 
#                                  'linestyle' : 'none', 'alpha' : 1.,
#                                  'markeredgewidth' : 0.5})
                            
H_cpxBasanite_Xia00 = Diffusivities(
                      description = ''.join(('H in basanite cpx\n3630 + 3500',
                            ' cm$^{-1}$\nf$_{02}=10^{-14}$ Xia et al. 2000')),
                      logDz = [np.log10(1.8E-12), np.log10(6.5E-12)],
                      celsius_z = [850, 950],
                      basestyle = {'color' : 'k', 'marker' : '>', 
                                   'markerfacecolor' : 'k',
                                   'markersize' :  8, 
                                   'linestyle' : 'none', 'alpha' : 1.,
                                   'markeredgewidth' : 1})

H_diopside_noFe = Diffusivities(
                  description = ''.join(('pure synth. diopside in air',
                                         '\nSundvall et al. 2009')),
                  logDx = [-12.3, -14.7], 
                  celsius_x = np.array([1273, 1073]) - 273.,
                  logDy = [-12.6, -13.9, -15.1],
                  celsius_y = np.array([1273, 1173, 1073]) - 273.,
                  basestyle = {'color' : 'b', 'marker' : '>', 
                               'markerfacecolor' : 'b', 'markersize' :  8, 
                               'linestyle' : 'none', 'alpha' : 1.,
                               'markeredgewidth' : 1})

H_diopside_Sundvall = Diffusivities(
                      description = ''.join(('synth. Fe-bearing diopside\n',
                                             'in air; Sundvall et al. 2009')),
                      logDy = [-13.7, -15.3, -15.9],
                      celsius_y = np.array([1000., 900., 800.]),)
#                      basestyle = {'color' : 'y', 'marker' : 'v', 
#                                   'markerfacecolor' : 'y', 'markersize' :  13, 
#                                   'linestyle' : 'none', 'alpha' : 1.,
#                                   'markeredgewidth' : 1})

#Fe_diopside = Diffusivities()
#Fe_diopside.description = 'Fe in diopside\nAzough & Freer 2000'
#Fe_diopside.celsius_unoriented = np.array([950., 1100.])
#Fe_diopside.activation_energy_kJmol = 161.5
#Fe_diopside.D0_m2s = 6.22E-15
#Fe_diopside.logD_unoriented = Fe_diopside.makelines_Ea_D0(Ea=161.5, log10D0=10.**(6.22E-15), 
#                                  celsius_list=Fe_diopside.celsius_unoriented)
#Fe_diopside.basestyle = {'color' : 'k', 'marker' : 'None', 
#                         'linestyle' : '-', 'alpha' : 1.}


# H in Kunlun diopside
# 1000 C, in whole-block paper, label: K5
Dx_wholeblock = ufloat(-13.2, 0.2)
Dy_wholeblock = ufloat(-13.4, 0.2)
Dz_wholeblock = ufloat(-13.6, 0.3)
Dy_slice_FTIR = ufloat(-13.1, 0.3)
Dz_slice_FTIR = ufloat(-13.1, 0.2)
Dy_slice_SIMS = ufloat(-13.3, 0.4)
Dz_slice_SIMS = ufloat(-13.2, 0.4)
Dy_ave = np.mean([Dy_slice_FTIR, Dy_slice_SIMS, Dy_wholeblock])
Dz_ave = np.mean([Dz_slice_FTIR, Dz_slice_SIMS, Dz_wholeblock])


Kunlun_bulkH = Diffusivities()
Kunlun_bulkH.description = 'Kunlun diopside\nbulk H, QFM'
Kunlun_bulkH.celsius_all = [1000.]
Kunlun_bulkH.logDx = [-13.0]
Kunlun_bulkH.logDy = [-13.5]
Kunlun_bulkH.logDz = [-13.5]
Kunlun_bulkH.basestyle = {'color' : 'black', 'marker' : 'D', 
                          'markersize' :  markersizefloat, 
                          'linestyle' : 'none', 'alpha' : 0.5,}
Kunlun_bulkH.logDx_error = [0.1]
Kunlun_bulkH.logDy_error = [0.1]
Kunlun_bulkH.logDz_error = [0.1]

Kunlun_peak3617 = Diffusivities()
Kunlun_peak3617.description = 'Kunlun diopside\nPeak at 3617 cm$^{-1}$'
Kunlun_peak3617.celsius_all = [1000.]
Kunlun_peak3617.logDx = [-13.4]
Kunlun_peak3617.logDy = [-12.3]
Kunlun_peak3617.logDz = [-13.5]
Kunlun_peak3617.logDx_error = [0.4]
Kunlun_peak3617.logDy_error = [0.03]
Kunlun_peak3617.logDz_error = [0.2]
Kunlun_peak3617.basestyle = dict(Kunlun_bulkH.basestyle.items())
Kunlun_peak3617.basestyle['color'] = 'red'

Kunlun_peak3540 = Diffusivities()
Kunlun_peak3540.description = 'Kunlun diopside\nPeak at 3540 cm$^{-1}$'
Kunlun_peak3540.celsius_all = [1000.]
Kunlun_peak3540.logDx = [-12.7]
Kunlun_peak3540.logDy = [-12.2]
Kunlun_peak3540.logDz = [-12.7]
Kunlun_peak3540.logDx_error = [0.3]
Kunlun_peak3540.logDy_error = [0.1]
Kunlun_peak3540.logDz_error = [0.3]
Kunlun_peak3540.basestyle = dict(Kunlun_bulkH.basestyle.items())
Kunlun_peak3540.basestyle['color'] = 'blue'

# Jaipur diopside, Woods et al. 2001
Jaipur_bulk = Diffusivities()
Jaipur_bulk.description = 'Jaipur diopside bulk H\n(Woods et al. 2001)'
Jaipur_bulk.celsius_x = [700, 750, 800, 800, 800, 850, 700, 850]
Jaipur_bulk.celsius_y = [750, 800, 800, 800, 850]
Jaipur_bulk.celsius_z = [700, 750, 800, 800, 800, 850]
Jaipur_bulk.logDx = np.log10([5e-13, 7e-12, 1.1e-11, 1.5e-11, 2e-11, 1.8e-11, 
                     3.5e-12, 3.5e-11])
Jaipur_bulk.logDy = np.log10([1.5e-12, 3e-12, 2.5e-12, 8e-13, 2e-12])
Jaipur_bulk.logDz = np.log10([2e-12, 7e-12, 1.5e-11, 1.5e-11, 6e-12, 3e-11])
Jaipur_bulk.basestyle = {'color' : 'green', 'marker' : 's',
                         'markersize' : markersizefloat, 
                         'linestyle' : 'none'}

O_diopside_self = Diffusivities()
O_diopside_self.description = 'O self-diffusion in diopside\nRyerson & McKeegan 1994'
O_diopside_self.logDz = [np.log10(1.96e-21), np.log10(2.04e-21), np.log10(1.48e-21),
                         np.log10(2.29e-21), np.log10(1.19e-20), np.log10(2.47e-20),
                         np.log10(3.63e-20), np.log10(1.9e-20), np.log10(3.17e-20),
                         np.log10(1.75e-20), np.log10(2.17e-20), np.log10(3.74e-20),
                         np.log10(8.94e-20), np.log10(9.64e-20), np.log10(1.22e-19)]
O_diopside_self.celsius_z = [1104, 1104, 1104, 1105, 1150, 1200, 1200, 1200, 
                             1200, 1202, 1202, 1202, 1251, 1251, 1251]
O_diopside_self.basestyle = {'color' : 'red', 'marker' : '^', 
                           'linestyle' : 'none',
                           'fillstyle' : 'right', 'alpha' : 0.5}                             

U_diopside = Diffusivities()
U_diopside.description = 'U in diopside\nVan Orman et al. 1998'
U_diopside.celsius_z = [1150, 1150, 1200, 1200, 1200, 1300, 1300]  
U_diopside.logDz = [np.log10(4.54e-22), np.log10(9.44e-22), np.log10(2.90e-21),
                    np.log10(3.95e-21), np.log10(2.35e-21), np.log10(2.35e-21),
                    np.log10(2.18e-20)]
U_diopside.basestyle = {'color' : 'orange', 'marker' : 'o', 
                           'linestyle' : 'none', 'markersize' : 7, 
                           'fillstyle' : 'right', 'alpha' : 0.5}                             

Th_diopside = Diffusivities()
Th_diopside.description = 'Th in diopside\nVan Orman et al. 1998'
Th_diopside.celsius_z = [1150, 1150, 1200, 1200, 1200, 1300, 1300]  
Th_diopside.logDz = [np.log10(1.59e-21), np.log10(1.2e-21), np.log10(5.26e-21),
                     np.log10(3.80e-21), np.log10(3.22e-21), np.log10(2.12e-20),
                     np.log10(2.75e-20)]
Th_diopside.basestyle = {'color' : 'yellow', 'marker' : 'o', 
                           'linestyle' : 'none', 'markersize' : 7, 
                           'fillstyle' : 'right', 'alpha' : 0.5}                             

Ce_diopside = Diffusivities()
Ce_diopside.description = 'Ce in diopside\nVan Orman et al. 2001'
Ce_diopside.celsius_z = [1300, 1275, 1250, 1250, 1225, 1200, 1200, 1200, 1200, 
                         1175, 1150]
Ce_diopside.logDz = [np.log10(31.9e-21), np.log10(25.0e-21), np.log10(11.5e-21),
                     np.log10(6.83e-21), np.log10(5.83e-21), np.log10(2.53e-21),
                     np.log10(4.45e-21), np.log10(4.01e-21), np.log10(4.53e-21),
                     np.log10(0.68e-21), np.log10(0.62e-21)]
Ce_diopside.basestyle = {'color' : 'lime', 'marker' : 'o', 
                           'linestyle' : 'none', 'markersize' : 7, 
                           'fillstyle' : 'right', 'alpha' : 0.5}                             

Al_diopside = Diffusivities()
Al_diopside.description = 'Al in diopside; Sautter\net al. EPSL 1988'
Al_diopside.logD_unoriented = [-18.495]
#Al_diopside.logDu_error = [0.66]
Al_diopside.celsius_unoriented = [1180.]
Al_diopside.basestyle = {'color' : 'g', 'marker' : 'x',
                         'fillstyle' : 'full'}

CaMg_diopside_2010 = Diffusivities()
CaMg_diopside_2010.description = 'Ca-Mg interdiffusion in di.\nZhang et al. 2010'
CaMg_diopside_2010.logDx = [-19.33, -19.37, -19.92, -19.55, -19.99, -19.97, 
                            -20.39, -20.41, -20.82, -21.18]
CaMg_diopside_2010.celsius_x = [1150.00, 1150.00, 1100.00, 1100.00, 1050.00, 
                                1050.00, 1000.00, 1000.00, 950.00, 950.00]
CaMg_diopside_2010.logDy = [-19.39, -19.38, -20.06, -20.07, -20.95, -21.08, 
                            -20.92, -21.18, -21.26, -21.36]
CaMg_diopside_2010.celsius_y = [1150.00, 1150.00, 1050.00, 1050.00, 1000.00, 
                                1000.00, 1000.00, 1000.00, 950.00, 950.00]
CaMg_diopside_2010.logDz = [-19.49, -19.44, -20.12, -20.24, -21.36, -21.38]
CaMg_diopside_2010.celsius_z = [1150.00, 1150.00, 1050.00, 1050.00, 950.00, 950.00]
                            
CaMg_diopside_2010.basestyle = {'color' : 'orange', 'marker' : 's', 
                                'linestyle' : 'none', 'markersize' : 7,
                                'fillstyle' : 'full', 'alpha' : 0.8}

CaMg_diopside_1983 = Diffusivities()
CaMg_diopside_1983.description = 'Ca-Mg interdiffusion in cpx\nBrady & McCallister 1983'
CaMg_diopside_1983.logD_unoriented = np.array([np.log10(1.4E-16), np.log10(6.9E-17), np.log10(5.6E-17),
                                 np.log10(2.0E-17), np.log10(5.7E-18), np.log10(5.6E-16),
                                 np.log10(2.8E-16), np.log10(1.7E-16), np.log10(6.9E-17),
                                 np.log10(8.3E-16), np.log10(4.2E-16), np.log10(5.6E-16)]) - 2.
CaMg_diopside_1983.celsius_all = [1100., 1100., 1100., 1100., 1100., 1150., 1150., 
                             1150., 1150., 1200., 1200., 1250]
CaMg_diopside_1983.basestyle = {'color' : 'orange', 'marker' : 's', 
                           'linestyle' : 'none', 'markersize' : 9,
                           'fillstyle' : 'none', 'alpha' : 0.5}
                           
FeMg_cpx_2013 = Diffusivities()
FeMg_cpx_2013.description = 'Fe-Mg interdiffusion in cpx\nMueller et al. 2013'
FeMg_cpx_2013.logDz = [-21, -20.8, -20.46, -18, -21.89, -21.05, -18.64, -19.75, 
                       -18.3, -17.52, -19.52, -20.4, -20.46, -19.41, -18.7, 
                       -20.4, -19.7, -18.7, -19.92, -19.52, -19.82, -19.92]
FeMg_cpx_2013.celsius_z = [850, 900, 950, 1150, 800, 905, 1106, 1006, 1154, 
                           1200, 1035, 924, 956, 1048, 1102, 945, 999, 1100, 
                           1000, 1007, 1007, 1007]
FeMg_cpx_2013.basestyle = {'color' : 'grey', 'marker' : 'h', 
                           'linestyle' : 'none', 'markersize' : 7,
                           'fillstyle' : 'right', 'alpha' : 0.8}
                         
FeMg_diopside = Diffusivities()
FeMg_diopside.description = 'Fe-Mg interdiffusion\nDimanov & Wiedenbeck 2006'
FeMg_diopside.logDz = np.array([-16.6649, -14.9891, -14.5297, -14.6098, -15.0794, 
                       -15.6405, -15.9297, -15.8744, -16.139, -15.224, 
                       -16.0353, -14.1226]) - 2. # original reported in log D (cm2/s)
FeMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
                           1050, 1150, 1100, 1100]
FeMg_diopside.basestyle = {'color' : 'green', 'marker' : 'o', 
                           'linestyle' : 'none',
                           'fillstyle' : 'right', 'alpha' : 0.5}

MnMg_diopside = Diffusivities()
MnMg_diopside.description = 'Mn-Mg interdiffusion\nDimanov & Wiedenbeck 2006'
MnMg_diopside.logDz = np.array([-16.6923, -15.0946, -14.8291, -14.8714, -15.223, 
                       -15.711, -15.9744, -16.0584, -16.2665, -15.8388, 
                       -16.2668, -14.3497]) - 2. # original reported in log D (cm2/s)
MnMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
                           1050, 1150, 1100, 1100]
MnMg_diopside.basestyle = {'color' : 'yellow', 'marker' : 'o', 
                           'linestyle' : 'none',
                           'fillstyle' : 'right', 'alpha' : 0.5}

FeMnMg_diopside = Diffusivities()
FeMnMg_diopside.description = '(Fe,Mn)Mg interdiffusion in di.\nDimanov & Wiedenbeck 2006'
FeMnMg_diopside.logDz = np.array([-16.6784, -15.0387, -14.6541, -14.7212, -15.1453, 
                         -15.6743, -15.9515, -15.9567, -16.1981, -15.4307, 
                         -16.1358, -14.2215]) - 2. # original reported in log D (cm2/s)
FeMnMg_diopside.celsius_z = [1000, 1150, 1100, 1190, 1100, 1100, 1100, 1100, 
                           1050, 1150, 1100, 1100]
FeMnMg_diopside.basestyle = {'color' : 'g', 'marker' : 's', 
                           'linestyle' : 'none',
                           'fillstyle' : 'right', 'alpha' : 0.5}

Ti_diopside = Diffusivities()
Ti_diopside.description = 'Ti in diopside\nCherniak & Liang 2012'
Ti_diopside.logDz = [-19.89, -19.90, -20.25, -20.43, -20.98, -20.92, -21.41, 
                     -21.05, -21.43, -21.58, -22.22, -22.27, -22.24, -22.77]
Ti_diopside.celsius_z = [1250, 1200, 1200, 1151, 1102, 1102, 1090, 1052, 999,
                         1000, 952, 954, 905, 905]
Ti_diopside.basestyle = {'color' : 'b', 'marker' : 'o', 
                           'linestyle' : 'none', 'markersize' : 7, 
                           'fillstyle' : 'right', 'alpha' : 0.5}
He_cpx = Diffusivities()
He_cpx.description = 'He in cpx\nTrull & Kurz 1993'
He_cpx.logD_unoriented = np.array([np.log10(1.26E-10), np.log10(1.34E-10), 
                                   np.log10(1.32E-10), np.log10(1.17E-10),
                                   np.log10(7.24E-10), np.log10(9.29E-10), 
                                   np.log10(8.04E-10), np.log10(6.52E-10),
                                   np.log10(1.47E-8), np.log10(1.48E-8), 
                                   np.log10(9.19E-9), np.log10(5.06E-9), 
                                   np.log10(4.12E-9)]) - 2.
He_cpx.celsius_unoriented = [965., 965., 965., 965., 1070., 1070., 1070., 1070.,
                             1170., 1170., 1170., 1170., 1170.]
He_cpx.basestyle = {'color' : 'lawngreen', 'marker' : '$\spadesuit$', 
                    'linestyle' : 'none', 
                    'markersize' : 10, 'fillstyle' : 'full', 'alpha' : 1.}

Li_cpx_interstitial = Diffusivities()
Li_cpx_interstitial.description = 'Li in cpx as interstitial\nRichter et al. 2014'
Li_cpx_interstitial.logD_unoriented = np.log10(np.array([9.2E-10, 1.6E-9, 1.6E-9, 
                                     3.3E-8, 3.7E-8, 1.2E-8, 4.6E-10, 6.4E-12]))
Li_cpx_interstitial.celsius_unoriented = [900.]*len(Li_cpx_interstitial.logD_unoriented)
Li_cpx_interstitial.basestyle = {'color' : 'darkorchid', 'marker' : '+', 
                                 'linestyle' : 'none', 
                                 'markersize' : 10, 'alpha' : 1.}
                             
Li_cpx_effective = Diffusivities()
Li_cpx_effective.description = 'effective Li in cpx\nRichter et al. 2014'
Li_cpx_effective.logD_unoriented = np.log10(np.array([1.8E-10, 2.6E-10, 2.5E-9, 2.1E-9]))
Li_cpx_effective.celsius_unoriented = [900.]*len(Li_cpx_effective.logD_unoriented)
Li_cpx_effective.basestyle = {'color' : 'darkorchid', 'marker' : '+', 
                                 'linestyle' : 'none', 'mew' : 3,
                                 'markersize' : 10, 'alpha' : 1.}


def fast_vs_slow(celsius, minutes, abc=[2000., 2000., 2000.], sample=None, 
                 printout=False):
    """Required input: (1) temperature in celsius, (2) time in minutes,
    (3) either sample, an instance of pynams.Samples (gets taken first) or
    dimensions in microns using keyword abc (defaults to 2x2x2 mm)
    Output: two lists of diffusivities: Kohlstedt and Mackwell 1998 
    fast and slow mechanisms for H diffusion in San Carlos olivine"""    
    if sample is not None:
        a, b, c = pynams.get_3thick(sample)
    else:
        a, b, c = abc

    if printout is True:
        print 'slower proton-vacancy pv mechanism diffusivities'
    D_slow = [KM98_slow.whatIsD(celsius, orient='x', printout=printout), 
              KM98_slow.whatIsD(celsius, orient='y', printout=printout),          
              KM98_slow.whatIsD(celsius, orient='z', printout=printout)]
    

    if printout is True:
        print
        print 'faster proton-polaron pp mechanism diffusivites'
    D_fast = [KM98_fast.whatIsD(celsius, orient='x', printout=printout), 
              KM98_fast.whatIsD(celsius, orient='y', printout=printout),          
              KM98_fast.whatIsD(celsius, orient='z', printout=printout)]
    
    fig, axes, v, x, y = diffusion.diffusion3Dnpi(lengths_microns=[a, b, c],
                                                  log10Ds_m2s=D_slow, 
                                                  time_seconds=minutes*60., 
                                                  initial=1., final=0.,
                                                  plot3=True, centered=False)

    v2, x2, y2 = diffusion.diffusion3Dnpi(lengths_microns=[a, b, c], 
                                          styles3=[st.style_1]*3,
                                          log10Ds_m2s=D_fast, figaxis3=axes,
                                          time_seconds=minutes*60., 
                                          initial=1., final=0., plot3=True, 
                                          centered=False)
                             
    return fig, axes
