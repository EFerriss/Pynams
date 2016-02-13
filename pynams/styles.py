# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:43:49 2015

@author: Ferriss
"""
import matplotlib.pyplot as plt
import numpy as np

style_points = {'color' : 'k', 'marker' : 'o', 'markersize' : 6,
                'fillstyle' : 'none', 'linestyle' : 'none',}
#                'markerfacecolor' : 'w'}
style_baseline = {'color' : 'k', 'linewidth' : 1, 'linestyle' :'-'}
style_spectrum = {'color' : 'b', 'linewidth' : 3}
style_spectrum_red = {'color' : 'r'}
style_fitpeak = {'color' : 'g', 'linewidth' : 1}
style_summed = {'color' : 'orangered', 'linewidth' : 2, 'linestyle' : '--'}
style_profile = {'markeredgecolor' : 'black', 'linestyle' : '', 'marker' : 'o', 
                 'markersize' : 10, 'markerfacecolor' : 'grey', 'alpha' : 0.5}

style_initial_point = {'color' : 'grey', 'label' : 'initial', 'linestyle' : 'none',
                       'marker' : 'o', 'alpha' : 0.5}
style_initial = {'color' : 'blue', 'label' : 'initial', 'linestyle' : '--'}                 
style_initialgrey = {'color' : 'grey', 'label' : 'initial', 'linestyle' : '-'}
style_1 = {'linestyle' : '-', 'color' : 'k', 'marker' : None, 'linewidth' : 1}
style_1a = {'linestyle' : '--', 'color' : 'k', 'marker' : None, 'linewidth' : 1}
style_2 = {'color':'green', 'linewidth': 1.5}
style_2a = {'color':'green', 'linewidth': 1.5, 'linestyle' : '--'}
style_3 = {'color':'blue', 'linewidth': 1.}
style_3a = {'color':'blue', 'linewidth': 1., 'linestyle' : '--'}
style_4 = {'color':'brown', 'linewidth':1.5, 'linestyle' : '--'}
style_4a = {'color':'brown', 'linewidth':1.5, 'linestyle' : '-.'}
style_5 = {'color':'sage', 'linewidth':2., 'linestyle' : '--'}
style_6 = {'color':'darkorchid', 'linewidth':2., 'linestyle' : '--'}

style_unoriented = {'fillstyle' : 'none'}
style_lightgreen = {'color' : 'lightgreen', 'linewidth' : 4}

style_diffusioncurve = {'color' : 'grey', 'linewidth' : 4, 
                        'linestyle' : '-', 'alpha' : 0.5}
style_errorbounds = {'color' : 'grey', 'linewidth' : 2, 
                     'linestyle' : '--', 'alpha' : 0.5}

# different profile directions
style_Dx = {'fillstyle' : 'left', 'color' : 'red', 'markerfacecolor' : 'red'}
style_Dy = {'fillstyle' : 'bottom', 'color' : 'green', 
            'markerfacecolor' : 'green' }
style_Dz = {'fillstyle' : 'right', 'color' : 'blue', 
            'markerfacecolor' : 'blue'}
style_Du = {'fillstyle' : 'none', 'color' : 'k', 
            'markerfacecolor' : 'white'}
style_Dx_line = {'linestyle' : '--', 'color' : 'red'}
style_Dy_line = {'linestyle' : '-.', 'color' : 'green'}
style_Dz_line = {'linestyle' : ':', 'color' : 'blue'}
style_Du_line = {'linestyle' : '-', 'color' : 'black'}
style_orient = [style_Dx, style_Dy, style_Dz, style_Du]            
style_orient_lines = [style_Dx_line, style_Dy_line, style_Dz_line, style_Du_line]            

# different ray paths
style_Rx = {'marker' : 'd'}
style_Ry = {'marker' : 'o'}
style_Rz = {'marker' : 's'}

style_bulk = {'marker' : 'o', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 7}

style_peak0 = {'marker' : 'x', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'r', 'markersize' : 6, 'label' : '3645 cm$^{-1}$'}

style_peak1 = {'marker' : '+', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'orange', 'markersize' : 6, 'label' : '3617 cm$^{-1}$'}

style_peak2 = {'marker' : 's', 'fillstyle' : 'full', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 6, 'markerfacecolor' : 'teal',
              'alpha' : 0.5, 'label' : '3540 cm$^{-1}$'}

style_peak3 = {'marker' : '_', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'g', 'markersize' : 6, 'mew' : 2,
              'label' : '3460 cm$^{-1}$'}

style_peak4 = {'marker' : '|', 'fillstyle' : 'none', 'linestyle' : 'none',
              'color' : 'b', 'markersize' : 6, 'mew' : 2,
              'label' : '3443 cm$^{-1}$'}

style_peak5 = {'marker' : 'd', 'fillstyle' : 'full', 'linestyle' : 'none',
              'color' : 'k', 'markersize' : 6, 'markerfacecolor' : 'violet',
              'alpha' : 0.5, 'label' : '3350 cm$^{-1}$'}

def plot_3panels_outline(style=None, top=1.2, figsize=(6.5, 2.5),
                         shrinker=0.1, heights_instead=False,
                         wholeblock=True):
    """Outline setup for 3 subplots for 3D profiles"""
    if style is None:
        style = style_lightgreen

    fig, axis3 = plt.subplots(nrows=1, ncols=3)
    fig.set_size_inches(figsize)

    for k in range(3):
        axis3[k].set_ylim(0, top)
#        axis3[k].grid()
        box = axis3[k].get_position()
        plt.setp(axis3[k].xaxis.get_majorticklabels(), rotation=45)
        axis3[k].set_position([box.x0 + box.width*shrinker, 
                               box.y0 + box.height*shrinker, 
                               box.width*(1.0-shrinker), 
                               box.height*(1.0-shrinker)])
                               
    if wholeblock is True:
        if heights_instead is False:
            axis3[0].set_ylabel('Area/Area$_0$')
        else:
            axis3[0].set_ylabel('Height/Height$_0$')
            
    else:
        if heights_instead is False:
            axis3[0].set_ylabel('Area (cm$^{-2}$)')
        else:
            axis3[0].set_ylabel('Height (cm$^{-1}$)')
        
    axis3[0].set_xlabel('|| a')
    axis3[1].set_xlabel('position ($\mu$m) || b')
    axis3[2].set_xlabel('|| c')
    plt.setp(axis3[1].get_yticklabels(), visible=False)
    plt.setp(axis3[2].get_yticklabels(), visible=False)
    return fig, axis3
    
def plot_3panels(positions_microns, area_profiles, lengths=None,
                 styles3=[None, None, None], top=1.2, figaxis3=None, 
                 show_line_at_1=True, init=1.,
                 centered=True, percent_error=3., xerror=50.,
                 heights_instead=False, wholeblock=True,
                 use_errorbar=False):
    """Make 3 subplots for 3D and 3DWB profiles. The position and area profiles
    are passed in lists of three lists for a, b, and c.
    Positions are assumed to start at 0 and then are centered.
    """
    if figaxis3 is None:
        fig, axis3 = plot_3panels_outline(top=top, wholeblock=wholeblock,
                                          heights_instead=heights_instead)
    else:
        axis3 = figaxis3

    if lengths is None:
        lengths = np.ones(3)
        for k in range(3):
            lengths[k] = max(positions_microns[k] - 
                            min(positions_microns[k]))

    for k in range(3):      
        x = positions_microns[k]
        y = area_profiles[k]
        
        if len(x) != len(y):
            print 'Problem in plot_3panels'
            print 'len(x):', len(x)
            print 'len(y):', len(y)

        a = lengths[k] / 2.
        axis3[k].set_xlim(-a, a)

        if show_line_at_1 is True:
            axis3[k].plot([-a, a], [init, init], '--k')
            
        if styles3[k] is None:
            styles3[k] = style_lightgreen
            
        if np.isnan(y).any():
            axis3[k].text(0, axis3[k].get_ylim()[1]/2., 
                         'nan values!\n\nProbably the\ninitial area was 0',
                         horizontalalignment='center', backgroundcolor='w',
                         verticalalignment='center')
                         
        elif np.isinf(y).any():
            infstring = ''.join(('inf values!\n\nProbably the',
                                '\ninitial area was 0\nand a peak grew'))
            axis3[k].text(0, axis3[k].get_ylim()[1]/2., infstring,
                         horizontalalignment='center', backgroundcolor='w',
                         verticalalignment='center')

        else:
            if use_errorbar is True:
                yerror = np.array(y) * percent_error/100.
                axis3[k].errorbar(x-a, y, xerr=xerror, yerr=yerror, 
                                **styles3[k])
            else:
                axis3[k].plot(x-a, y, **styles3[k])

    if figaxis3 is None:
        return fig, axis3   
