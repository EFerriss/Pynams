# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:43:49 2015

@author: Ferriss

Contains my most commonly used plotting style dictionaries (e.g., blue dots)
and some frequently used plotting setups, e.g., 3 subplots

"""
from __future__ import print_function, division, absolute_import
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import numpy as np

style_points = {'color' : 'b', 'marker' : 'o', 'markersize' : 6,
                'fillstyle' : 'full', 'linestyle' : 'none',}
style_lightgreen = {'color' : 'lightgreen', 'linewidth' : 4}
style_blue = {'color' : 'blue', 'linestyle' : '-', 'marker': None}
style_points0 = {'color' : 'black', 'marker' : 'D', 'markersize' : 6,
                 'fillstyle' : 'full', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial',
                 'markerfacecolor' : 'k'}
style_points1 = {'color' : 'red', 'marker' : '^', 'markersize' : 6,
                 'fillstyle' : 'full', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial',
                 'markerfacecolor' : 'red'}
style_points2 = {'color' : 'indigo', 'marker' : 'o', 'markersize' : 6,
                 'fillstyle' : 'none', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial',
                 'markerfacecolor' : 'w'}
style_points3 = {'color' : 'blue', 'marker' : 's', 'markersize' : 4,
                 'fillstyle' : 'none', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial'}
style_points4 = {'color' : 'green', 'marker' : 'd', 'markersize' : 5,
                 'fillstyle' : 'none', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial'}
style_points5 = {'color' : 'darkgoldenrod', 'marker' : 'D', 'markersize' : 4,
                 'fillstyle' : 'full', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial',
                 'markerfacecolor' : 'y'}
style_points6 = {'color' : 'orangered', 'marker' : 'o', 'markersize' : 5,
                 'fillstyle' : 'none', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial'}
style_points7 = {'color' : 'violet', 'marker' : 'v', 'markersize' : 7,
                 'fillstyle' : 'full', 'linestyle' : 'none', 
                 'linewidth' : 1, 'alpha' : 1, 'label' : 'initial'}
style_baseline = {'color' : 'k', 'linewidth' : 1, 'linestyle' :'-'}
style_spectrum = {'color' : 'b', 'linewidth' : 3}
style_spectrum_red = {'color' : 'r'}
style_fitpeak = {'color' : 'g', 'linewidth' : 1}
style_summed = {'color' : 'orangered', 'linewidth' : 2, 'linestyle' : '--'}
style_profile = {'markeredgecolor' : 'black', 'linestyle' : '', 'marker' : 'o', 
                 'markersize' : 10, 'markerfacecolor' : 'grey', 'alpha' : 0.5}
style_initial = {'color' : 'blue', 'label' : 'initial', 'linestyle' : '--'}                 
style_1a = {'linestyle' : '--', 'color' : 'k', 'marker' : None, 'linewidth' : 1}
style_1 = {'linestyle' : '-', 'color' : 'k', 'marker' : None, 'linewidth' : 1}
style_2 = {'color':'red', 'linewidth': 2.5, 'linestyle' : '-.'}
style_2a = {'color':'green', 'linewidth': 1.5, 'linestyle' : '--'}
style_3 = {'color':'orange', 'linewidth': 2., 'linestyle' : '--'}
style_3a = {'color':'blue', 'linewidth': 1., 'linestyle' : '--'}
style_4 = {'color':'yellow', 'linewidth':2.5, 'linestyle' : '-'}
style_4a = {'color':'brown', 'linewidth':1.5, 'linestyle' : '-.'}
style_5 = {'color':'green', 'linewidth':2.5, 'linestyle' : '-.'}
style_6 = {'color':'cyan', 'linewidth':2., 'linestyle' : '--'}
style_7 = {'color':'steelblue', 'linewidth':3., 'linestyle' : '-.'}
style_8 = {'color':'violet', 'linewidth':2., 'linestyle' : '-'}
style_grey = {'color' : 'grey', 'linewidth':4, 'linestyle' : '-'}

## different profile directions
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


def get_iorient(orient):
    """Converts x, y, z, u and a, b, c, u to 0, 1, 2, 3. 
    This is a helper function for bound methods in class diffusivitiy 
    and for determining thickness from raypath for wholeblocks"""
    if orient == 'x' or orient == 'a':
        iorient = 0
    elif orient == 'y' or orient == 'b':
        iorient = 1
    elif orient == 'z' or orient == 'c':
        iorient = 2
    elif orient == 'u' or orient == None:
        iorient = 3
    else:
        iorient = orient
    return iorient

        
def ylim_picker(spectrum, wn_xlim_left=4000, wn_xlim_right=3000, pad_top=0.1, 
                pad_bot=0., raw_data=False):
    """
    Takes a Spectrum object and returns reasonable min and max values for 
    y-axis of plots based on the absorbance values for the specified wavenumber
    range and padded top and bottom with pad variable
    """
    try:
        if spectrum.thickness_microns is None:
            absorbance = spectrum.abs_raw
        else:
            spectrum.start_at_zero(wn_xlim_left=wn_xlim_left,
                               wn_xlim_right=wn_xlim_right)
            absorbance = spectrum.abs_full_cm
    except AttributeError:
        absorbance = spectrum.abs_raw
        
    idx_lo = (np.abs(spectrum.wn_full-wn_xlim_right)).argmin()
    idx_hi = (np.abs(spectrum.wn_full-wn_xlim_left)).argmin()
    
    y = absorbance[idx_lo:idx_hi]

    bottom = min(y) 
    top = max(y)
    ylow = bottom - pad_bot
    yhigh = top + pad_top

    return ylow, yhigh
          

def make_line_style(direction, style_marker):
    """Take direction and marker style and return line style dictionary
    that reflects the direction (x, y, z, or u for unoriented) with the 
    color of the base style"""
    if direction == 'x':
        d = style_Dx_line
    if direction == 'y':
        d = style_Dy_line
    if direction == 'z':
        d = style_Dz_line
    if direction == 'u':
        d = style_Du_line
    d.update({'linewidth' : 2})
    return d


def plot_spectrum_outline(size_inches=(6, 6), shrinker=0.15,
                          figaxis=None, wn_xlim_left=4000., 
                          wn_xlim_right=3000., pad_top=0.1, 
                          pad_bot=0., raw_data=False):
    """
    Makes a standard figure outline for plotting FTIR spectra.
    Returns the figure and axis handles.
    """
    if figaxis is None:
        f, ax = plt.subplots(figsize=size_inches)
    else:
        ax = figaxis
    ax.set_xlabel('Wavenumber (cm$^{-1})$')
    ax.set_ylabel('Absorbance (cm$^{-1})$')
    ax.set_xlim(wn_xlim_left, wn_xlim_right)
    ax.grid()    
    box = ax.get_position()
    ax.set_position([box.x0 + box.width*shrinker, 
                     box.y0 + box.height*shrinker, 
                     box.width*(1.0-shrinker), 
                     box.height*(1.0-shrinker)])

    plt.setp(ax.get_xticklabels(), rotation=45)
    
    if figaxis is None:
        return f, ax
    

def plot_area_profile_outline(centered=True, peakwn=None,
                              set_size=(6.5, 4), ytop=1.2, 
                              wholeblock=False, heights_instead=False,
                              show_water_ppm=True):
    """
    Set up area profile outline and style defaults. 
    Default is for 0 to be the middle of the profile (centered=True).
    """
    fig = plt.figure(figsize=set_size)
    ax = SubplotHost(fig, 1,1,1)
    fig.add_subplot(ax)

    ax_ppm = ax.twinx()
    ax_ppm.axis["top"].major_ticklabels.set_visible(False)
    
    if show_water_ppm is True:
        pass
    else:
        ax_ppm.axis["right"].major_ticklabels.set_visible(False)    
    
    ax.set_xlabel('Position ($\mu$m)')
    
    # Set y-label
    if wholeblock is True:
        if heights_instead is False:
            ax.set_ylabel('Area/Area$_0$')
        else:
            ax.set_ylabel('Height/Height$_0$')            
    else:
        if heights_instead is False:
            ax.set_ylabel('Area (cm$^{-2}$)')
        else:
            ax.set_ylabel('Height (cm$^{-1}$)')

    ax.set_ylim(0, ytop)

    ax.grid()
    return fig, ax, ax_ppm


def plot_3panels_outline(style=None, ytop=1.2, figsize=(6.5, 2.5),
                         shrinker=0.1, heights_instead=False,
                         wholeblock=True, unit='microns'):
    """Outline setup for 3 subplots for 3D profiles"""
    if style is None:
        style = style_lightgreen

    fig, axis3 = plt.subplots(nrows=1, ncols=3)
    fig.set_size_inches(figsize)

    for k in range(3):
        axis3[k].set_ylim(0, ytop)
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
        
    axis3[0].set_xlabel('|| x')
    if unit == 'microns':
        axis3[1].set_xlabel('position ($\mu$m) || y')
    elif unit == 'mm':
        axis3[1].set_xlabel('position (mm) || y')
    else:
        print('unit must = microns or mm')
    axis3[2].set_xlabel('|| z')
    plt.setp(axis3[1].get_yticklabels(), visible=False)
    plt.setp(axis3[2].get_yticklabels(), visible=False)
    return fig, axis3
    

def plot_3panels(positions_microns, area_profiles, lengths=None,
                 styles3=[None, None, None], ytop=1.2, figaxis3=None, 
                 show_line_at_1=True, init=1., 
                 centered=True, unit='microns',
                 percent_error=3., xerror=50., yerror=None,
                 heights_instead=False, wholeblock=True,
                 use_errorbar=False, scale=1.):
    """
    Make 3 subplots for 3D and 3DWB profiles. The position and area profiles
    are passed in lists of three lists for a, b, and c.
    """
    if figaxis3 is None:
        fig, axis3 = plot_3panels_outline(ytop=ytop, wholeblock=wholeblock,
                                          heights_instead=heights_instead,
                                          unit=unit)
    else:
        axis3 = figaxis3

    if lengths is None:
        lengths = np.ones(3)
        for k in range(3):
            lengths[k] = max(positions_microns[k] - min(positions_microns[k]))

    for k in range(3): 
        x = positions_microns[k]
        if unit == 'mm':
            x = np.array(x) / 1000.
        y = np.array(area_profiles[k])

        if len(x) != len(y):
            print('Problem in plot_3panels')
            print('len(x):', len(x))
            print('len(y):', len(y))

        a = lengths[k] / 2.
        pos = x 
        
        current_length = axis3[k].get_xlim()[1]
        if centered is True:
            if current_length < a:
                axis3[k].set_xlim(-a, a)
        else:
            if current_length < lengths[k]:
                axis3[k].set_xlim(0., lengths[k])            

        if show_line_at_1 is True:
            axis3[k].plot([-a, lengths[k]], [init, init], '--k')
            
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
                if yerror is None:
                    yerrorplot = np.array(y*scale) * percent_error/100.
                else:
                    yerrorplot = np.ones_like(pos) * yerror
                axis3[k].errorbar(pos, y*scale, 
                                  xerr=xerror, yerr=yerrorplot, **styles3[k])
            else:
                axis3[k].plot(pos, y*scale, **styles3[k])

    if figaxis3 is None:
        return fig, axis3   

        