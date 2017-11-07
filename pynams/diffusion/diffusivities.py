# -*- coding: utf-8 -*-
"""
Created on Tue May 05 08:16:08 2015
@author: Ferriss

class Diffusivities() groups together temperatures and diffusivities
for use in plotting directly onto Arrhenius diagrams and solving 
for activation energies and pre-exponential components

Also provides basic plotting and functions for creating Arrhenius diagrams and 
handling groups of diffusivities that get plotted on them, including 
determining activation energies and pre-exponential terms and calculating
the expected diffusivity at a given temperature.

"""
from __future__ import print_function, division, absolute_import
from pynams import styles
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from uncertainties import ufloat
import sys

GAS_CONSTANT = 0.00831 # kJ/mol K


def solve_Ea_D0(log10D_list, celsius_list, printout=True):
    """
    Takes lists of diffusivities as log10 in m2/s and associated 
    temperatures in celsius. Returns activation energy Ea in kJ/mol K and D0 
    in m2/s. The errors on the individual diffusivities are not included.
    """
    T = np.array(celsius_list) + 273.15
    x = 1.E4 / T
    y = np.array(log10D_list)

    try:
        # If I don't add in a very low weighted extra number, the covariance
        # matrix, and hence the error, comes out as infinity. The actual
        # fitting results don't change. The very low number in this case I'm 
        # getting from sys.float_info.
        x_extra = np.concatenate((x, x[-1:]), axis=0)
        y_extra = np.concatenate((y, y[-1:]), axis=0)    
        weights = list(np.ones(len(x))) + [sys.float_info.epsilon]
        try:
            fit_extra, cov_extra = np.polyfit(x_extra, 
                                              y_extra, 1, 
                                              w=weights, 
                                              cov=True)
        except np.linalg.linalg.LinAlgError as err:
            return 0, 0                            
            
        if cov_extra[0][0] > 0:
            Ea_error = cov_extra[0][0]
        else:
            Ea_error = 0.
        if cov_extra[1][1] > 0:
            D0_error = cov_extra[1][1]
        else:
            D0_error = 0.
            
        Ea = -ufloat(fit_extra[0], Ea_error) * 2.303 * GAS_CONSTANT * 1.E4
        D0 = 10.**ufloat(fit_extra[1], D0_error)
        
    except ValueError:
        fit_extra = np.polyfit(x_extra, y_extra, 1)
        Ea = fit_extra[0] * 2.303*GAS_CONSTANT*1.E4
        D0 = 10.**fit_extra[1]
    except IndexError:
        # probably no points or unequal number of points
        Ea = 0
        D0 = 0
    except TypeError:
        Ea = 0
        D0 = 0
    
    if printout is True:
        print(Ea, D0)
        
    return Ea, D0


def whatIsD(Ea, D0, celsius, printout=True):
    """
    Takes activation energy in kJ/mol, D0 in m2/s and 
    temperature in celsius. Returns log10 diffusivity in m2/s
    """
    T = celsius + 273.15
    try:
        D0 = D0.n
    except AttributeError:
        pass
    
    try:
        Ea = Ea.n
    except AttributeError:
        pass
    
    D = D0 * np.exp(-Ea / (GAS_CONSTANT * T))
    D = np.log10(D)
        
    if printout is True:
        print('log10 D at ', celsius, 'C: ', '{:.1f}'.format(D), ' in m2/s')
    return D


class Diffusivities():
    def __init__(self, 
                description=None, 
                celsius = [[], [], [], []],
                log10D = [[], [], [], []],
                activation_energy_kJmol = [],
                D0_m2s = [],
                sample=None,  
                ):
        """
        Handy groupings of experimentally determined diffusivities.
        
        Data is grouped as lists in order or orientation:
        parallel a, parallel b, parallel c, and unoriented.
        """
        self.description = description
        self.celsius = celsius
        self.log10D = log10D
        self.sample = sample
        self.activation_energy_kJmol = activation_energy_kJmol
        self.D0_m2s = D0_m2s


    def fill_in_data(self, df, mech, percentpv, paper=None):
        """
        See usage example in literaturevalues.py in this folder.
        
        Input: grouped pandas dataframe, 
               mechanism name (bulk, [Mg], etc.)
               the paper name used to make the pandas groupby
               default assumes paper = self.description
        
        Fills in temperature and diffusivities from a spreadsheet.
        """
        self.log10D = [[], [], [], []]
        self.celsius = [[], [], [], []]
        if paper is None:
            paper = self.description
        for idx, orient in enumerate(['a', 'b', 'c', 'u']):
            group = (paper, mech, float(percentpv), orient)
            try:
                self.log10D[idx] = list(df.get_group(group)['log10D'])                
            except KeyError:
                self.log10D[idx] = []
            
            try:
                self.celsius[idx] = list(df.get_group(group)['celsius'])
            except KeyError:
                self.celsius[idx] = []
       
        
    def solve_Ea_D0(self, printout=True):
        """
        Returns and saves as attributes the
        activation energy in kJ/mol and D0 (log10) in m2/s 
        """ 
        self.activation_energy_kJmol = []
        self.D0_m2s = []
        for idx in range(4):
            logD = self.log10D[idx]
            celsius = self.celsius[idx]
            Ea, D0 = solve_Ea_D0(logD, celsius, printout=printout)
            self.activation_energy_kJmol.append(Ea)
            self.D0_m2s.append(D0)


    def whatIsD(self, celsius, orient='ALL', printout=True):
        """ 
        Takes temperature in celsius. 
        
        Returns log10 diffusivity in m2/s.
        """ 
        D = [] 
        for idx, direction in enumerate(['a', 'b', 'c', 'not oriented']):
            try:
                Ea = self.activation_energy_kJmol[idx]
                D0 = self.D0_m2s[idx]
                
            except IndexError:
                print('\nThese should have 4 items')
                print('Remember 0 placeholder for no measurements')
                print('activation energy')
                print(self.activation_energy_kJmol)
                print('D0')
                print(self.D0_m2s)
                return None
            
            try:
                if (D0 != 0) and (Ea != 0):
                    xD = whatIsD(Ea, D0, celsius, printout=False)
                else: 
                    xD = 'D0 or Ea = 0'
            except TypeError:
                xD = 'Type Error'
                
            if (printout is True) and (orient == 'ALL'):
                print(xD, '||', direction)
            D.append(xD)
        
        orient_list = ['a', 'b', 'c', 'u', 'not oriented', 'ALL']
        if orient in orient_list:
            idx = orient_list.index(orient)
            if idx == 4:
                idx = 3
            if (printout is True) and (orient != 'ALL'):
                print(D[idx], '||', orient)
        else:
            print('Accepted orient values:')
            print('ALL default', orient_list)
        return D
        

    def plotD(self, axes, orient='ALL', plotdata=True,
              offset_celsius=0, plotline=True, extrapolate_line=False,
              legend_add=False, 
              legend_handle=None, style=styles.style_points.copy(), 
              error_color=None, lower_legend_by=-2.0,
              style_line=styles.style_1, label=None):
        """
        Takes axes handle for Arrhenius diagram created by 
        Arrhenius_outline() below and plots data (plotdata=True) and 
        best-fit line (plotline=True) for specified orientations in
        orient (default='ALL'). 
        
        If extrapolate_line=True, the line will be extended all the way 
        to the edges of the plot. 
        
        You can also add it to the legend (default) or not (legend_add=False).
        If you are adding to the legend, you *must* also pass in the 
        legend_handle, which was returned by Arrhenius_outline()
        
        """
        if orient == 'ALL':
            orient_list = list(range(4))
        else:
            iorient = styles.get_iorient(orient)    
            orient_list = [iorient]
                  
        for iorient in orient_list:
            celsius = self.celsius[iorient]
            logD = self.log10D[iorient]

            if orient == 'ALL':
                label = None
                
            if label is None:
                if iorient == 0:
                    label = '|| [100]'
                elif iorient == 1:
                    label = '|| [010]'
                elif iorient == 2:
                    label = '|| [001]'
                elif iorient == 3:
                    label = 'not oriented'                   
  
            if legend_add is True and legend_handle is None:
                print(self.description)
                print('Need legend_handle for legend')
                return
               
            if (len(celsius) == 0) and (self.celsius_all is not None):
                celsius = self.celsius_all

            if (len(celsius) == 0):
                continue
    
            if logD is None:
                continue
    
            if len(celsius) != len(logD):
                print('\n', self.description)
                print('temp and logD not the same length')
                continue
    
            # change temperature scale                   
            x = []
            for k in range(len(celsius)):
                x.append(1.0e4 / (celsius[k] + offset_celsius + 273.15))

            if plotline is True:    
                Tmin = min(x)
                Tmax = max(x)
                T = [Tmin, Tmax]
                p = np.polyfit(x, logD, 1)
                
                if extrapolate_line is True:
                    extrap = axes.get_xlim()
                    axes.plot(extrap,np.polyval(p, extrap), **style_line)
                else:
                    axes.plot(T,np.polyval(p, T), **style_line)
                   
            if label is not None:
                style['label'] = label
            if plotdata is True:
                axes.plot(x, logD, **style)
            
        if legend_add is True:
            self.add_to_legend(axes, legend_handle, style=style,
                               style_line=style_line, plotline=plotline,
                               lower_legend_by=lower_legend_by)

            
    def add_to_legend(self, axes, legend_handle, lower_legend_by=-2.0,
                      orient=None, plotline=False,
                      ncol=2, style=styles.style_points.copy(), 
                      style_line=styles.style_1, label=None):
        """
        Take a figure axis and its list of legend handles 
        and adds information to it
        """
        if label is None:
            if self.description is None:
                print('Need label or self.description to make legend')
                return
            else:
               descript = self.description
        else:
            descript = label
        style['label'] = descript

        add_marker = mlines.Line2D([], [], **style) 
        
        legend_handle.append(add_marker)
        
        low = axes.get_xlim()[0]
        high = axes.get_xlim()[1]
        bottom = axes.get_ylim()[0]
        main_legend = plt.legend(handles=legend_handle, 
                                 numpoints=1, ncol=ncol, 
                                 bbox_to_anchor=(low, bottom, high-low, 
                                                 lower_legend_by),
                                 bbox_transform=axes.transData, 
                                 mode='expand')
        plt.gca().add_artist(main_legend)
        return legend_handle

generic = Diffusivities()
generic.basestyle = {'marker' : 's', 'color' : 'black', 'alpha' : 0.5,
                     'markersize' : 8, 'linestyle': 'none'}


def Arrhenius_outline(xlow=6., xhigh=11., ybottom=-18., ytop=-8.,
                      celsius_labels = np.arange(0, 2000, 100),
                      shrink_axes_to_fit_legend_by = 0.3, make_legend=False,
                      lower_legend_by=-2., ncol=2):
    """
    Make Arrhenius diagram outline. 
    
    Returns figure, axis, legend handle.
    
    low, high, top, and bottom set the x and y axis limits. 

    celsius_labels sets where to make the temperature tick marks.
    
    If you have issues with the legend position or overlap with main diagram,
    play with the numbers for shrink_legend_by and lower_legend_by
    
    ncol sets the number of columns in the legend.
    """
    fig = plt.figure()
    ax = SubplotHost(fig, 1,1,1)
    ax_celsius = ax.twin()
    parasite_tick_locations = 1e4/(celsius_labels + 273.15)
    ax_celsius.set_xticks(parasite_tick_locations)
    ax_celsius.set_xticklabels(celsius_labels)
    fig.add_subplot(ax)
    ax.axis["bottom"].set_label("10$^4$/Temperature (K$^{-1}$)")
    ax.axis["left"].set_label("log$_{10}$diffusivity (m$^{2}$/s)")
    ax_celsius.axis["top"].set_label("Temperature ($\degree$C)")
    ax_celsius.axis["top"].label.set_visible(True)
    ax_celsius.axis["right"].major_ticklabels.set_visible(False)
    ax.set_xlim(xlow, xhigh)
    ax.set_ylim(ybottom, ytop)
    ax.grid()
    
    # main legend below
    if make_legend is True:
        legend_handles_main = []
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height*shrink_axes_to_fit_legend_by, 
                         box.width, box.height*(1.0-shrink_axes_to_fit_legend_by)])
        main_legend = plt.legend(handles=legend_handles_main, numpoints=1, 
                                 ncol=ncol, 
                                 bbox_to_anchor=(xlow, ybottom, xhigh-xlow, 
                                                 lower_legend_by),
                                 bbox_transform=ax.transData, mode='expand')
        plt.gca().add_artist(main_legend)
    else:
        legend_handles_main = None
    return fig, ax, legend_handles_main


def Arrhenius_add_line(fig_ax, Ea, D0, xlow=6.0, xhigh=10.0, 
                       style={'color' : 'k', 'linestyle' : '-'}):
    """
    Takes figure axis from above, Ea activation energy in kJ/mol, D0 in m2/s
    Plots Arrhenius line from 1E4/T = xlow to xhigh
    """
    T = 1E4 / np.linspace(xlow, xhigh) 
    log10D = np.log10(D0) - (Ea/(2.303 * GAS_CONSTANT * T))
    fig_ax.plot(1E4 / T, log10D, **style)

