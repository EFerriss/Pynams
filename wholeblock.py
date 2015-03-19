# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:42:04 2015

@author: Ferriss

For applying the whole-block method of Ferriss et al. 2015 American 
Mineralogist to obtain FTIR concentration profiles, diffusivities, and 
internal concentration distributions following diffusion experiments.

Also includes 1D and 3D slice diffusion calculations for a rectangular
parallelepiped

Input requires profile and spectrum classes set up in pynams.py

MIT license
"""
import numpy as np
import scipy
import pynams
import math
import matplotlib.pyplot as plt
import uncertainties
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.transforms as mtransforms

#%% Generate 3D whole-block area and water profiles
def plot_wbprofile_outline(profile):
    """Makes outline of figure for whole-block profile"""
    f, ax = profile.plot_area_profile_outline()
    ax.set_ylim(0, 1.2)
    ax.set_ylabel('Final area / Initial area')
    return f, ax

def profile_area_3DWB(initial_profile, final_profile, show_plot=True):
    """Take initial area profile and final area profile and returns
    a profile of the ratio of the two (A/Ao). Defaults to making a plot"""
    init = initial_profile
    fin = final_profile

    # Check lengths and area lists exist
    for x in [init, fin]:
        if x.len_microns is None:
            x.set_len()
        if len(x.areas_list) == 0:
            print 'making area list!'
            x.make_area_list(show_plot=False)
        if len(x.positions_microns) == 0:
            print 'No position information!!'
            return False
    
    if init.len_microns != fin.len_microns:
        print 'initial and final lengths must be the same!'
        return False
    else:
        leng = init.len_microns

    # Use best-fit line through initial values to normalize final data
    p = np.polyfit(initial_profile.positions_microns, 
                   initial_profile.areas_list, 1)
    normalizeto = np.polyval(p, fin.areas_list)
    wb_areas = fin.areas_list / normalizeto
    
    if show_plot is True:
        f, ax = plot_wbprofile_outline(initial_profile)
        style = fin.choose_marker_style()
        ax.plot([-leng/2.0, leng/2.0], [1, 1], **pynams.style_1)
        ax.plot(fin.positions_microns - (leng/2.0), wb_areas, **style)        
    return wb_areas

def profile_water_3DWB(initial_profile, final_profile, 
                       water_ppmH2O_initial=None,
                       show_plot=True):
    """Take initial and final profiles and initial water content.
    Returns the whole-block water concentration profile. Default makes
    a plot showing A/Ao and water on parasite y-axis"""
    init = initial_profile
    fin = final_profile
    
    if water_ppmH2O_initial is not None:
        w0 = water_ppmH2O_initial
    else:    
        if init.sample.initial_water is not None:
            w0 = init.sample.initial_water
        elif fin.sample.initial_water is not None:
            w0 = fin.sample.initial_water
        else:
            print 'Need initial water content.'
            return False
        
    
    wb_areas = profile_area_3DWB(init, fin, show_plot=False)
    water = wb_areas * w0
    
    if show_plot is True:
        # Use a parasite y-axis to show water content
        fig = plt.figure()
        ax_areas = SubplotHost(fig, 1,1,1)
        fig.add_subplot(ax_areas)
        area_tick_marks = np.arange(0, 100, 0.2)
        ax_areas.set_yticks(area_tick_marks)
        ax_water = ax_areas.twin()
        ax_water.set_yticks(area_tick_marks)
        if isinstance(w0, uncertainties.Variable):
            ax_water.set_yticklabels(area_tick_marks*w0.n)
        else:
            ax_water.set_yticklabels(area_tick_marks*w0)
        ax_areas.axis["bottom"].set_label('Position ($\mu$m)')
        ax_areas.axis["left"].set_label('Final area / Initial area')
        ax_water.axis["right"].set_label('ppm H$_2$O')
        ax_water.axis["top"].major_ticklabels.set_visible(False)
        ax_water.axis["right"].major_ticklabels.set_visible(True)
        ax_areas.set_ylim(0, 1.2)
        leng = init.len_microns
        ax_areas.set_xlim(-leng/2.0, leng/2.0)
        ax_areas.grid()
            
        style = fin.choose_marker_style()
        ax_areas.plot([-leng/2.0, leng/2.0], [1, 1], **pynams.style_1)
        ax_areas.plot(fin.positions_microns-leng/2.0, wb_areas, **style)
        return water, fig, ax_areas
    else:
        return water

#%%  Generate diffusion profiles in 1D, 3D, and 3DWB
def diffusion_1D(length_microns, time_seconds, log10D_m2s, 
                 initial_value=None, in_or_out='out', 
                 equilibrium_value=None, show_plot=True, fig_axis=None,
                 points=100, infinity=5000):
    """Takes length of profile (microns), time (s), and the 
    log base 10 of diffusivity (m2/s) and returns diffusion profile.
    Optional input: 
    - initial concentration (default 1), 
    - whether diffusion is in or out of sample (default out), 
    - whether to use error functions or infinite sums (default erf), 
    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
    - whether to plot results (default False, so no plot)
    - which figure axis to plot onto
    - points sets how many points to calculate in profile. Default is 100.
    """
    if time_seconds < 0:
        print 'no negative time'
        return
        
#    if equilibrium_value is not None:
#        # add test to ensure equilib is a number
#        equilib = equilibrium_value
#    elif in_or_out == 'in':
#        equilib = 1.0
#    elif in_or_out == 'out':
#        equilib = 0.0
#    else: 
#        print 'Not sure what the equilibrium value is'
#
#    if initial_value is not None:
#        pass
#    elif in_or_out == 'in':
#        initial_value = 0.
#    elif in_or_out == 'out':
#        initial_value = 1.0
#    else:
#        print 'Not sure what the initial value is'

    # change to meters for calculation or else numbers too large --> nan
    twoA = length_microns / 1e6       # length in meters
    a = twoA / 2.
    x = np.linspace(-a, a, points)  # positions in meters
    x_microns = x * 1e6
    D = (10.0**log10D_m2s)            # diffusivity in meters^2/second
    t = time_seconds                  # time in seconds    
       
    xsum = np.zeros_like(x)    
    for n in range(infinity):
       # positive number that converges to 1
        xadd1 = ((-1.)**n) / ((2.*n)+1.)        
        # time conponent
        xadd2 = np.exp(
                        (-D * (((2.*n)+1.)**2.) * (np.pi**2.) * t) / 
                        (twoA**2.) 
                        )                        
        # There the position values come in to create the profile
        xadd3 = np.cos(
                        ((2.*n)+1.) * np.pi * x / twoA
                        )        
        xadd = xadd1 * xadd2 * xadd3
        xsum = xsum + xadd
        
    profile = xsum * 4. / np.pi
    
    if show_plot is True:
        if fig_axis is None:
            fig, ax = plt.subplots()
        else:
            ax = fig_axis
        ax.plot(x_microns, profile)
        ax.grid()
        ax.set_ylim(0, 1.2)
        xside_microns = length_microns/2.
        ax.set_xlim(-xside_microns, xside_microns)
        ax.set_ylabel('C/C0')
        ax.set_xlabel('position ($\mu$m)')
        
    return profile

def diffusion_3D(list_of_3_lengths, time_seconds, log10D_m2s,
                 initial_value=None, in_or_out='out', erf_or_infsum='erf',
                 equilibrium_value=None, make_plot=True, points=100):
    """ Takes list of three lengths (all in microns), time (s), and 
    diffusivity (m2/s) and returns 3D matrix of internal concentrations.
    This is the 3D-slice approach described in Ferriss et al. 2015
    Optional input: 
    - initial concentration (default 1), 
    - whether diffusion is in or out of sample (default out), 
    - whether to use error functions or infinite sums (default erf), 
    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
    - whether to plot results as three profiles through center (default True)
    """
    profiles = []
    for length in list_of_3_lengths:
        prof = diffusion_1D(length, time_seconds, log10D_m2s, initial_value, 
                     in_or_out, erf_or_infsum, equilibrium_value, 
                     False, None, points)
        profiles.append(prof)
    
    v = np.ones_like(profiles[0])
    for d in range(0, points):
        v[d] = profiles[0][d]
#            for f=1:points+1
#                v(d,e,f)=xerf(d)*yerf(e)*zerf(f);

    fig, ax = plt.subsubplots()
    x = np.linspace(-length, length, points)
    ax.plot(x, v)

    return v

def diffusion_3DWB(list_of_3_lengths, time_seconds, log10D_m2s, 
                 list_of_3_raypaths=['c', 'b', 'b'],
                 initial_value=None, in_or_out='out', erf_or_infsum='erf',
                 equilibrium_value=None, make_plot=True):
    """ Takes list of three lengths (x, y, and z; all in microns), time (s), 
    diffusivity (m2/s), and list of three ray paths (e.g., ['b', 'a', 'a']). 
    Returns 3 profiles of whole-block concentrations,
    the concentrations measured through the center of the block.
    This is the 3D-WB approach described in Ferriss et al. 2015
    Optional input: 
    - initial concentration (default 1), 
    - whether diffusion is in or out of sample (default out), 
    - whether to use error functions or infinite sums (default erf), 
    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
    - whether to plot results as three profiles through center (default True)
    """
    pass

#%% Solve for diffusivity in 1D, 3D, or 3D-WB
def getD_1D(profile, time_seconds, initial_value=None,
                         in_or_out='out', erf_or_infsum='erf',
                         equilibrium_value=None):
    """ Take a profile, time in seconds, and optional other information
    Return diffusivity using least-squares solution.
    For whole-block 1D-approximation, set initial_value = plateau value
    """
    pass

def getD_3Dslice(list_of_3_profiles, time_seconds, initial_value=None,
                         in_or_out='out', erf_or_infsum='erf',
                         equilibrium_value=None):
    """Take a list of 2 or 3 profiles measured through the center of the 
    sample, time in seconds, and optional other information.
    Return diffusivities from least-squares solution.
    """
    pass

def getD_3DWB(list_of_3_wb_profiles, time_seconds, make_plot=True,
                           initial_value=None, in_or_out='out', 
                           erf_or_infsum='erf', equilibrium_value=None):
    """Takes a list of three whole-block concentration profiles (either A/Ao 
    or water ok but must be consistent for all three) in three orthogonal 
    directions and experiment time in seconds. 
    Returns diffusivities in each direction. 
    Default to plotting."""
    pass


#%% Group profiles (to include: maps!) together as whole-block unit
class WholeBlock():
    profiles = None
    # generated by check_and_setup(self)
    directions = None
    raypaths = None
    initial_profiles = None
    # optional for diffusion work
    time_seconds = None
    diffusivity_log10_m2s = None
    solubility = None
            
    def setupWB(self):
        """Sets up and checks WholeBlock instance
        - Check that wb_profile_list contains a list of three (3) profiles
        - Generate list of initial profiles
        - Generate list of profile directions
        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
        - Generate list of ray paths
        - Verify three ray path directions are compatible with directions list
        """
        pass
    
    def plot_3panels(self):
        """Plot whole-block profiles values"""
        pass
    
    def diffuse(self, time_seconds=None, log10D_m2s=None, 
                initial_value=None, in_or_out='out', erf_or_infsum='erf', 
                equilibrium_value=None, make_plot=True):
        """Applies diffusion_3DWB function to WholeBlock instance"""
        pass
#        diffusion_3DWB(self.wb_profile_list, time_seconds, log10D_m2s, 
#                 self.wb_raypath_list, initial_value, in_or_out, 
#                 erf_or_infsum, equilibrium_value, make_plot)

    def solveD(self, time_seconds=None, make_plot=True, initial_value=None, 
               in_or_out='out', erf_or_infsum='erf', equilibrium_value=None):
        """Applies getD_3DWB function to WholeBlock instance"""
        pass

    def invert(self, grid_xyz, symmetry_constraint=True, 
               smoothness_constraint=True, rim_constraint=True, 
               rim_value=None, weighting_factor_lambda=0.2, 
               show_residuals_plot=True):
        """Takes a list of three whole-block concentration profiles (either A/Ao 
        or water ok but must be consistent for all three) in three orthogonal 
        directions and list of three integers to indicate number of divisions
        in each direction. Returns matrix of values in each grid cell. 
        Default plot showing residuals for how well results match the whole-block
        observations."""
        pass
