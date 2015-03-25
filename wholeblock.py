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
def profile_area_3DWB(final_profile, initial_profile=None, 
                      initial_area_list=None, 
                      initial_area_positions_microns=None,
                      show_plot=True, top=1.2, fig_axis=None):
    """Take final 3D-wholeblock FTIR profile and returns
    a profile of the ratio of the two (A/Ao). Requires information about 
    initial state, so either an initial profile (best) set up either here
    or as the profile.wb_initial profile, or else an initial area list passed
    in with its position. Defaults to making a plot with max y-value 'top'.
    """
    fin = final_profile

    # initial checks
    if len(fin.positions_microns) == 0:
        print 'Need position information'
        return False
    if fin.len_microns is None:
        check = fin.set_len()
        if check is False:
            print 'Need more info to set profile length'
            return False
    if fin.len_microns is None:
        fin.set_len()
    if len(fin.areas_list) == 0:
        print 'making area list for profile'
        fin.make_area_list(show_plot=False)

    # What to normalize to? Priority given to self.wb_initial_profile, then
    # initial_profile passed in here, then initial area list passed in here.
    if fin.wb_initial_profile is not None:
        init = fin.wb_initial_profile
        if init.len_microns != fin.len_microns:
            print 'initial and final lengths must be the same!'
            return False
        # Make sure area lists are populated
        for profile in [init, fin]:
            if len(profile.areas_list) == 0:
                print 'making area list for profile'
                profile.make_area_list(show_plot=False)
        A0 = init.areas_list
        positions0 = init.positions_microns

    elif initial_profile is not None:
        init = initial_profile
        if isinstance(init, pynams.Profile) is False:
            print 'initial_profile argument must be a pynams Profile.'
            print 'Consider using initial_area_list and positions instead'
            return False
        # Make sure area lists are populated
        for profile in [init, fin]:
            if len(profile.areas_list) == 0:
                print 'making area list for profile'
                profile.make_area_list(show_plot=False)
        A0 = init.areas_list
        positions0 = init.positions_microns

    elif initial_area_list is not None:
        if initial_area_positions_microns is None:
            print 'Need initial_area_positions_microns for initial_area_list'
            return False
        A0 = initial_area_list
        positions0 = initial_area_positions_microns
        if len(fin.areas_list) == 0:
            print 'making area list for final profile'
            fin.make_area_list(show_plot=False)
    else:
        print 'Need some information about initial state'
        return False

    # More initial checks
    if len(fin.areas_list) != len(fin.positions_microns):
        print 'area and position lists do not match'
        print 'length areas_list:', len(fin.areas_list)
        print 'length positions list:', len(fin.positions_microns)
        return False    
    if len(A0) < 1:        
        print 'Nothing in initial area list'
        return False
    if len(positions0) < 1:
        print 'Nothing in initial positions list'
        return False
    if len(A0) == 1:
        print 'Using single point to generate initial line'
        A0.extend([A0[0], A0[0]])
        positions0.extend([0, fin.len_microns])
        
    # Use best-fit line through initial values to normalize final data
    p = np.polyfit(positions0, A0, 1)
    normalizeto = np.polyval(p, fin.areas_list)
    wb_areas = fin.areas_list / normalizeto
    leng = fin.len_microns
    
    if show_plot is True:
        if fig_axis is None:
            f, ax = final_profile.plot_area_profile_outline()
        else:
            ax = fig_axis
        ax.set_ylim(0, top)
        ax.set_ylabel('Final area / Initial area')
        style = fin.choose_marker_style()
        ax.plot([-leng/2.0, leng/2.0], [1, 1], **pynams.style_1)
        ax.plot(fin.positions_microns - (leng/2.0), wb_areas, **style)
        if fig_axis is None:
            return wb_areas, f, ax
        else:
            return wb_areas
    else:
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
                 equilibrium_value=None, show_plot=True,
                 erf_or_sum = 'erf', points=100, infinity=5000, 
                 style=None, fig_axis=None):
    """Takes length of profile (microns), time (s), and the 
    log base 10 of diffusivity (m2/s) and returns 1D unit diffusion profile.
    Optional input: 
    - initial concentration (default 1), 
    - whether diffusion is in or out of sample (default out), 
    - whether to use error functions or infinite sums (default erf), 
    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
    - whether to plot results (default False, so no plot)
    - which figure axis to plot onto
    - points sets how many points to calculate in profile. Default is 100.
    - what 'infinity' is if using infinite sum approximation
    - curve plotting style dictionary
    """
    if time_seconds < 0:
        print 'no negative time'
        return
        
    # change to meters for calculation or else numbers too large --> nan
    twoA = length_microns / 1e6       # length in meters
    a = twoA / 2.
    x = np.linspace(-a, a, points)  # positions in meters
    x_microns = x * 1e6
    D = (10.0**log10D_m2s)            # diffusivity in meters^2/second
    t = time_seconds                  # time in seconds    
    sqrtDt = (D*t)**0.5

    if erf_or_sum == 'sum':       
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
        unit_profile = xsum * 4. / np.pi
        
    elif erf_or_sum == 'erf':
        unit_profile = ((scipy.special.erf((a+x)/(2*sqrtDt))) + 
                   (scipy.special.erf((a-x)/(2*sqrtDt)))
                   - 1) 
    else:
        print ('erf_or_sum must be set to either "erf" for python built-in ' +
               'error function approximation (defaul) or "sum" for infinite ' +
               'sum approximation with infinity=whatever, defaulting to ' + 
               str(infinity))
        return
    
    if show_plot is True:
        if fig_axis is None:
            fig, ax = plt.subplots()
            ax.grid()
            ax.set_ylim(0, 1.2)
            xside_microns = length_microns/2.
            ax.set_xlim(-xside_microns, xside_microns)
            ax.set_ylabel('C/C0')
            ax.set_xlabel('position ($\mu$m)')
        else:
            ax = fig_axis
            
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}

        ax.plot(x_microns, unit_profile, **style)

        if fig_axis is None:
            return unit_profile, fig, ax
    else:
        return unit_profile

def diffusion_3D(list_of_3_lengths, time_seconds, list_of_log10D_m2s, 
                 initial_value=None, in_or_out='out', 
                 equilibrium_value=None, show_plot=True,
                 erf_or_sum = 'erf', points=100, infinity=5000, style=None,
                 show_1Dplots=False):
    """ Takes list of three lengths (all in microns), time (s), and 
    diffusivity (m2/s) and returns 3D matrix of internal concentrations.
    Similar options as diffusion_1D.
    """
    # First create 3 1D profiles, 1 in each direction
    profiles = []
    for k in range(3):
        prof = diffusion_1D(list_of_3_lengths[k], time_seconds, 
                            list_of_log10D_m2s[k], initial_value, 
                            in_or_out, equilibrium_value, show_1Dplots,
                            erf_or_sum, points, infinity, style)
        profiles.append(prof)
    # Then multiply them together to get a 3D matrix
    v = np.ones((points, points, points))
    for d in range(0, points):
        for e in range(0, points):
            for f in range(0, points):
                v[d][e][f] = profiles[0][d]*profiles[1][e]*profiles[2][f]
                
    if show_plot is True:
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}
        fig, axis = plt.subplots(nrows=1, ncols=3)

        mid = points/2        
        aslice = v[:, mid][:, mid]
        bslice = v[mid][:, mid]
        cslice = v[mid][mid]
        sliceprofiles = [aslice, bslice, cslice]

        for k in range(3):
            a = list_of_3_lengths[k]
            x = np.linspace(-a, a, points)            
            axis[k].plot(x, sliceprofiles[k], **style)
            axis[k].set_xlim(-a, a)
            axis[k].grid()
            axis[k].set_ylim(0, 1.2)

        axis[0].set_ylabel('C/C$_0$')
        axis[1].set_xlabel('position along slice ($\mu$m)')
        plt.setp(axis[1].get_yticklabels(), visible=False)
        plt.setp(axis[2].get_yticklabels(), visible=False)
        plt.tight_layout
        fig.autofmt_xdate()
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
