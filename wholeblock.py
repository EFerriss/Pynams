# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:42:04 2015

@author: Ferriss

For applying the whole-block method of Ferriss et al. 2015 American 
Mineralogist to obtain FTIR concentration profiles, diffusivities, and 
internal concentration distributions following diffusion experiments.

Also includes 1D and 3D slice diffusion calculations for a rectangular
parallelepiped

Requires profile classes set up in pynams.py

MIT license
"""

#%% Generate 3D whole-block area and water profiles
def wb_area_profile(initial_profile, final_profile, make_plot=True):
    """Take initial area profile and final area profile and returns
    a profile of the ratio of the two (A/Ao). Defaults to making a plot"""
    pass

def wb_water_profile(initial_profile, final_profile, water_ppmH2O_initial,
                     make_plot=True):
    """Take initial and final profiles and initial water content.
    Returns the whole-block water concentration profile. Default makes
    a plot showing A/Ao and water on parasite y-axis"""
    pass

#%%  Generate diffusion profiles in 1D, 3D, and 3D through the whole-block
def diffusion_1D(length, time_seconds, D_m2_per_s, 
                 initial_value=None, in_or_out='out', erf_or_infsum='erf',
                 equilibrium_value=None, make_plot=False):
    """Takes length of profile (microns), time (s), and 
    diffusivity (m2/s) and returns diffusion profile.
    Optional input: 
    - initial concentration (default 1), 
    - whether diffusion is in or out of sample (default out), 
    - whether to use error functions or infinite sums (default erf), 
    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
    - whether to plot results (default False, so no plot)
    """
    pass

def diffusion_3D(list_of_3_lengths, time_seconds, D_m2_per_s,
                 initial_value=None, in_or_out='out', erf_or_infsum='erf',
                 equilibrium_value=None, make_plot=True):
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
    pass

def diffusion_3DWB(list_of_3_lengths, time_seconds, D_m2_per_s, 
                 list_of_3_raypaths=['z', 'y', 'y'],
                 initial_value=None, in_or_out='out', erf_or_infsum='erf',
                 equilibrium_value=None, make_plot=True):
    """ Takes list of three lengths (x, y, and z; all in microns), time (s), 
    diffusivity (m2/s), and list of three ray paths (e.g., ['y', 'x', 'x']). 
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

# Inversion of whole-block data to obtain concentration information
def invert_wholeblock(list_of_wb_profiles, grid_xyz,
              symmetry_constraint=True, smoothness_constraint=True,
              rim_constraint=True, rim_value=None,
              weighting_factor_lambda=0.2, show_residuals_plot=True):
    """Takes a list of three whole-block concentration profiles (either A/Ao 
    or water ok but must be consistent for all three) in three orthogonal 
    directions and list of three integers to indicate number of divisions
    in each direction. Returns matrix of values in each grid cell. 
    Default plot showing residuals for how well results match the whole-block
    observations."""
    pass