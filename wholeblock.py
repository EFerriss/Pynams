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
from scipy.optimize import curve_fit
import pynams
import math
import matplotlib.pyplot as plt
import uncertainties
from uncertainties import ufloat
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.transforms as mtransforms
import lmfit

# - plotting profiles in three panels
# - Generating whole-block area and water profiles
# - Diffusion in 1D, 3D, and 3D-WB
# - Forward model curve-fitting to diffusion profiles
#    - Write objective functions that work with lmfit
#    - Call those objective functions with, e.g., getD_1D()

#%% 3D Plot setup
def plot_3panels_outline(style=None, top=1.2):
    """Outline setup for 3 subplots for 3D profiles"""
    if style is None:
        style = {'color' : 'lightgreen', 'linewidth' : 4}
    fig, axis3 = plt.subplots(nrows=1, ncols=3)
    for k in range(3):
        axis3[k].set_ylim(0, top)
        axis3[k].grid()
    axis3[0].set_ylabel('C/C$_0$')
    axis3[1].set_xlabel('position along slice ($\mu$m)')
    plt.setp(axis3[1].get_yticklabels(), visible=False)
    plt.setp(axis3[2].get_yticklabels(), visible=False)
    plt.tight_layout
    fig.autofmt_xdate()
    return fig, axis3
    
def plot_3panels(positions_microns, area_profiles, lengths, 
                 style=None, top=1.2, figaxis3=None, show_line_at_1=True,
                 centered=True):
    """Make 3 subplots for 3D and 3DWB profiles. The position and area profiles
    are passed in lists of three lists for a, b, and c.
    Positions are assumed to start at 0 and then are centered.    
    """
    if figaxis3 is None:
        fig, axis3 = plot_3panels_outline(style, top)
    else:
        axis3 = figaxis3        

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
            axis3[k].plot([-a, a], [1., 1.], '-k')
            
        if style is None:
            if len(x) > 45:
                axis3[k].plot(x-a, y)
            else:
                axis3[k].plot(x-a, y, 'o')
        else:
            axis3[k].plot(x-a, y, **style)

    if figaxis3 is None:
        return fig, axis3

#%% Generate 3D whole-block area and water profiles
def make_3DWB_area_profile(final_profile, initial_profile=None, 
                           initial_area_list=None, 
                           initial_area_positions_microns=None,
                           show_plot=True, top=1.2, fig_ax=None):
    """Take final 3D-wholeblock FTIR profile and returns
    a profile of the ratio of the two (A/Ao). A/A0 is also saved in the 
    profile attribute wb_areas.
    Requires information about initial state, so either an initial 
    profile (best) as the profile's initial_profile attribute,
    or else an initial area list passed in with its position. 
    Defaults to making a plot with max y-value 'top'.
    Note initial_area_positions_microns is assumed to start at 0 and then 
    gets shifted for the fit.
    """
    fin = final_profile
    leng = fin.set_len()

    # If whole-block areas are already made, use them. 
    # Otherwise make them.
    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
        wb_areas = fin.wb_areas
    else:
        print fin.wb_areas
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
        if fin.initial_profile is not None:
            
            init = fin.initial_profile
            if init.len_microns != fin.len_microns:
                print 'initial and final lengths must be the same!'
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas_list) == 0:
                    print profile.profile_name
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
        p = np.polyfit(positions0-(leng/2.), A0, 1)
        normalizeto = np.polyval(p, fin.areas_list)
        wb_areas = fin.areas_list / normalizeto
         
        # Save whole-block areas as part of profile
        fin.wb_areas = wb_areas    
    
    if show_plot is True:
        if fig_ax is None:
            f, ax = final_profile.plot_area_profile_outline()
        else:
            ax = fig_ax
        ax.set_ylim(0, top)
        ax.set_ylabel('Final area / Initial area')
        style = fin.choose_marker_style()
        ax.plot([-leng/2.0, leng/2.0], [1, 1], **pynams.style_1)
        ax.plot(fin.positions_microns - (leng/2.0), wb_areas, **style)
        if fig_ax is None:
            return wb_areas, f, ax
        else:
            return wb_areas
    else:
        return wb_areas


def make_3DWB_water_profile(final_profile, water_ppmH2O_initial=None,
                            initial_profile=None, 
                            initial_area_list=None, 
                            initial_area_positions_microns=None,
                            show_plot=True, top=1.2, fig_ax=None):
    """Take a profile and initial water content.
    Returns the whole-block water concentration profile based on
    the profile's attribute wb_areas. If wb_areas have not been made, 
    some initial profile information and various options are passed
    to make_3DWB_area_profile().
    Default makes a plot showing A/Ao and water on parasite y-axis
    """
    fin = final_profile
    init = initial_profile

    # Set initial water
    if water_ppmH2O_initial is not None:
        w0 = water_ppmH2O_initial
    else:
        if fin.sample is not None:
            if fin.sample.initial_water is not None:
                w0  = fin.sample.initial_water
        elif init is not None:
            if init.sample is not None:
                if init.sample.initial_water is not None:
                    w0 = init.sample.initial_water
        else:
            print 'Need initial water content.'
            return False
    
    # Set whole-block areas
    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
        wb_areas = fin.wb_areas
    else:  
        wb_areas = make_3DWB_area_profile(fin, initial_profile, 
                                          initial_area_list, 
                                          initial_area_positions_microns)
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
        ax_areas.grid()
        ax_areas.set_ylim(0, 1.2)
        if fin.len_microns is not None:
            leng = fin.len_microns
        else:
            leng = fin.set_len()
        ax_areas.set_xlim(-leng/2.0, leng/2.0)
            
        style = fin.choose_marker_style()
        ax_areas.plot([-leng/2.0, leng/2.0], [1, 1], **pynams.style_1)
        ax_areas.plot(fin.positions_microns-leng/2.0, wb_areas, **style)
        return water, fig, ax_areas
    else:
        return water

#%% 
#
#
# Writing diffusion equations as suitable ojective functions that can be 
# called by lmfit.minimize.
#
#
def diffusion_1D(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', show_plot=True, fig_axis=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=50):
    """Function set up to follow lmfit fitting requirements.
    Requires input as lmfit parameters value dictionary 
    passing in key information as 'length_microns',
    'time_seconds', 'log10D_m2s', and 'initial_unit_value'. 
    
    Here is an example of the setup:
    params = lmfit.Parameters()
    #           (Name,                 Value,   Vary,   Min,  Max,  Expr)
    params.add_many(('length_microns',  1000,   False,  0.0, None,  None),
                ('time_seconds',      3600*10., False,  0.0,  None,  None),
                ('log10D_m2s',         -12.,    True,   None, None,  None),
                ('initial_unit_value',  1.,     False,  0.0,   1.0, None))

    Fitting is a profile function fitDiffusivity().
    For fitting, include data_x_microns and data_y_unit_areas. 
    If need_to_center_x_data is True (default), x data starts at 0 and
    will get moved. Also, lmfit screws up if you try setting a minimum
    on log10D_m2s.

    Return values
    If data are None (default), returns 1D unit diffusion profile 
    x_microns and y as vectors of length points (default 50). 

    Optional keywords:    
     - erf_or_sum: whether to use python's error functions (default) 
       or infinite sums
     - show_plot: whether to plot results (default True, so plot)
     - fig_ax: which figure axis to plot onto (default None makes new fig)
     - style: curve plotting style dictionary
#    - whether diffusion is in or out of sample (default out), 
#    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
     - points sets how many points to calculate in profile. Default is 50.
     - what 'infinity' is if using infinite sum approximation
    """
    # extract important values from parameter dictionary passed in
    p = params.valuesdict()
    L_meters = p['length_microns'] / 1e6
    t = p['time_seconds']
    D = 10.**p['log10D_m2s']
         # I will make this keyword optional at some point
    initial_value = p['initial_unit_value']
    a_meters = L_meters / 2.

    if t < 0:
        print 'no negative time'
        return           

    # Fitting to data or not? Default is not
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        if len(data_x_microns) == len(data_y_unit_areas):
            fitting = True
        else:
            print 'x and y data must be the same length'
            print 'x', len(data_x_microns)
            print 'y', len(data_y_unit_areas)
        
    # x is in meters and assumed centered around 0
    if fitting is True:
        # Change x to meters and center it
        x = np.array(data_x_microns) / 1e6
        if need_to_center_x_data is True:
            x = x - a_meters
    else:
        x = np.linspace(-a_meters, a_meters, points)
    
    if erf_or_sum == 'infsum':
        xsum = np.zeros_like(x)
        for n in range(infinity):
            xadd =  ((((-1.)**n) / ((2.*n)+1.))    + 
                    (np.exp((-D * (((2.*n)+1.)**2.) * (np.pi**2.) * t) / 
                        (L_meters**2.)))    + 
                    (np.cos(((2.*n)+1.) * np.pi * x / L_meters)))
            xsum = xsum + xadd        
        model = xsum * 4. / np.pi        

    elif erf_or_sum == 'erf':
        sqrtDt = (D*t)**0.5
        model = ((scipy.special.erf((a_meters+x)/(2*sqrtDt))) + 
                   (scipy.special.erf((a_meters-x)/(2*sqrtDt))) - 1) 

    else:
        print ('erf_or_sum must be set to either "erf" for python built-in ' +
               'error function approximation (defaul) or "sum" for infinite ' +
               'sum approximation with infinity=whatever, defaulting to ' + 
               str(infinity))
        return False

    # Revisit for in_or_out
    model = model * initial_value

    x_microns = x * 1e6

    if show_plot is True:
        a_microns = a_meters * 1e6
        if fig_axis is None:
            fig, fig_axis = plt.subplots()
            fig_axis.grid()
            fig_axis.set_ylim(0, 1.2)
            fig_axis.set_xlim(-a_microns, a_microns)
            fig_axis.set_ylabel('C/C$_0$')
            fig_axis.set_xlabel('position ($\mu$m)')

        if style is None:
            if fitting is True:
                style = {'linestyle' : 'none', 'marker' : 'o'}
            else:
                style = {'color' : 'lightgreen', 'linewidth' : 4}
            
        fig_axis.plot([-a_microns, a_microns], [initial_value, initial_value], 
                      '-k')
        fig_axis.plot(x_microns, model, **style)

    # If not including data, just return the model values
    # With data, return the residual for use in fitting.
    if fitting is False:
        return x_microns, model
    return model-data_y_unit_areas


def diffusion_3D(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', show_plot=True, fig_ax=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=50, show_1Dplots=False):
    """ Diffusion in 3 dimensions without path integration.
    lmfit parameter list input must include:
    length_microns_a, length_microns_b, length_microns_c, 
    time_seconds, initial unit value, and 
    diffusivities in log10 m2/s logDx, logDy, logDz
    
    Example parameter setup for input here:
    params = lmfit.Parameters()
    #  (Name,           Value,      Vary,   Min,  Max,  Expr)
    params.add('microns_twoa', L3[0], False, 0.0, None, None)
    params.add('microns_twob', L3[1], False, 0.0, None, None)
    params.add('microns_twoc', L3[2], False, 0.0, None, None)
    params.add('initial_unit_value_a', 1., False, 0.0, None, None)
    params.add('initial_unit_value_b', 1., True, 0.0, None, None)
    params.add('initial_unit_value_c', 1., False, 0.0, None, None)
    params.add('log10Dx', D3[0], True, 0.0, None, None)
    params.add('log10Dy', D3[1], True, 0.0, None, None)
    params.add('log10Dz', D3[2], False, 0.0, None, None)
    params.add('time_seconds', t, False, 0.0, None, None)    
    """
    # Fitting to data or not? Default is not
    # Add appropriate x and y data to fit
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_data = np.array(data_x_microns)
        y_data = np.array(data_y_unit_areas)
        if np.shape(x_data) == np.shape(y_data):
            fitting = True
            print 'fitting to data'
        else:
            print 'x and y data must be the same shape'
            print 'x', np.shape(x_data)
            print 'y', np.shape(y_data)

    p = params.valuesdict()
    L3_microns = np.array([p['microns_twoa'], p['microns_twob'], 
                          p['microns_twoc']])
    t = p['time_seconds']
    init = [p['initial_unit_value_a'], 
            p['initial_unit_value_b'],
            p['initial_unit_value_c']]
    vary_init = [params['initial_unit_value_a'].vary, 
                 params['initial_unit_value_b'].vary,
                 params['initial_unit_value_c'].vary]
    log10D3 = [p['log10Dx'], p['log10Dy'], p['log10Dz']]
    vary_D = [params['log10Dx'].vary, 
              params['log10Dy'].vary, 
              params['log10Dz'].vary]

    # First create 3 1D profiles, 1 in each direction
    xprofiles = []    
    yprofiles = []
    kwdict = {'show_plot' : show_1Dplots, 'points' : points}
    
    for k in range(3):
        p1D = lmfit.Parameters()
        p1D.add('length_microns', L3_microns[k], False)
        p1D.add('time_seconds', t, params['time_seconds'].vary)
        p1D.add('log10D_m2s', log10D3[k], vary_D[k])
        p1D.add('initial_unit_value', init[k], vary_init[k])
        
        x, y = diffusion_1D(p1D, **kwdict)
        xprofiles.append(x)
        yprofiles.append(y)
                                      
    # Then multiply them together to get a 3D matrix
    v = np.ones((points, points, points))
    for d in range(0, points):
        for e in range(0, points):
            for f in range(0, points):
                v[d][e][f] = yprofiles[0][d]*yprofiles[1][e]*yprofiles[2][f]

    mid = points/2        
    aslice = v[:, mid][:, mid]
    bslice = v[mid][:, mid]
    cslice = v[mid][mid]
    sliceprofiles = [aslice, bslice, cslice]
                
    if show_plot is True:
        if fig_ax is None:
            f, fig_ax = plot_3panels_outline()
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}
        positions = []
        for k in range(3):
            a = L3_microns[k] / 2.
            x = np.linspace(0, a*2., points)
            positions.append(x)
        plot_3panels(positions, sliceprofiles, L3_microns, style=style,
                     figaxis3=fig_ax)
            
    # Returning full matrix and 
    # slice profiles in one long list for use in fitting
    sliceprofiles = [aslice, bslice, cslice]
    
    if fitting is False:
        return v, sliceprofiles
    else:
        ### Still need to set up residuals! ###
        residuals = np.zeros_like(sliceprofiles)
        return residuals


def diffusion_3DWB(params, data_x_microns=None, data_y_unit_areas=None, 
                 raypaths=None, erf_or_sum='erf', show_plot=True, 
                 fig_ax=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=200, show_1Dplots=False):
    """ Diffusion in 3 dimensions with path integration.
    lmfit parameter list input must include:
    length_microns_a, length_microns_b, length_microns_c, 
    time_seconds, initial unit value, 
    diffusivities in log10 m2/s logDx, logDy, logDz, 
    
    Also must pass in a keyword list of raypaths in an order consistent
    with the length directions as, e.g., ['c', 'c', 'b']
    
    Example parameter setup for input here:
    params = lmfit.Parameters()
    #  (Name,           Value,      Vary,   Min,  Max,  Expr)
    params.add('microns_twoa', L3[0], False, 0.0, None, None)
    params.add('microns_twob', L3[1], False, 0.0, None, None)
    params.add('microns_twoc', L3[2], False, 0.0, None, None)
    params.add('initial_unit_value_a', 1., False, 0.0, None, None)
    params.add('initial_unit_value_b', 1., True, 0.0, None, None)
    params.add('initial_unit_value_c', 1., False, 0.0, None, None)
    params.add('log10Dx', D3[0], True, 0.0, None, None)
    params.add('log10Dy', D3[1], True, 0.0, None, None)
    params.add('log10Dz', D3[2], False, 0.0, None, None)
    params.add('time_seconds', t, False, 0.0, None, None)
    """
    if raypaths is None:
        print 'raypaths must be in the form of a list of three abc directions'
        return

    # v is the model 3D array of internal concentrations
    ### Need to add in all the keywords ###
    v, sliceprofiles = diffusion_3D(params, show_plot=False, points=points)

    # Fitting to data or not? Default is not
    # Add appropriate x and y data to fit
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_array = np.array(data_x_microns)
        y_array = np.array(data_y_unit_areas)
        if np.shape(x_array) == np.shape(y_array):
            print 'fitting to data'
            fitting = True
        else:
            print 'x and y data must be the same shape'
            print 'x', np.shape(x_array)
            print 'y', np.shape(y_array)
            
    # Whole-block measurements can be obtained through any of the three 
    # planes of the whole-block, so profiles can come from one of two ray path
    # directions. These are the planes.
    raypathA = v.mean(axis=0)
    raypathB = v.mean(axis=1)
    raypathC = v.mean(axis=2)

    # Specify whole-block profiles in model
    mid = points/2
    if raypaths[0] == 'b':
        wbA = raypathB[:, mid]
    elif raypaths[0] == 'c':
        wbA = raypathC[:, mid]       
    else:
        print 'raypaths[0] for profile || a must be "b" or "c"'
        return
        
    if raypaths[1] == 'a':
        wbB = raypathA[:, mid]
    elif raypaths[1] == 'c':
        wbB = raypathC[mid]       
    else:
        print 'raypaths[1] for profile || b must be "a" or "c"'
        return

    if raypaths[2] == 'a':
        wbC = raypathA[mid]
    elif raypaths[2] == 'b':
        wbC = raypathB[mid]       
    else:
        print 'raypaths[2] for profile || c must be "a" or "b"'
        return

    p = params.valuesdict()
    L3 = [p['microns_twoa'], p['microns_twob'], p['microns_twoc']]
#    if fitting is True:
#        wb_profiles = data_y_unit_areas
#        wb_positions = np.array(data_x_microns)
#    else:
    wb_profiles = [wbA, wbB, wbC]
    wb_positions = []
    for k in range(3):
        a = L3[k] / 2.
        x_microns = np.linspace(0., 2.*a, points)
        wb_positions.append(x_microns)
        
    if show_plot is True:
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}

        if fig_ax is None:
            f, fig_ax = plot_3panels(wb_positions, wb_profiles, L3, style)
        else:
            plot_3panels(wb_positions, wb_profiles, L3, style, 
                         figaxis3=fig_ax)                         

    if fitting is False:        
        return wb_positions, wb_profiles
    
    if fitting is True:
        # Return residuals 
        y_model = []
        y_data = []
        residuals = []
        for k in range(3):
            for pos in range(len(x_array[k])):
                # wb_positions are centered, data are not
                microns = x_array[k][pos]
                # Find the index of the full model whole-block value 
                # closest to the data positions
                idx = (np.abs(wb_positions[k]-microns).argmin())
                
                model = wb_profiles[k][idx]
                data = y_array[k][pos]
                res = model - data
                
                y_model.append(model)
                y_data.append(data)
                residuals.append(res)                
        return residuals

#%% Group profiles together as whole-block unit
class WholeBlock():
    profiles = []
    # generated by setupWB below
    directions = None
    raypaths = None
    initial_profiles = None
    lengths = None
    # optional for diffusion work and plotting
    style_base = None
    time_seconds = None
    diffusivities_log10_m2s = None
                
    def setupWB(self):
        """Sets up and checks WholeBlock instance
        - Check that profiles list contains a list of three (3) profiles
        - Generate list of initial profiles
        - Generate list of profile directions
        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
        - Generate list of ray paths
        - Verify three ray path directions are compatible with directions list
        """
        if len(self.profiles) != 3:
            print 'For now, only a list of 3 profiles is allowed'
            return    
        d = []
        r = []
        ip = []
        L = []
        for prof in self.profiles:
            if isinstance(prof, pynams.Profile) is False:
                print 'Only profiles objects allowed in profile list!'
                return
            d.append(prof.direction)
            r.append(prof.raypath)
            ip.append(prof.initial_profile)
            L.append(prof.set_len())
        self.directions = d
        self.raypaths = r
        self.initial_profiles = ip 
        self.lengths = L
        return

    def plot_3panels(self, fig_ax3=None):
        """Plot data on three panels"""
        if ((self.directions is None) or (self.raypaths is None) or
            (self.initial_profiles is None) or (self.lengths is None)):
                self.setupWB()

        positions = []
        areas = []
        for prof in self.profiles:
            positions.append(prof.positions_microns)
            if prof.wb_areas is None:
                make_3DWB_area_profile(prof, show_plot=False)
            areas.append(prof.wb_areas)

        if fig_ax3 is None:
            f, ax = plot_3panels(positions, areas, self.lengths, 
                             style=self.style_base)
            return f, ax
        else:
            plot_3panels(positions, areas, self.lengths, figaxis3=fig_ax3,
                             style=self.style_base)

    def show_diffusion(self, time_seconds=None, list_of_log10D_m2s=None, 
                initial_value=None, in_or_out='out', erf_or_sum='erf', 
                equilibrium_value=None, show_plot=True, 
                show_slice=False, style=None, wb_or_3Dnpi='wb', 
                fig_ax=None, points=100):
        """Applies 3-dimensionsal diffusion equations to instance shape.
        Requires lengths, time in seconds, and three diffusivities either
        explicitly passed here or as attributes of the WholeBlock object.
        Assuming whole-block diffusion (wb_or_3Dnpi='wb') but could also 
        do 3D non-path-integrated ('npi')
        """
        # initial checks
        if ((self.directions is None) or (self.raypaths is None) or
            (self.initial_profiles is None)):
                self.setupWB()

        if time_seconds is None:
            if self.time_seconds is not None:
                time_seconds = self.time_seconds
            else:
                print 'Need time information'
                return        

        if list_of_log10D_m2s is None:
            if self.diffusivities_log10_m2s is not None:
                list_of_log10D_m2s = self.diffusivities_log10_m2s
            else:
                print 'Need list of three log10 diffusivities in m2/s'
                return

        if (wb_or_3Dnpi != 'npi') and (wb_or_3Dnpi != 'wb'):
            print 'wb_or_3Dnpi only takes "wb" or "npi"'
            return

        if self.lengths is None:
            self.setupWB()
            
        if self.raypaths is None:
            self.setupWB
        # end initial checks and setup

        # Set up parameters to pass into equations
        L3 = self.lengths
        D3 = list_of_log10D_m2s        
        params = lmfit.Parameters()
        # (Name, Value, Vary, Min, Max, Expr)
        params.add('microns_twoa', L3[0], False, None, None, None)
        params.add('microns_twob', L3[1], False, None, None, None)
        params.add('microns_twoc', L3[2], False, None, None, None)
        params.add('initial_unit_value_a', 1., False, None, None, None)
        params.add('initial_unit_value_b', 1., False, None, None, None)
        params.add('initial_unit_value_c', 1., False, None, None, None)
        params.add('log10Dx', D3[0], True, None, None, None)
        params.add('log10Dy', D3[1], True, None, None, None)
        params.add('log10Dz', D3[2], True, None, None, None)
        params.add('time_seconds', time_seconds, False, None, None, None)

        # Set up the plot
        if show_plot is True:
            if fig_ax is None:
                fig, fig_ax = plot_3panels_outline()
            
        kws = {'show_plot' : show_plot, 'fig_ax' : fig_ax,
               'style' : style, 'points' : points}

        # send to diffusion equation
        # diffusion line plots get made here
        if wb_or_3Dnpi == 'npi':
            diffusion_3D(params, **kws)
        else:
            R3 = self.raypaths
            kws['raypaths'] = R3
            diffusion_3DWB(params, **kws)

        # Add data on top of diffusion line in plot
        if show_plot is True:
            self.plot_3panels(fig_ax3=fig_ax)

        return params
       
    def fitD(self, initial_unit_values=[1., 1., 1.], 
             vary_initials=['False', 'False', 'False'], 
             guesses=[-13., -13., -13.], show_plot=True,
             vary_D=['True', 'True', 'True'],
             approx_with_1D=False, wb_or_3Dnpi='wb', polyorder=1,
             show_initial_guess=True, style_initial=None,
             style_final={'color' : 'red'}, points=200):
        """Forward modeling to determine diffusivities in three dimensions 
        from whole-block data. 
        """
        if wb_or_3Dnpi != 'wb' and wb_or_3Dnpi != 'npi':
            print 'wb_or_3Dnpi can only be wb or npi'
            return
            
        # Plot setup and initial guess lines
        if show_plot is True:
            fig, ax3 = plot_3panels_outline()
        else: 
            ax3 = None

        # Set up parameters and keywords to go into diffusion equations
        # This will plot the initial guess diffusivities if show_initial_guess
        # is True
        dict_plotting = {'show_plot' : show_initial_guess, 'fig_ax' : ax3,
                 'wb_or_3Dnpi' : wb_or_3Dnpi, 'style' : style_initial,
                 'points' : points}
        dict_fitting = {'show_plot' : False, 'points' : points}

        params = self.show_diffusion(**dict_plotting)

        # x and y are the data that will be fit to
        x = []
        y = []
        for prof in self.profiles:
            if prof.positions_microns is None:
                print ''
                print prof.profile_name
                print 'Need to set profile positions'
                return
            if prof.wb_areas is None:
                prof.make_3DWB_water_list(polyorder)
            x.append(prof.positions_microns)
            y.append(prof.wb_areas)

        
        if wb_or_3Dnpi == 'wb':
            dict_fitting['raypaths'] = self.raypaths            
            lmfit.minimize(diffusion_3DWB, params, args=(x, y), 
                           kws=dict_fitting)
            
        elif wb_or_3Dnpi == 'npi':
            diffusion_3D(params, x, y, **dict_fitting)
        
#            lmfit.minimize(diffusion_3D, params, args=(x, y), kws=kwdict)

        # results        
        best_Dx = ufloat(params['log10Dx'].value, 
                         params['log10Dx'].stderr)
        best_Dy = ufloat(params['log10Dy'].value, 
                         params['log10Dy'].stderr)
        best_Dz= ufloat(params['log10Dz'].value, 
                         params['log10Dz'].stderr)
        best_init_a = ufloat(params['initial_unit_value_a'].value, 
                             params['initial_unit_value_a'].stderr)
        best_init_b = ufloat(params['initial_unit_value_b'].value, 
                             params['initial_unit_value_b'].stderr)
        best_init_c = ufloat(params['initial_unit_value_c'].value, 
                             params['initial_unit_value_c'].stderr)      

        self.diffusivities_log10_m2s = [best_Dx.n, best_Dy.n, best_Dz.n]

        # and then after
        if show_plot is True:
            dict_plotting['show_plot'] = True
            dict_plotting['style'] = style_final
            self.show_diffusion(**dict_plotting)
            
        print '\ntime in hours:', params['time_seconds'].value / 3600.
        print '\ninitial unit values:'
        print best_init_a
        print best_init_b
        print best_init_c
        print '\nbestfit log10D in m2/s:'
        print best_Dx
        print best_Dy
        print best_Dz

        return params
    
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
