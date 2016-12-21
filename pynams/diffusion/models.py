# -*- coding: utf-8 -*-
"""
Created on Tue May 05 08:16:08 2015
@author: Ferriss

Diffusion in 1 and 3 dimensions with and without path-integration

THIS MODULE ASSUMES THAT ALL CONCENCENTRATIONS ARE ALREADY NORMALIZED TO 1.

### 1-dimensional diffusion ###
Simplest function call is diffusion1D(length, diffusivity, time)
    Step 1. Create lmfit parameters with params = params_setup1D
            (Here is where you set what to vary when fitting)
    Step 2. Pass these parameters into diffusion1D_params(params)
    Step 3. Plot with plot_diffusion1D
With profiles in styles, use profile.plot_diffusion() and fitD()

### 3-dimensional diffusion without path integration: 3Dnpi ###
Simplest: diffusion3Dnpi(lengths, D's, time) to get a figure
    Step 1. Create parameters with params = params_setup3D
    Step 2. Pass parameters into diffusion3Dnpi_params(params) to get profiles
            Returns full 3D matrix v, sliceprofiles, then slice positions
    Step 3. Plot with styles.plot_3panels(slice positions, slice profiles)

### Whole-Block: 3-dimensional diffusion with path integration: 3Dwb ###
    Step 1. Create parameters with params = params_setup3D 
            Same as for non-path integrated 3D.
    Step 2. Pass parameters into diffusion3Dwb(params)
    Requires raypath information, unlike 3Dnpi

"""
import pynams.styles as styles
import numpy as np
import pynams.diffusion.lmfit as lmfit
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

#%% 1D Diffusion in a thin slab; Eq 4.18 in Crank, 1975
def diffusionThinSlab(log10D_m2s, thickness_microns, max_time_hours=2000, 
                      infinity=200, timesteps=300):
    """ Takes log10 of the diffusivity D in m2/s, thickness in microns,
    and maximum time in hours, and returns time array in hours and corresponding 
    curve of C/C0, the concentration divided by the initial concentration
    """
    t_hours = np.linspace(0., max_time_hours, timesteps)
    t_seconds = t_hours * 3600.
    L_meters = thickness_microns / 2.E6
    cc = np.zeros_like(t_seconds)
    D_m2s = 10.**log10D_m2s

    for idx in range(len(t_seconds)):
        infsum = 0.
        for n in range(infinity):
            exponent_top = D_m2s * (((2.*n)+1.)**2.) * (np.pi**2.) * t_seconds[idx] 
            exponent_bot = 4. * (L_meters**2)
            exponent = np.exp(-1. * exponent_top / exponent_bot)
            addme = 8. * (1. / ((((2.*n)+1)**2.) * (np.pi**2.))) * exponent
            infsum = infsum + addme
        cc[idx] = infsum
    
    return t_hours, cc
    

#%% 1D diffusion profiles
def params_setup1D(microns, log10D_m2s, time_seconds, init=1., fin=0.,
                   vD=True, vinit=False, vfin=False):
    """Takes required info for diffusion in 1D - length, diffusivity, time,
    and whether or not to vary them - vD, vinit, vfin. 
    Return appropriate lmfit params to pass into diffusion1D_params"""
    params = lmfit.Parameters()
    params.add('microns', microns, False, None, None, None)
    params.add('log10D_m2s', log10D_m2s, vD, None, None, None)
    params.add('time_seconds', time_seconds, False, None, None, None)
    params.add('initial_unit_value', init, vinit, None, None, None)
    params.add('final_unit_value', fin, vfin, None, None, None)
    return params

def diffusion1D_params(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', need_to_center_x_data=True,
                 infinity=100, points=50):
    """Function set up to follow lmfit fitting requirements.
    Requires input as lmfit parameters value dictionary 
    passing in key information as 'length_microns',
    'time_seconds', 'log10D_m2s', and 'initial_unit_value'. 

    If data are None (default), returns 1D unit diffusion profile 
    x_microns and y as vectors of length points (default 50). 

    Optional keywords:    
     - erf_or_sum: whether to use python's error functions (default) 
       or infinite sums
     - whether to center x data
     - points sets how many points to calculate in profile. Default is 50.
     - what 'infinity' is if using infinite sum approximation
     
    If not including data, returns the x vector and model y values.
    With data, return the residual for use in fitting.

    Visualize results with plot_diffusion1D
    """
    # extract important values from parameter dictionary passed in
    p = params.valuesdict()
    L_meters = p['microns'] / 1e6
    t = p['time_seconds']
    D = 10.**p['log10D_m2s']
    initial_value = p['initial_unit_value']
    final_value = p['final_unit_value']

    if initial_value > final_value:
        going_out = True
        solubility = initial_value
        minimum_value = final_value
    else:
        going_out = False        
        solubility = final_value
        minimum_value = initial_value
    
    a_meters = L_meters / 2.
    twoA = L_meters

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

    if going_out is False:
        model = np.ones_like(model) - model

    concentration_range = solubility - minimum_value
    model = (model * concentration_range) + minimum_value

    x_microns = x * 1e6

    # If not including data, just return the model values
    # With data, return the residual for use in fitting.
    if fitting is False:
        return x_microns, model
    return model-data_y_unit_areas

def plot_diffusion1D(x_microns, model, initial_value=None,
                     fighandle=None, axishandle=None, top=1.2,
                     style=None, fitting=False, show_km_scale=False,
                     show_initial=True):
    """Takes x and y diffusion data and plots 1D diffusion profile input"""
    a_microns = (max(x_microns) - min(x_microns)) / 2.
    a_meters = a_microns / 1e3
    
    if fighandle is None and axishandle is not None:
        print 'Remember to pass in handles for both figure and axis'
    if fighandle is None or axishandle is None:
        fig = plt.figure()          
        ax  = SubplotHost(fig, 1,1,1)
        ax.grid()
        ax.set_ylim(0, top)
    else:
        fig = fighandle
        ax = axishandle

    if style is None:
        if fitting is True:
            style = {'linestyle' : 'none', 'marker' : 'o'}
        else:
            style = styles.style_lightgreen

    if show_km_scale is True:
        ax.set_xlabel('Distance (km)')
        ax.set_xlim(0., 2.*a_meters/1e3)
        x_km = x_microns / 1e6
        ax.plot((x_km) + a_meters/1e3, model, **style)
    else:                
        ax.set_xlabel('position ($\mu$m)')
        ax.set_xlim(-a_microns, a_microns)
        ax.plot(x_microns, model, **style)

    if initial_value is not None and show_initial is True:
        ax.plot(ax.get_xlim(), [initial_value, initial_value], '--k')

    ax.set_ylabel('Unit concentration or final/initial')
    fig.add_subplot(ax)

    return fig, ax

def diffusion1D(length_microns, log10D_m2s, time_seconds, init=1., fin=0.,
                erf_or_sum='erf', show_plot=True, 
                fighandle=None, axishandle=None,
                style=None, need_to_center_x_data=True,
                infinity=100, points=100, top=1.2, show_km_scale=False):
    """Simplest implementation.
    Takes required inputs length, diffusivity, and time 
    and plots diffusion curve on new or specified figure. 
    Optional inputs are unit initial value and final values. 
    Defaults assume diffusion 
    out, so init=1. and fin=0. Reverse these for diffusion in.
    Returns figure, axis, x vector in microns, and model y data."""
    params = params_setup1D(length_microns, log10D_m2s, time_seconds, 
                            init, fin,
                            vD=None, vinit=None, vfin=None)
                            
    x_microns, model = diffusion1D_params(params, None, None, 
                                          erf_or_sum, need_to_center_x_data, 
                                          infinity, points)

    fig, ax = plot_diffusion1D(x_microns, model, initial_value=init, 
                               fighandle=fighandle, axishandle=axishandle,
                               style=style, fitting=False, 
                               show_km_scale=show_km_scale)
    
    return fig, ax, x_microns, model


#%% 3-dimensional diffusion parameter setup
def params_setup3D(microns3, log10D3, time_seconds, 
                   initial=1., final=0., isotropic=False, slowb=False,
                   vD=[True, True, True], vinit=False, vfin=False):
    """Takes required info for diffusion in 3D without path averaging and 
    return appropriate lmfit params.
    
    Returning full matrix and 
    slice profiles in one long list for use in fitting

    """
    params = lmfit.Parameters()
    params.add('microns3', microns3, False, None, None, None)
    params.add('log10Dx', log10D3[0], vD[0], None, None, None)
    params.add('time_seconds', time_seconds, False, None, None, None)
    params.add('initial_unit_value', initial, vinit, None, None, None)
    params.add('final_unit_value', final, vfin, None, None, None)

    if isotropic is True:
        params.add('log10Dy', log10D3[1], vD[1], None, None, 'log10Dx')
        params.add('log10Dz', log10D3[2], vD[2], None, None, 'log10Dx')
    elif slowb is True:
        params.add('log10Dy', log10D3[1], vD[1], None, None, 'log10Dx - 1.')
        params.add('log10Dz', log10D3[2], vD[2], None, None, 'log10Dx')
    else:            
        params.add('log10Dy', log10D3[1], vD[1], None, None, None)
        params.add('log10Dz', log10D3[2], vD[2], None, None, None)

    return params

def diffusion3Dnpi_params(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', centered=True, 
                 infinity=100, points=50):
    """ Diffusion in 3 dimensions in a rectangular parallelipiped.
    Takes params - Setup parameters with params_setup3D.
    General setup and options similar to diffusion1D_params.
    
    Returns complete 3D concentration
    matrix v, slice profiles, and 
    positions of slice profiles.
    
    ### NOT COMPLETELY SET UP FOR FITTING JUST YET ###
    """   
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
    L3_microns = np.array(p['microns3'])
    t = p['time_seconds']
    init = p['initial_unit_value']
    vary_init = [params['initial_unit_value'].vary]
    fin = p['final_unit_value']
    vary_fin = [params['final_unit_value'].vary]
    log10D3 = [p['log10Dx'], p['log10Dy'], p['log10Dz']]
    vary_D = [params['log10Dx'].vary, 
              params['log10Dy'].vary, 
              params['log10Dz'].vary]

    # If initial values > 1, scale down to 1 to avoid blow-ups later
    going_out = True
    scale = 1.
    if init > 1.0:
        scale = init
        init = 1.
    if init < fin:
        going_out = False
    
    if init > fin:
        minimum_value = fin
    else:
        minimum_value = init
    
    if going_out is False:        
        # I'm having trouble getting diffusion in to work simply, so this
        # is a workaround. The main effort happens as diffusion going in, then
        # I subtract it all from 1.
        init, fin = fin, init
        
    # First create 3 1D profiles, 1 in each direction
    xprofiles = []    
    yprofiles = []
    kwdict = {'points' : points}
    
    for k in range(3):
        p1D = lmfit.Parameters()
        p1D.add('microns', L3_microns[k], False)
        p1D.add('time_seconds', t, params['time_seconds'].vary)
        p1D.add('log10D_m2s', log10D3[k], vary_D[k])
        p1D.add('initial_unit_value', init, vary_init)
        p1D.add('final_unit_value', fin, vary_fin)
        
        x, y = diffusion1D_params(p1D, **kwdict)

        xprofiles.append(x)
        yprofiles.append(y)
                                      
    # Then multiply them together to get a 3D matrix
    # I should figure out how to do this without the for-loops
    v = np.ones((points, points, points))
    for d in range(0, points):
        for e in range(0, points):
            for f in range(0, points):
                v[d][e][f] = yprofiles[0][d]*yprofiles[1][e]*yprofiles[2][f]

    v = v * scale

    if going_out is False:
        v = np.ones((points, points, points)) - v
        v = v + np.ones_like(v)*minimum_value

    mid = int(points/2.)
    
    aslice = v[:, mid][:, mid]
    bslice = v[mid][:, mid]
    cslice = v[mid][mid]
    sliceprofiles = [aslice, bslice, cslice]

    slice_positions_microns = []
    for k in range(3):
        a = L3_microns[k] / 2.
        if centered is False:            
            x = np.linspace(0, a*2., points)
        else:
            x = np.linspace(-a, a, points)
        slice_positions_microns.append(x)
          
    # Returning full matrix and 
    # slice profiles in one long list for use in fitting
    sliceprofiles = [aslice, bslice, cslice]
    
    if fitting is False:
        return v, sliceprofiles, slice_positions_microns
    else:
        ### Still need to set up residuals! ###
        residuals = np.zeros_like(sliceprofiles)
        return residuals

def diffusion3Dnpi(lengths_microns, log10Ds_m2s, time_seconds, points=50,
                    initial=1, final=0., top=1.2, plot3=True, centered=True,
                    styles3=[None]*3, figaxis3=None):
        """
        Required input:
        list of 3 lengths, list of 3 diffusivities, and time 
        
        Optional input:
        initial concentration (1), final concentration (0), 
        whether to plot output (plot3=True), maximum y limit on plot (top=1.2),
        and number of points to use during calculation (points=50)
        
        If plot3=True (default), plots results and returns:
        1. f, figure of plot of 3D non-path-averaged diffusion profiles.
        2. ax, 3 axes of figure
        3. v, 3D matrix of diffusion
        4. x, list of 3 sets of x values plotted
        5. y, list of 3 sets of y values plotted
        
        If plot3=False, returns only v, x, y
        """
        params = params_setup3D(lengths_microns, log10Ds_m2s, time_seconds,
                                initial=initial, final=final)
                                                                
        v, y, x = diffusion3Dnpi_params(params, points=points, centered=False)

        if centered is True:
            for idx in xrange(3):
                x[idx] = x[idx] - (lengths_microns[idx] / 2.)

        if plot3 is True:
            try:
                if figaxis3 is None:
                    f, ax = styles.plot_3panels(x, y, figaxis3=figaxis3, 
                                                top=top, centered=centered,
                                                styles3=styles3, 
                                                lengths=lengths_microns)
                    return f, ax, v, x, y
                else:
                    styles.plot_3panels(x, y, top=top, centered=centered, 
                                        styles3=styles3, figaxis3=figaxis3,
                                        lengths=lengths_microns)
                    return v, x, y
            except(TypeError):
                print
                print 'TypeError: problem in plot_3panels()'
                
        else:
            return v, x, y
            
#%% 3D whole-block: 3-dimensional diffusion with path integration
def diffusion3Dwb_params(params, data_x_microns=None, data_y_unit_areas=None, 
                          raypaths=None, erf_or_sum='erf', show_plot=True, 
                          fig_ax=None, style=None, need_to_center_x_data=True,
                          infinity=100, points=50, show_1Dplots=False):
    """ Diffusion in 3 dimensions with path integration.
    Requires setup with params_setup3Dwb
    """
    if raypaths is None:
        print 'raypaths must be in the form of a list of three abc directions'
        return

    # v is the model 3D array of internal concentrations
    ### Need to add in all the keywords ###
    v, sliceprofiles, slicepositions = diffusion3Dnpi_params(params, 
                    points=points, erf_or_sum=erf_or_sum,
#                    need_to_center_x_data=need_to_center_x_data
                    )

    
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
    L3 = p['microns3']
    
    wb_profiles = [wbA, wbB, wbC]
    wb_positions = []
    for k in range(3):
        a = L3[k] / 2.
        x_microns = np.linspace(0., 2.*a, points)
        wb_positions.append(x_microns)
        
    if show_plot is True:
        if style is None:
            style = [None, None, None]
            for k in range(3):
                style[k] = styles.style_lightgreen

        if fig_ax is None:
            f, fig_ax = styles.plot_3panels(wb_positions, wb_profiles, L3, style)
        else:
            styles.plot_3panels(wb_positions, wb_profiles, L3, style, 
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

def diffusion3Dwb(lengths_microns, log10Ds_m2s, time_seconds, raypaths,
                   initial=1., final=0., top=1.2, points=50., show_plot=True,
                   figax=None, isotropic=False):
        """Takes list of 3 lengths, list of 3 diffusivities, and time.
        Returns plot of 3D path-averaged (whole-block) diffusion profiles"""
        params = params_setup3D(lengths_microns, log10Ds_m2s, time_seconds,
                                initial=initial, final=final)

        return params
        
        x, y = diffusion3Dwb_params(params, raypaths=raypaths, show_plot=False,
                                    points=points)

        if show_plot is True:
            if figax is None:
                f, ax = styles.plot_3panels(x, y, top=top)
                return f, ax, x, y
            else: 
                styles.plot_3panels(x, y, top=top, figaxis3=figax)
        return x, y

        
