# -*- coding: utf-8 -*-
"""
Created on Tue May 05 08:16:08 2015
@author: Ferriss

Diffusion in 1 and 3 dimensions with and without path-integration

### 1-dimensional diffusion ###
diffusionThinSlab and diffusion1D for profiles
With pynams profiles, use profile.plot_diffusion() and fitD()

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
from __future__ import print_function, division, absolute_import
import pynams.styles as styles
import numpy as np
import lmfit
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost


#%% 1D Diffusion in a thin slab; Eq 4.18 in Crank, 1975
def diffusionThinSlab(log10D_m2s, thickness_microns, max_time_hours=2000, 
                      infinity=200, timesteps=300):
    """ 
    Eq 4.18 in Crank, 1975
    Takes log10 of the diffusivity D in m2/s, thickness in microns,
    and maximum time in hours, and returns time array in hours and 
    corresponding curve of C/C0, the concentration divided by the 
    initial concentration.
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
                   vD=True, vinit=False, vfin=False, vTime=False):
    """
    Takes required info for diffusion in 1D - length, diffusivity, time,
    and whether or not to vary them - vD, vinit, vfin, vTime.
    
    Return appropriate lmfit params to pass into diffusion1D_params
    See https://lmfit.github.io/lmfit-py/parameters.html
    """
    params = lmfit.Parameters()
    params.add('microns', microns, False, None, None, None)
    params.add('log10D_m2s', log10D_m2s, vD, None, None, None)
    params.add('time_seconds', time_seconds, vTime, None, None, None)
    params.add('initial_unit_value', init, vinit, None, None, None)
    params.add('final_unit_value', fin, vfin, None, None, None)
    return params


def diffusion1D_params(params, 
                       data_x_microns=None, 
                       data_y_unit_areas=None, 
                       erf_or_sum='erf', 
                       centered=True, 
                       symmetric=True,
                       infinity=100, 
                       points=50):
    """
    Function set up to follow lmfit fitting requirements.
    See https://lmfit.github.io/lmfit-py/parameters.html
    
    Requires:
        lmfit parameters value dictionary 
        passing in key information as 'length_microns',
        'time_seconds', 'log10D_m2s', and 'initial_unit_value'. 
        See function params_setup1D

    If data are None (default), returns 1D unit diffusion profile 
    x_microns and y as vectors of length points (default 50). 

    Optional keywords:    
     - erf_or_sum: whether to use python's error functions (default) 
       or infinite sums
     - whether to center x data
     - whether the profile is symmetric or not. If not, change is to the right.
     - points sets how many points to calculate in profile. Default is 50.
     - what 'infinity' is if using infinite sum approximation
     
    If not including data, returns the x vector and model y values.
    With data, return the residual for use in fitting.
    """

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
        return           

    # Fitting to data or not? Default is not
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        if len(data_x_microns) == len(data_y_unit_areas):
            fitting = True
        else:
            print('x and y data must be the same length')
            print('x', len(data_x_microns))
            print('y', len(data_y_unit_areas))
        
    # x is in meters and assumed centered around 0
    if fitting is True:            
        # Change x to meters and center it
        x = np.array(data_x_microns) / 1e6
        x = x - a_meters
    else:
        x = np.linspace(-a_meters, a_meters, points)
    
    # Make the infinite sum
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
            # Where the position values come in to create the profile
            xadd3 = np.cos(
                            ((2.*n)+1.) * np.pi * x / twoA
                            )        
            xadd = xadd1 * xadd2 * xadd3
            xsum = xsum + xadd
            
        model = xsum * 4. / np.pi
    
    # or make the error function
    elif erf_or_sum == 'erf':
        if symmetric is False:
            a_meters = a_meters*2
            x = x*2
            
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
    return model- data_y_unit_areas


def diffusion1D(length_microns, log10D_m2s, time_seconds, init=1., fin=0.,
                erf_or_sum='erf', show_plot=True, 
                style=styles.style_blue, infinity=100, points=100, 
                centered=True, axes=None, symmetric=True,
                maximum_value=1.):
    """
    Simplest implementation for 1D diffusion.
    
    Takes required inputs length, diffusivity, and time 
    and plots diffusion curve on new or specified figure. 
    Optional inputs are unit initial value and final values. 
    Defaults assume diffusion out, so init=1. and fin=0. 
    Reverse these for diffusion in.
    
    Change scale of y-values with maximum_value keyword.
    
    Returns figure, axis, x vector in microns, and model y data.
    """    
    if symmetric is True:
        params = params_setup1D(length_microns, log10D_m2s, time_seconds,
                                init=init, fin=fin)
        x_diffusion, y_diffusion = diffusion1D_params(params, points=points)
        if centered is False:
            a_length = (max(x_diffusion) - min(x_diffusion)) / 2
            x_diffusion = x_diffusion + a_length
    else:
        # multiply length by two
        params = params_setup1D(length_microns*2, log10D_m2s, time_seconds,
                                init=init, fin=fin)
        x_diffusion, y_diffusion = diffusion1D_params(params, points=points)        

        # divide elongated profile in half
        x_diffusion = x_diffusion[int(points/2):]
        y_diffusion = y_diffusion[int(points/2):]
        if centered is True:
            a_length = (max(x_diffusion) - min(x_diffusion)) / 2
            x_diffusion = x_diffusion - a_length 

    if show_plot is True:
        if axes is None:
            fig = plt.figure()          
            ax  = SubplotHost(fig, 1,1,1)
            ax.grid()
            ax.set_ylim(0, maximum_value)
            ax.set_xlabel('position ($\mu$m)')
            ax.set_xlim(min(x_diffusion), max(x_diffusion))
            ax.plot(x_diffusion, y_diffusion*maximum_value, **style)
            ax.set_ylabel('Unit concentration or final/initial')
            fig.add_subplot(ax)
        else:
            axes.plot(x_diffusion, y_diffusion*maximum_value, **style)
            fig = None
            ax = None            
    else:
        fig = None
        ax = None
    
    return fig, ax, x_diffusion, y_diffusion


#%% 3-dimensional diffusion parameter setup
def length_checker(microns3):
    """
    Checks lengths passed in, and returns possibly modified lengths
    """
    if isinstance(microns3, int):
        microns3 = [float(microns3)]*3
    elif isinstance(microns3, float):
        microns3 = [microns3]*3
    elif isinstance(microns3, list):
        if len(microns3) == 1:
            microns3 = microns3*3
        elif len(microns3) == 2:
            new_length = np.mean(microns3)
            print('Warning: assuming 3rd length=', new_length)
            microns3 = microns3.append(new_length)
    else:
        print('only int, float, or list allowed for lengths_microns')
        return False
    return microns3


def D_checker(log10D3):
    """
    Checks diffusivities passed in and returns possibly modified diffusivities.
    """
    if isinstance(log10D3, int):
        log10D3 = [float(log10D3)]*3
    elif isinstance(log10D3, float):
        log10D3 = [log10D3]*3
    elif isinstance(log10D3, list):
        if len(log10D3) == 1:
            log10D3 = log10D3*3
        elif len(log10D3) == 2:
            new_D = np.mean(log10D3)
            print('Warning: assuming 3rd D = log10', new_D)
            log10D3 = log10D3.append(new_D)
    else:
        print('only int, float, or list allowed for log10Ds_m2s')
        return False
    return log10D3


def params_setup3D(microns3, log10D3, time_seconds, 
                   initial=1., final=0., isotropic=False, slowb=False,
                   vD=[True, True, True], vinit=False, vfin=False):
    """
    Requires:
        microns3: a list of 3 lengths in microns
        log10D3: a list of 3 diffusivities in m2/s and converted to log10
        time_seconds: the time in seconds
        
    Optional input:
        initial and final unit values
        whether to force diffusive anisotropy (default: isotropic=False)
        slowb: whether to force [010] to be an order of magnitude slower
        vD: a list of whether to vary diffusivities and which ones
        vinit: whether to vary the initial concentration
        vfin: whether to vary the final concentration
    
    Returns:
        an lmfit parameters object for use in 3D fitting functions
        See https://lmfit.github.io/lmfit-py/parameters.html

    """
    params = lmfit.Parameters()
    params.add('microns3', microns3, False, None, None, None)
    params.add('log10Dx', float(log10D3[0]), vD[0], None, None, None)
    params.add('time_seconds', float(time_seconds), False, None, None, None)
    params.add('initial_unit_value', float(initial), vinit, None, None, None)
    params.add('final_unit_value', float(final), vfin, None, None, None)

    if isotropic is True:
        params.add('log10Dy', float(log10D3[1]), vD[1], None, None, 'log10Dx')
        params.add('log10Dz', float(log10D3[2]), vD[2], None, None, 'log10Dx')
    elif slowb is True:
        params.add('log10Dy', float(log10D3[1]), vD[1], None, None,\
                   'log10Dx - 1.')
        params.add('log10Dz', float(log10D3[2]), vD[2], None, None, 'log10Dx')
    else:            
        params.add('log10Dy', float(log10D3[1]), vD[1], None, None, None)
        params.add('log10Dz', float(log10D3[2]), vD[2], None, None, None)

    return params


def diffusion3Dnpi_params(params, 
                          data_x_microns=None, 
                          data_y_unit_areas=None, 
                          erf_or_sum='erf', 
                          centered=True, 
                          infinity=100, 
                          points=50):
    """ 
    Diffusion in 3 dimensions in a rectangular parallelipiped.
    Requires:
        params - lmfit parameters. Setup parameters with params_setup3D.
        For more on parameters, see 
        https://lmfit.github.io/lmfit-py/parameters.html 
    
    Option inpur aimilar to diffusion1D_params.
    
    Returns:
        v: complete 3D concentration matrix
        y: slice profiles
        x: positions of slice profiles
    
    ### NOT SET UP FOR FITTING YET ###
    """   
    
    fitting = False
    
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_data = np.array(data_x_microns)
        y_data = np.array(data_y_unit_areas)
        
        if np.shape(x_data) == np.shape(y_data):
            fitting = True
            print('fitting to data')
        else:
            print('x and y data must be the same shape')
            print('x', np.shape(x_data))
            print('y', np.shape(y_data))

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

    # going in or out?
    if init < fin: 
        going_out = False
        init, fin = fin, init
    else:
        going_out = True
       
    # First create 3 1D profiles, 1 in each direction
    xprofiles = []    
    yprofiles = []
    kwdict = {'points' : points}
    
    for microns, D, vD in zip(L3_microns, log10D3, vary_D):
        p1D = lmfit.Parameters()
        p1D.add('microns', float(microns), False)
        p1D.add('time_seconds', t, params['time_seconds'].vary)
        p1D.add('log10D_m2s', D, vD)
        p1D.add('initial_unit_value', 1., vary_init)
        p1D.add('final_unit_value', 0., vary_fin)
        
        try:
            x, y = diffusion1D_params(p1D, **kwdict)
        except TypeError:
            print('Problem with x values after diffusion1D')
            print(diffusion1D_params(p1D, **kwdict))
            return

        xprofiles.append(x)
        yprofiles.append(y)
                                      
    # Then multiply them together to get a 3D matrix
    # I should figure out how to do this without the for-loops
    v = np.ones((points, points, points))
    for d in range(0, points):
        for e in range(0, points):
            for f in range(0, points):
                v[d][e][f] = yprofiles[0][d]*yprofiles[1][e]*yprofiles[2][f]

    if going_out is False:
        v = np.ones((points, points, points)) - v

    scale = np.abs(fin-init)
    minimum_value = min([fin, init])
    v = (v * scale) + minimum_value

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
                   initial=1, final=0., ytop=None, show_plot=True, 
                   centered=True, styles3=[None]*3, axes=None, 
                   show_line_at_initial=True):
        """
        Required input:
        list of 3 lengths, list of 3 diffusivities, and time 
        
        Optional input:
        initial unit concentration (default initial=1), 
        final unit concentration (default final=0), 
        whether to plot output (default show_plot=True), 
        maximum y limit on plot (ytop),
        and number of points to use during calculation (default points=50)
        
        If show_plot=True (default), plots results and returns:
        1. f, figure of plot of 3D non-path-averaged diffusion profiles.
        2. ax, 3 axes of figure
        3. v, 3D matrix of diffusion
        4. x, list of 3 sets of x values plotted
        5. y, list of 3 sets of y values plotted
        
        If show_plot=False, returns only v, x, y
        """       
        lengths_microns = length_checker(lengths_microns)
        log10Ds_m2s = D_checker(log10Ds_m2s)
        params = params_setup3D(lengths_microns, log10Ds_m2s, time_seconds,
                                initial=initial, final=final)
        
        v, y, x = diffusion3Dnpi_params(params, points=points, 
                                        centered=centered)

        if show_plot is True:
            try:
                if axes is None:
                    f, ax = styles.plot_3panels(x, y, figaxis3=axes, 
                                        ytop=ytop, centered=centered,
                                        styles3=styles3, init=initial,
                                        lengths=lengths_microns,
                                        show_line_at_1=show_line_at_initial)
                    return f, ax, v, x, y
                else:
                    styles.plot_3panels(x, y, ytop=ytop, centered=centered, 
                                        styles3=styles3, figaxis3=axes,
                                        lengths=lengths_microns,
                                        show_line_at_1=show_line_at_initial,
                                        init=initial)
                    return v, x, y
            except(TypeError):
                print
                print('TypeError: problem in plot_3panels()')
                return v, x, y
                
        else:
            return v, x, y
            
        
#%% 3D whole-block: 3-dimensional diffusion with path integration
def diffusion3Dwb_params(params, data_x_microns=None, data_y_unit_areas=None, 
                          raypaths=None, erf_or_sum='erf', show_plot=True, 
                          fig_ax=None, style=None, need_to_center_x_data=True,
                          infinity=100, points=50, show_1Dplots=False):
    """ 
    Diffusion in 3 dimensions with path integration.
    Requires setup with params_setup3Dwb
    """
    if raypaths is None:
        print('raypaths must be in the form of a list of three abc directions')
        return

    v, y, x = diffusion3Dnpi_params(params, points=int(points), 
                                    centered=False)
    
    # Fitting to data or not? Default is not
    # Add appropriate x and y data to fit
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_array = np.array(data_x_microns)
        y_array = np.array(data_y_unit_areas)
        if np.shape(x_array) == np.shape(y_array):
            print('fitting to data')
            fitting = True
        else:
            print('x and y data must be the same shape')
            print('x', np.shape(x_array))
            print('y', np.shape(y_array))
            
    # Whole-block measurements can be obtained through any of the three 
    # planes of the whole-block, so profiles can come from one of two ray path
    # directions. These are the planes.
    raypathA = v.mean(axis=0)
    raypathB = v.mean(axis=1)
    raypathC = v.mean(axis=2)
    
    # Specify whole-block profiles in model
    mid = int(points/2)
    if raypaths[0] == 'b':
        wbA = raypathB[:, mid]
    elif raypaths[0] == 'c':
        wbA = raypathC[:, mid]       
    else:
        print('raypaths[0] for profile || a must be "b" or "c"')
        return
        
    if raypaths[1] == 'a':
        wbB = raypathA[:, mid]
    elif raypaths[1] == 'c':
        wbB = raypathC[mid]       
    else:
        print('raypaths[1] for profile || b must be "a" or "c"')
        return

    if raypaths[2] == 'a':
        wbC = raypathA[mid]
    elif raypaths[2] == 'b':
        wbC = raypathB[mid]       
    else:
        print('raypaths[2] for profile || c must be "a" or "b"')
        return

    p = params.valuesdict()
    L3 = p['microns3']
    
    wb_profiles = [wbA, wbB, wbC]
    wb_positions = []
    for length in L3:
        x_microns = np.linspace(0., float(length), points)
        wb_positions.append(x_microns)
        
    if show_plot is True:
        if style is None:
            style = [None, None, None]
            for k in range(3):
                style[k] = styles.style_lightgreen

        if fig_ax is None:
            f, fig_ax = styles.plot_3panels(wb_positions, wb_profiles, 
                                            L3, style)
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
                # closest to the data position
                idx = (np.abs(wb_positions[k]-microns).argmin())
                
                model = wb_profiles[k][idx]
                data = y_array[k][pos]
                res = model - data
                
                y_model.append(model)
                y_data.append(data)
                residuals.append(res)                
        return np.array(residuals)


def diffusion3Dwb(lengths_microns, log10Ds_m2s, time_seconds, raypaths,
                   initial=1., final=0., ytop=1.2, points=50, show_plot=True,
                   axes=None, isotropic=False, centered=True,
                   show_line_at_initial=True, styles3=[None]*3):
        """
        Takes list of 3 lengths, list of 3 diffusivities, and time.
        Returns plot of 3D path-averaged (whole-block) diffusion profiles.
        
        Requirements and keywords similar to diffusion.3Dnpi except a list
        of raypaths (thickness directions as 'a', 'b', or 'c')
        """
        lengths_microns = length_checker(lengths_microns)
        log10Ds_m2s = D_checker(log10Ds_m2s)
        params = params_setup3D(lengths_microns, log10Ds_m2s, time_seconds,
                                initial=initial, final=final)
        
        x, y = diffusion3Dwb_params(params, raypaths=raypaths, show_plot=False,
                                    points=int(points))

        if centered is True:
            for idx in range(3):
                x[idx] = x[idx] - (lengths_microns[idx] / 2.)
                        
        if show_plot is True:
            if axes is None:
                f, ax = styles.plot_3panels(x, y, figaxis3=axes, 
                                        ytop=ytop, centered=centered,
                                        styles3=styles3, init=initial,
                                        lengths=lengths_microns,
                                        show_line_at_1=show_line_at_initial)

                return f, ax, x, y
            else: 
                styles.plot_3panels(x, y, figaxis3=axes, 
                                    ytop=ytop, centered=centered,
                                    styles3=styles3, init=initial,
                                    lengths=lengths_microns,
                                    show_line_at_1=show_line_at_initial)

        return x, y