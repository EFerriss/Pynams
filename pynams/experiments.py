# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 09:33:13 2015

@author: Ferriss

Functions to help with performing and interpreting experiments on 
nominally anhydrous minerals. 


"""
from __future__ import print_function, division, absolute_import
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from scipy import constants
from scipy.interpolate import interp1d
import numpy as np

GAS_CONSTANT = constants.physical_constants['molar gas constant'][0] # J/molK
FARADAY_CONSTANT = constants.physical_constants['Faraday constant'][0]

#%% water concentration as a function of water fugacity   
def solubility_of_H_in_olivine(Celsius, water_fugacity_GPa=None, 
                               pressure_GPa=1, author='Mosenfelder',
                               printout=True):
    """
    Returns fully saturated solubility of olivine in H/10^6 Si and ppm H2O.

    Required: temperature in Celsius and either the pressure in GPa (default=1)
    or the water fugacity in GPa. Setting the water fugacity overrides the
    water fugacity calculated based on the pressure and an assumed temperature
    of 1100 C.
        
    Options for volume change with pressure dV are set by kwarg author.
    Use 'Mosenfelder' (default) for dV from Mosenfelder et al. 2006, 
    'Kohlstedt' for Kohlstedt et al. 1996, or 'Zhao' for Zhao et al. 2004
    """
    Kelvin = Celsius + 273.15
    
    if water_fugacity_GPa is None and pressure_GPa is None:
        print('Required: either water_fugacity in GPa or pressure in GPa')
        return
    
    if author == 'Mosenfelder':
        # Mosenfelder et al. 2006
        dV = 10.2e-6 # m^3/mol
    elif author == 'Kohlstedt':
        # Kohlstedt et al. 1996
        dV = 10.0e-6 # m^3/mol
    elif author == 'Zhao':
        # Zhao et al. 2004
        dV = 10.6e-6 # m^3/mol
    else:
        print('author options are Mosenfelder, Kohlstedt, or Zhao')
        print('defaulting to Mosenfelder et al. 2006')
        dV = 10.2e-6 # m^3/mol

    # water fugacity as a function of pressure at constant temp of 1100 C
    # www.esci.umn.edu/people/researchers/withe012/fugacity.htm
    GPa = [0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    fugacity = [0.00997, 0.098, 0.6, 2, 12.13, 55.3, 218.66, 784.98, 2625, 
                8305.9, 25092.63, 72904., 204764.6, 558205.715, 1481665.06]

    if pressure_GPa is None:
        f_spline_fit = interp1d(np.log10(fugacity), GPa, kind='cubic')
        pressure_GPa = f_spline_fit(np.log10(water_fugacity_GPa))

    if water_fugacity_GPa is None:
        f_spline_fit_p2f = interp1d(GPa, fugacity, kind='cubic')
        water_fugacity_GPa = f_spline_fit_p2f(pressure_GPa)

    # The exponent part keeps going to 1 because all the numbers are so small.
    pressure_Nm2 = pressure_GPa * 1e9
    exponent_part = np.exp((-1.0*pressure_Nm2*dV) / (GAS_CONSTANT*Kelvin))
    A = 2.45 # H/10^6 Si / GPa
    n = 1
    C = A * (water_fugacity_GPa**n) * exponent_part
    ppm = C * 60
    if printout is True:
        print('Solubility of H in olivine:', ppm, 'ppm H2O')
    return ppm


def convertHunit(conc, from_unit='H/10^6 Si', to_unit='ppm H2O', phase='Fo90',
             printout=True):
    """
    Convert hydrogen concentrations to/from H/10^6 Si and ppm H2O.
    Based on Table 3 of Denis et al. 2013
    """
    if phase == 'Fo90':
        H_to_1_ppm = 16.35
    elif phase == 'opx':
        H_to_1_ppm = 11.49
    elif phase == 'cpx':
        H_to_1_ppm = 11.61
    else:
        print('Valid options for phase are Fo90, opx, and cpx')
        return
      
    if from_unit == 'H/10^6 Si':
        if to_unit == 'ppm H2O':
            new_conc = conc / H_to_1_ppm
        elif to_unit == 'per m3':
            new_conc = conc * (1.0/308.67) * (1e30)
        else:
            print('only going to units "ppm H2O" and "per m3"')
            return
        
    elif from_unit == 'ppm H2O':
        if to_unit == 'H/10^6 Si':
            new_conc = conc * H_to_1_ppm
        elif to_unit == 'per m3':
            new_conc = (conc * H_to_1_ppm) * (1.0/308.67) * (1e30)
        else:
            print('only going to "H/10^6 Si" or "per m3"')
            return
            
    elif from_unit == 'per m3':
        if to_unit == 'H/10^6 Si':
            new_conc = conc / ((1.0/308.67) * (1e30))
        elif to_unit == 'ppm H2O':
            new_conc = (conc / ((1.0/308.67) * (1e30))) / H_to_1_ppm
        else:
            print('only going to "H/10^6 Si" or "ppm H2O"')
            return
        
    else:
        print('Only going from H/10^6 Si, ppm H2O, and per m3 for now')
        return
        
    if printout is True:
        output = ' '.join(('{:.2f}'.format(conc), from_unit, '=', 
                           '{:.2f}'.format(new_conc), to_unit, 'for', phase))
        print(output)
    return new_conc


def bubble_tower(panel='middle', minor_setting=40, 
                 major_setting=155., major_gas='CO2',
                 target_log10fO2=-14., minor_gas='CO',
                 furnace_diameter_inches=2.25,
                 max_flow_rate=4.0, figsize=(6,5)):
    """
    For use estimating gas mixing ratios and settings for gas mixing 
    furnace in Dave Walker's lab. 
    
    Required input: 
        panel (left or middle - default),
        minor gas (usually CO) minor_setting (default=40)
        major gas (usually CO2) major_setting (default=155)
        
    Flow rate measurements are included in dictionaries.
    The dictionary keys are the needle valve settings, and
    the values are lists of flow rate measurements in mL / seconds
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Needle valve setting')
    ax.set_ylabel('Flow rate mL/s')
    ax.set_ylim(0, max_flow_rate)

    monoxide = dict()
    dioxide = dict()

    if panel == 'left':
        ax.set_xlim(0, 175)

        if major_gas == 'CO2':
            monoxide[30] = [0.]
            monoxide[40] = [0.5/13.5, 1./29.1]
            monoxide[50] = [1./12.6]
            monoxide[60] = [3./23.0]
            monoxide[135] = [50./29.6]
            monoxide[150] = [50./25.2]
            monoxide[168] = [50./19.7]
                    
            dioxide[100] = [50./78.4]
            dioxide[125] = [50./28.3]
            dioxide[150] = [50./17.9]
            dioxide[170] = [50./14.1]
        elif major_gas == 'CO':
            monoxide[100] = [10./15.2]
            monoxide[145] = [50./27.6]
            monoxide[173] = [50./18.4]

            dioxide[60] = [5./17.2]
            dioxide[70] = [10./21.8]
            dioxide[80] = [20./29.0]
            dioxide[90] = [20./21.3]
            dioxide[100] = [30./23.8]
            dioxide[115] = [30./18.0]
            dioxide[116] = [40./32.1, 30./17.7]
            dioxide[117] = [30./17.1, 30./17.7]
            dioxide[119] = [50./29.4]
            dioxide[168] = [50./14.5]
            
        else:
            print('major_gas either CO2 (default) or CO')
            return

    elif panel == 'middle':
        ax.set_xlim(0, 150)

        if major_gas == 'CO2':        
            monoxide[7] = [1./13.2] # jerky
            monoxide[19] = [2./15.3]        
            monoxide[30] = [5./25.1]
            monoxide[40] = [5./18.8]
            
            dioxide[70] = [20./26.0, 20./25.6]
            dioxide[100] = [40./27.55]
            dioxide[145] = [50./18.4]
            dioxide[150] = [50./14.4]
        else:
            print('Only CO2 as major gas on middle panel')
            return
    else:
        print('Valid entries for "panel" are left and middle')
        return
            
    # Change from list of flow rates to a single mean value
    dioxide.update((x, np.mean(y)) for x, y in dioxide.items())
    monoxide.update((x, np.mean(y)) for x, y in monoxide.items())
    
    p_dioxide = np.polyfit(dioxide.keys(), dioxide.values(), 1)
    p_monoxide = np.polyfit(monoxide.keys(), monoxide.values(), 1)

    # plot flow rate measurements
    ax.plot(dioxide.keys(), dioxide.values(), 'o', markerfacecolor='w',
            markeredgecolor='b', markeredgewidth=2, label='CO$_2$')
    ax.plot(monoxide.keys(), monoxide.values(), '^', markerfacecolor='w',
            markeredgecolor='r', markeredgewidth=2, label='CO')
    x = ax.get_xlim()
    ax.plot(x, np.polyval(p_dioxide, x), '-', color='b')
    ax.plot(x, np.polyval(p_monoxide, x), '-', color='r')
    
    # minimum acceptable flow rate for good mixing 
    # ~ 1/10 cross-sectional area of furance
    radius_cm = (furnace_diameter_inches / 2.) * 2.54
    cross_section_area_cm2 = constants.pi * radius_cm * radius_cm
    minflow = cross_section_area_cm2 / 10.
    ax.plot(x, [minflow, minflow], '-g')
    ax.text(x[0]+5, minflow+0.05, 'minimum flow rate for good mixing')
    ax.legend(loc=2)
    
    if major_gas == 'CO2':
        CO2_setting = major_setting
        CO_setting = minor_setting       
    else:
        CO2_setting = minor_setting
        CO_setting = major_setting
        minor_gas = 'CO2'

    yCO2 = np.polyval(p_dioxide, CO2_setting)
    yCO = np.polyval(p_monoxide, CO_setting)

    if major_gas == 'CO2':
        major_flowrate = yCO2
        minor_flowrate = yCO
    else:
        minor_gas = 'CO2'
        major_flowrate = yCO
        minor_flowrate = yCO2

    # Arrows on data
    ax.text(CO2_setting, yCO2, '$\leftarrow*$', fontsize=18, 
            rotation=90, ha='center', va='bottom')
    ax.text(CO_setting, yCO, '$\leftarrow*$', fontsize=18, 
            rotation=90, ha='center', va='bottom')
    
    ax.text(10., 2.0, ''.join(('major gas ', major_gas, '\nsetting ', 
                               '{:.0f}'.format(major_setting), '\n',
                               '{:.2f}'.format(major_flowrate), ' mL/s')))
    
    ax.text(10., 1.3, ''.join(('minor gas ', minor_gas, '\nsetting ', 
                               '{:.0f}'.format(minor_setting), '\n',
                               '{:.2f}'.format(minor_flowrate), ' mL/s')))

    ax.set_title(' '.join((panel, 'panel')))
    percentCO2 = 100. * yCO2 / (yCO2 + yCO)
    ax.text(10., 0.9, ''.join(('{:.2f}'.format(percentCO2), '% CO$_2$')))


def fO2(celsius, bars=1., buffer_curve='QFM'):
    """ 
    Input:
        Temperature in Celcius
        Pressure in bars (default is 1)
        Buffer curve. Options are QFM (default), NNO, and IW
            QFM = quartz - fayalite - magnetite buffer (default)
            NNO = Ni-NiO
            IW = iron-wustite; Fe-FeO

    Output is log10 oxygen fugacity in bars for buffer_curve

    Regression of data from O’Neill (1987b) from Herd, 2008.
    """
    Kelvin = celsius + 273.15

    if buffer_curve == 'QFM':    
        A = -24935.0
        B = 8.489
        C = 0.0
    elif buffer_curve == 'NNO':
        A = -24525.4
        B = 8.944
        C = 0.0
    elif buffer_curve == 'IW':
        A = -27489.
        B = 6.702
        C = 0.055
    else: 
        print('Only QFM, IW, and NNO supported for now')
        return False
        
    logfO2 = ((A / Kelvin) + B + (C * (bars-1.0) / Kelvin))
    return logfO2
    

def furnace_calibration(reading_celsius):
    """Calibration for furnace 4 in Dave Walker's lab at LDEO
    based on gold calibration"""
    real_temp = reading_celsius - (6. / (reading_celsius/1064.))
    print('{:.1f}'.format(real_temp), 'degrees C')
    return real_temp


def log10fO2_from_V(volts, celsius, buffermix='CO-CO2'):
    """
    Input:
        Volts reading from pO2 meter 
        Temperature in Celsius 
        
        
    Returns log10 oxygen fugacity in bars 
    
    Uses the Nernst equation:
    E(Volts) = -(RT/zF)lnQ 
    and for now assumes a mixture of CO and CO2, so Q ~= fO2 ^(1/2) and 
    number of electrons z = 2
    """
    Kelvin = celsius + 273.15
    z = 2.
    exponent_in_Q = -0.5
    
    my_constant = z * FARADAY_CONSTANT / (exponent_in_Q * GAS_CONSTANT * 2.303) 
    logfO2 = -1. * my_constant * volts / (Kelvin)
    
    print('{:.1f}'.format(logfO2), 'log10 oxygen fugacity')
    return logfO2

    
def V_from_log10fO2(log10fO2, celsius, buffermix='CO-CO2'):
    """
    Reverse of log10fO2_from_V
    
    Input:
        Target log10fO2 for furnace 
        Temperature in Celsius 
        
        
    Returns target V reading on pO2 sensor using the Nernst equation:
    E(Volts) = -(RT/zF)lnQ 
    and for now assuming a mixture of CO and CO2, so Q ~= fO2 ^(1/2) and 
    number of electrons z = 2
    """
    Kelvin = celsius + 273.15
    z = 2.
    exponent_in_Q = -0.5    
    my_constant = z * FARADAY_CONSTANT / (exponent_in_Q * GAS_CONSTANT * 2.303) 
    Volts = -1. * log10fO2 * Kelvin / my_constant
    print('\n', '{:.3f}'.format(Volts), 'target Volts on pO2 meter\n')

    
def make_capsule_shape(x, y, height, outerD, innerD, shape='regular'):
    """
    Makes and returns path representing the capsule shape for use in 
    pressure vessel schematics as in pressure_design()
    """
    thick = (outerD - innerD) / 2.
    if shape == 'regular':
        verts = [(x + thick*2 + innerD, y),
                 (x, y),
                 (x, y + height),
                 (x + thick, y + height),
                 (x + thick, y + innerD/2.),
                 (x + thick + innerD/2., y + thick),
                 (x + thick + innerD, y + innerD/2.),
                 (x + thick + innerD, y + height),
                 (x + thick*2 + innerD, y + height),
                 (0., 0.)]
        codes = [Path.MOVETO] + ([Path.LINETO] * 8) + [Path.CLOSEPOLY]
    elif shape == 'suaged':
        th_flap = thick/2.
        verts = [(x + thick*2 + innerD, y),
                 (x, y),
                 (x, y + height + th_flap),
                 (x + thick + innerD/2., y + height + th_flap),
                 (x + thick + innerD/2., y + height),
                 (x + thick, y + height),
#                 (x + thick, y + thick),
                 (x + thick, y + innerD/2.),
                 (x + thick + innerD/2., y + thick),
                 (x + thick + innerD, y + innerD/2.), 
#                 (x + thick + innerD, y + thick),
                 (x + thick + innerD, y + height),
                 (x + thick + innerD/2., y + height),
                 (x + thick + innerD/2., y + height + th_flap),
                 (x + thick*2 + innerD, y + height + th_flap),
                 (x + thick*2 + innerD, y + height - th_flap),
                 (0., 0.)]
        codes = [Path.MOVETO] + ([Path.LINETO] * (len(verts)-2)) + [Path.CLOSEPOLY]
    path = Path(verts, codes)
    return path

style_pressure_medium = {'hatch' : '..', 'facecolor' : 'lightgrey'}
style_graphite = {'hatch' : '/', 'facecolor' : 'dimgrey'}
style_MgO = {'hatch' : None, 'facecolor' : 'white', 'edgecolor' : 'k'}
style_pyrophyllite = {'hatch' : 'xx', 'facecolor' : 'hotpink'}
style_capsule = {'facecolor' : 'orange', 'edgecolor' : 'k'}
style_buffer = {'facecolor' : 'w', 'hatch' : 'xxxxx', 'edgecolor' : 'k'}

       
def pressure_design(capsule_material = 'copper',
                    pressure_medium_material='BaCO$_3$',
                    sleeve_material='pyrophyllite',
                    buffer_material='H$_2$O, Ni,\nNiO, SC ol.,\nSC enst.',
                    lid_shape = 'bevel', # or flat or suaged
                    h_graphite_button=1.5, 
                    h_pressure_medium=33.35,
                    h_graphite_cylinder=33.35,
                    h_sleeve=11.5,
                    h_sleeve_bottom=1.,
                    h_capsule = 9.7,
                    h_lid=1., # total 
                    h_MgO_base=10.6,
                    h_MgO_wafer=1.5,
                    h_MgO_top=10.6,
                    od_pressure_medium=18.9,
                    od_graphite_cylinder=11.5,
                    id_graphite_cylinder=10.0,
                    id_sleeve=8.7,
                    id_capsule=6.7,
                    legend_on=True, 
                    figsize=(3., 3),):
    """Creates and returns figure and axis
    showing experimental setup for a high pressure 
    experiment used to hydrate olivine and potentially
    other nominally anhydrous minerals in the green 4-post piston cylinder 
    at Lamont. All dimensions are input in mm. h=height, d=diameter,
    od=outer diameter, id=inner diameter
    """
    # reset all style labels to None for legend 
    style_capsule['color'] = None
    style_buffer['label'] = buffer_material
    for style in [style_graphite, style_MgO, style_capsule]:
        style['label'] = None
    
    d_graphite_button = od_pressure_medium
    od_MgO_base = id_graphite_cylinder
    od_sleeve = id_graphite_cylinder

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_xlim(0, od_pressure_medium)
    ax.set_xlabel('(mm)')
    ax.set_ylabel('(mm)')
    h_guts = h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top
    highest_point = max(h_pressure_medium, h_graphite_cylinder, h_guts)
    ax.set_ylim(0., h_graphite_button + highest_point + 2.)
    plt.tick_params(axis='x', top='off')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    th_gc = (od_graphite_cylinder - id_graphite_cylinder) / 2.
    xgc = (od_pressure_medium - od_graphite_cylinder) / 2.

    style_pressure_medium['label'] = pressure_medium_material
    pressure_medium = patches.Rectangle((0., h_graphite_button), # (x,y)
                                        od_pressure_medium, # width
                                        h_pressure_medium, # height
                                        **style_pressure_medium)

    graphite_button = patches.Rectangle((0., 0.), d_graphite_button,
                                        h_graphite_button, **style_graphite)

    style_graphite['label'] = 'graphite'
    graphite_cylinder = patches.Rectangle((xgc, h_graphite_button), 
                             od_graphite_cylinder, h_graphite_cylinder, 
                             **style_graphite)

    the_guts = patches.Rectangle((xgc + th_gc, h_graphite_button), 
                                 id_graphite_cylinder,
                                 h_graphite_cylinder, facecolor='w')
                             
    MgO_base = patches.Rectangle((xgc+th_gc, h_graphite_button),
                                 od_MgO_base, h_MgO_base, **style_MgO)
   
    # sleeve around capsule
    sleeve_path = make_capsule_shape(x=xgc + th_gc, 
                                     y=h_graphite_button + h_MgO_base,
                                     height=h_sleeve, 
                                     outerD=od_sleeve, 
                                     innerD=id_sleeve)
    if sleeve_material == 'MgO':
        style_sleeve = style_MgO.copy()
    elif sleeve_material == 'pyrophyllite':
        style_sleeve = style_pyrophyllite.copy()
        style_sleeve['label'] = 'pyrophyllite'
    else: 
        print('unknown sleeve material. Assuming pyrophyllite')
        style_sleeve = style_pyrophyllite.copy()
        style_sleeve['label'] = 'pyrophyllite'

    sleeve = patches.PathPatch(sleeve_path, **style_sleeve)

    # capsule
    th_sleeve = (od_sleeve - id_sleeve) / 2.
    if (lid_shape == 'bevel') or (lid_shape == 'flat'):
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule)
    elif lid_shape == 'suaged':
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule,
                       shape='suaged')  
    else:
        print('valid entries for lid_shape are flat, bevel, and suaged')
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule)                          
                       
    if capsule_material == 'copper':
        style_capsule['label'] = 'copper'
        style_capsule['facecolor'] = 'orange'
    elif capsule_material == 'silver':
        style_capsule['label'] = 'silver'
        style_capsule['facecolor'] = 'silver'
    elif capsule_material == 'gold':
        style_capsule['label'] = 'gold'
        style_capsule['facecolor'] = 'gold'
    elif capsule_material == 'platinum':
        style_capsule['label'] = 'platinum'
        style_capsule['facecolor'] = 'lightblue'
    elif capsule_material == 'nickel':
        style_capsule['label'] = 'nickel'
        style_capsule['facecolor'] = 'lightsage'
    else:
        print('unknown capsule material')
        style_capsule['label'] = 'capsule'        
    capsule = patches.PathPatch(capsule_path, **style_capsule)

    # MgO on top
    MgO_wafer = patches.Rectangle((xgc + th_gc, 
                                   h_graphite_button + h_MgO_base + h_sleeve),
                                   od_MgO_base, h_MgO_wafer, **style_MgO)

    style_MgO['label'] = 'MgO'
    MgO_top = patches.Rectangle((xgc + th_gc, 
                    h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer),
                    od_MgO_base, h_MgO_top, **style_MgO)

    thermocouple = patches.Rectangle((od_pressure_medium/2.-0.5, 
                    h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer),
                    1., h_MgO_top, facecolor='w')
    ax.plot([od_pressure_medium/2-0.15, od_pressure_medium/2.-0.15],
            [h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top + 2.,
             h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer], 
             color='r', linewidth=1)
    ax.plot([od_pressure_medium/2+0.15, od_pressure_medium/2.+0.15],
            [h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top + 2.,
             h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer], 
             color='b', linewidth=1)

    # buffer
    th_capsule = (id_sleeve - id_capsule) / 2.
    buffer_inside = patches.Rectangle((xgc + th_gc + th_sleeve + th_capsule,
                    h_graphite_button + h_MgO_base + h_sleeve_bottom + th_capsule),
                    id_capsule, h_capsule - th_capsule, **style_buffer)

    # capsule lid
    del style_capsule['label'] # so it doesn't appear twice in the legend        
    if lid_shape == 'flat':
        lid = patches.Rectangle((xgc + th_gc + th_sleeve, 
                    h_graphite_button + h_MgO_base + h_sleeve_bottom + h_capsule),
                    id_sleeve, h_lid, **style_capsule)
    elif lid_shape == 'bevel':
        x = xgc + th_gc + th_sleeve
        y = h_graphite_button + h_MgO_base + h_sleeve_bottom + h_capsule
        th_lid = h_lid / 2.
        th_capsule = (id_sleeve - id_capsule) / 2.
        lid_verts = [(x + th_capsule, y), 
                     (x, y),
                     (x, y + th_lid), 
                     (x + id_sleeve, y + th_lid), 
                     (x + id_sleeve, y), 
                     (x + id_sleeve - th_capsule, y), 
                     (x + id_sleeve - th_capsule, y - th_lid), 
                     (x + th_capsule, y - th_lid), 
                     (0., 0.)]
        lid_codes = [Path.MOVETO] + ([Path.LINETO] * 7) + [Path.CLOSEPOLY]
        lid_path = Path(lid_verts, lid_codes)
        lid = patches.PathPatch(lid_path, **style_capsule)
    elif lid_shape == 'suaged':
        th_flap = th_capsule / 2.
        x = xgc + th_gc + th_sleeve + th_flap
        ystart = h_graphite_button + h_MgO_base + h_sleeve_bottom
        y = ystart + h_capsule - th_flap*2
        th_lid = h_lid / 2.
        th_capsule = (id_sleeve - id_capsule) / 2.
        lid_verts = [(x + th_capsule - th_flap, y), 
                     (x, y),
                     (x, y + th_lid), 
                     (x + id_sleeve - th_capsule, y + th_lid), 
                     (x + id_sleeve - th_capsule, y), 
                     (x + id_sleeve - th_capsule - th_flap, y), 
                     (x + id_sleeve - th_capsule - th_flap, y - th_lid), 
                     (x + th_capsule - th_flap, y - th_lid), 
                     (0., 0.)]
        lid_codes = [Path.MOVETO] + ([Path.LINETO] * 7) + [Path.CLOSEPOLY]
        lid_path = Path(lid_verts, lid_codes)
        lid = patches.PathPatch(lid_path, **style_capsule)
    else:
        print('valid entries for lid_shape are flat, bevel, and suaged')
        lid = patches.Rectangle((xgc + th_gc + th_sleeve, 
                    h_graphite_button + h_MgO_base + th_sleeve + h_capsule),
                    id_sleeve, h_lid, **style_capsule)

    ax.add_patch(pressure_medium)
    ax.add_patch(graphite_button)
    ax.add_patch(graphite_cylinder)
    ax.add_patch(the_guts)
    ax.add_patch(MgO_base)
    ax.add_patch(sleeve)
    ax.add_patch(buffer_inside)
    ax.add_patch(MgO_wafer)
    ax.add_patch(MgO_top)
    ax.add_patch(thermocouple)
    ax.add_patch(capsule)
    ax.add_patch(lid)
      
    fig.tight_layout()
    if legend_on is True:
        plt.subplots_adjust(right=0.55, left=0.17, bottom=0.15, top=0.9)
        ax.legend(bbox_to_anchor=(2.25, 0.8), frameon=False)
    else:
        plt.subplots_adjust(right=0.9, left=0.17, bottom=0.15, top=0.9)
    return fig, ax
    
