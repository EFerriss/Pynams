"""
Code for transforming FTIR peak areas into water concentrations and 
handling individual peaks. 

The real meat of the code is over in spectra.py and profiles.py, 
but they both use functions stored here.

"""
from __future__ import print_function, division, absolute_import
from .uncertainties import ufloat
import numpy as np

def absorption_coefficients(phase, calibration):
    """
    Input the phase (olivine or cpx) and calibration reference
    (Bell or Withers), and it returns the absorption coefficient.
    This is so you don't have to look it up over and over again.
    """
    if (calibration == 'Bell') and (phase=='cpx'):
        absorption_coeff = 1.0 / ufloat(7.09, 0.32)

    # Bell et al. 2003        
    elif (calibration == 'Bell') and (phase=='olivine'):
        absorption_coeff = ufloat(0.188, 0.012)
        
    # Withers et al. 2012
    elif (calibration == 'Withers') and (phase=='olivine'):
        absorption_coeff = ufloat(0.119, 0.006)
        
    else:
        print('Calibrations supported so far:')
        print('   Bell et al. 1995 for cpx')
        print('   Bell et al. 2003 for olivine')
        print('   Withers et al. 2012 for olivine')
        return

    return absorption_coeff        
    
def area2water(area_cm2, phase='cpx', calibration='Bell', peak_idx=None):
    """
    Takes the peak area in cm-2, multiplies by the absorption coefficient, 
    and return the water concentration in ppm H2O
    """
    absorption_coeff = absorption_coefficients(phase=phase, 
                                               calibration=calibration,
                                               peak_idx=peak_idx)
    w = absorption_coeff * area_cm2
    return w
        
def make_gaussian(pos, h, w, x=np.linspace(3000, 4000, 150)):
    """
    Make and return Gaussian curve over specified range. This is used
    for peak fitting.
    """
    y = h * np.e**(-((x-pos) / (0.6005615*w))**2)
    return y

