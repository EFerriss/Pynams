# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 14:34:00 2016

Define attributes and functions related to physical samples

@author: Ferriss
"""
import numpy as np
from .uncertainties import ufloat

print '\nIf you see this message and a +/- number, my test worked!'
y = ufloat(3, 0.5)
print y

class Sample():
    """
    Stores basic sample geometry. Most sample info should be stored with 
    the IGSN, the international geosample number, available at geosamples.org.
    """
    def __init__(self, IGSN=None, length_a_microns = [],
                 length_b_microns = [], length_c_microns = [], 
                 length_unoriented_microns=None,
                 mineral_name=None, initial_water=None, ):
        self.IGSN = IGSN
        self.mineral_name = mineral_name
        self.initial_water = initial_water
        self.length_a_microns = length_a_microns
        self.length_b_microns = length_b_microns
        self.length_c_microns = length_c_microns

        if len(length_a_microns) > 0:
            twoA = np.mean(self.length_a_microns)
        else:
            twoA = None
        if len(length_b_microns) > 0:
            twoB = np.mean(self.length_b_microns)
        else:
            twoB = None        
        if len(length_c_microns) > 0:
            twoC = np.mean(self.length_c_microns)
        else:
            twoC = None
            
        self.thickness_microns = [twoA, twoB, twoC, length_unoriented_microns]
