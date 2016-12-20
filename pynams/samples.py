# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 14:34:00 2016

Define attributes and functions related to physical samples

@author: Ferriss
"""
import numpy as np

class Sample():
    """
    Stores basic sample geometry. Most sample info should be stored with 
    the IGSN, the international geosample number, available at geosamples.org.

    For a thin slab where you're not worried about the direction, 
    input thickness under thickness_thinslab_microns, e.g., 
    my_thin_slab = Sample(thickness_thinslab_microns=200.)

    or with multiple measurements, use a list of numbers, e.g.,
    my_thin_slab = Sample(thickness_thinslab_microns=[212., 206., 230.])
    
    For a sample block, measure or guess the orientation and input one or 
    more measurements under length_a/b/c_microns, e.g.,
    my_block = Sample(length_a_microns=1000,
                      length_b_microns=[2012, 2020],
                      length_c_microns=[2410, 2400])

    All thicknesses input should be in microns.
    
    The Sample will automatically generate an attribute thickness_microns that
    is a list of the average thickness || a, b, c, and unoriented.
    e.g., print my_block.thickness_microns produces
    [1000.0, 2016.0, 2405.0, None]
    
    You can add the initial water in ppm H2O if you have a clue what it is.
    
    Fe (total), Fe2+, Fe3_, Mg, Al, and Ti concentratrations can be added 
    in as atoms per formula unit (apfu): Sample(Fe2=0.17)
    
    You can have it determine the Mg# from
    by these usinig Sample.get_MgNum()
    """
    def __init__(self, IGSN=None, 
                 thickness_thinslab_microns=[],
                 length_a_microns = [],
                 length_b_microns = [], 
                 length_c_microns = [],                  
                 Fe=None, Fe2=None, Fe3=None, Mg=None, Al=None, Ti=None,
                 initial_water_ppmH2O=None):
        self.IGSN = IGSN
        self.initial_water_ppmH2O = initial_water_ppmH2O
        self.length_a_microns = length_a_microns
        self.length_b_microns = length_b_microns
        self.length_c_microns = length_c_microns
        self.Fe2 = Fe2
        self.Fe3 = Fe3
        self.Mg = Mg
        self.Al = Al
        self.Ti = Ti
        if (Fe2 is not None) and (Fe3 is not None) and (Fe is None):
            self.Fe = Fe2 + Fe3
            

        if isinstance(length_a_microns, float):
            twoA = length_a_microns
        elif isinstance(length_a_microns, int):
            twoA = float(length_a_microns)
        elif len(length_a_microns) > 0:
            twoA = np.mean(self.length_a_microns)            
        else:
            twoA = None

        if isinstance(length_b_microns, float):
            twoB = length_b_microns
        elif isinstance(length_b_microns, int):
            twoB = float(length_b_microns)
        elif len(length_b_microns) > 0:
            twoB = np.mean(self.length_b_microns)
        else:
            twoB = None        
            
        if isinstance(length_c_microns, float):
            twoC = length_c_microns            
        elif isinstance(length_c_microns, int):
            twoC = float(length_c_microns)
        elif len(length_c_microns) > 0:
            twoC = np.mean(self.length_c_microns)
        else:
            twoC = None
            
        if isinstance(thickness_thinslab_microns, float):
            thinness = thickness_thinslab_microns
        elif isinstance(thickness_thinslab_microns, int):
            thinness = float(thickness_thinslab_microns)
        if len(thickness_thinslab_microns) > 0:
            thinness = np.mean(thickness_thinslab_microns)
        else:
            thinness = None

        self.thickness_microns = [twoA, twoB, twoC, thinness]
        
    def get_MgNumber(self):
        try:
            MgNum = 100. * self.Mg / (self.Fe + self.Mg)
        except TypeError:
            print self.description
            print 'Check Mg and Fe are not None'
        else:
            self.MgNum = MgNum
            return MgNum
