# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 14:34:00 2016

Define attributes and functions related to physical samples

@author: Ferriss
"""
from __future__ import print_function, division, absolute_import
import numpy as np

class Sample():
    """
    Stores basic sample geometry. Most sample info should be stored with 
    the IGSN, the international geosample number, available at geosamples.org.
    You can include the IGSN and/or a description of the sample when 
    setting up the sample as well, e.g.,
    my_sample = Sample(IGSN=IEFERJAI4, description='San Carlos olivine')

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
    
    Fe (total), Fe2+, Fe3_, Mg, Al, and Ti concentratrations can be added 
    in as atoms per formula unit (apfu), and initial water content can be
    added in ppm H2O.
    Sample(Fe2=0.17, initial_water_ppmH2O=5)
    
    You can have it determine the Mg# from the Fe and Mg
    by these using Sample.get_MgNum()
    """
    def __init__(self, IGSN=None, 
                 description=None,
                 thickness_thinslab_microns=[],
                 length_a_microns = [],
                 length_b_microns = [], 
                 length_c_microns = [],                  
                 Fe=0, Fe2=0, Fe3=0., Mg=0, Al=0, Ti=0,
                 initial_water_ppmH2O=None):
        self.IGSN = IGSN
        self.description = description
        self.initial_water_ppmH2O = initial_water_ppmH2O
        self.length_a_microns = length_a_microns
        self.length_b_microns = length_b_microns
        self.length_c_microns = length_c_microns
        self.thicknesses_microns = [length_a_microns, 
                                    length_b_microns,
                                    length_c_microns]
        self.Fe2 = Fe2
        self.Fe3 = Fe3
        self.Mg = Mg
        self.Al = Al
        self.Ti = Ti
        if (Fe2 is not None) and (Fe3 is not None) and (Fe is None):
            self.Fe = Fe2 + Fe3
        else:
            self.Fe = Fe


        def floatify(lengthy):
            """
            Takes a length or list of lengths and makes it into a single 
            floating point number for condensing down lists of thicknesses.
            """
            if lengthy is None:
                floater = None
            elif isinstance(lengthy, list) or isinstance(lengthy, np.ndarray):
                if len(lengthy) > 0:
                    floater = np.mean(lengthy)
                else:
                    floater = None
            else:
                try:
                    floater = float(lengthy)
                except TypeError:
                    floater = lengthy
                    print('length should be integer, float, or list')
            return floater

        twoA = floatify(length_a_microns)
        twoB = floatify(length_b_microns)
        twoC = floatify(length_c_microns)
        thinness = floatify(thickness_thinslab_microns)
        self.thickness_microns = [twoA, twoB, twoC, thinness]
        self.lengths_microns = [twoA, twoB, twoC]

        
    def get_MgNumber(self):
        """ 
        Takes the Fe (or Fe2 + Fe3 if Fe=0) and Mg of a sample assumed to be 
        in atoms per formula unit (apfu) and calculates the Mg#
        """
        if self.Fe == 0:
            self.Fe = self.Fe2 + self.Fe3
        if self.Fe == 0 and self.Mg == 0:
            print("You don't seem to have any Fe or Mg")
            return
        if self.Fe == 0:
            print("You don't seem to have any Fe")
        if self.Mg == 0:
            print("You don't seem to have any Mg")            
        try:
            MgNum = 100. * self.Mg / (self.Fe + self.Mg)
            print('\nThe Mg#:')
            print(MgNum)
            self.MgNum = MgNum
        except TypeError:
            print('\nCheck your Mg and Fe values')