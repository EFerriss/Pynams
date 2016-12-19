# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 14:34:00 2016

Define attributes and functions related to physical samples

@author: Ferriss
"""
import numpy as np

#%% 
class Sample():    
    def __init__(self, thickness_microns=None, IGSN=None, 
                 mineral_name=None, initial_water=None, length_a_microns = [],
                 length_b_microns = [], length_c_microns = []):
        self.thickness_microns = thickness_microns
        self.IGSN = IGSN
        self.mineral_name = mineral_name
        self.initial_water = initial_water
        self.length_a_microns = length_a_microns
        self.length_b_microns = length_b_microns
        self.length_c_microns = length_c_microns
        
    def get_thicknesses(self):
        twoA = np.mean(self.length_a_microns)
        twoB = np.mean(self.length_b_microns)
        twoC = np.mean(self.length_c_microns)
        self.thickness_microns = [twoA, twoB, twoC]
        return [twoA, twoB, twoC]