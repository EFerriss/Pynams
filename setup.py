# -*- coding: utf-8 -*-
"""
Created on Wed Nov 04 17:50:54 2015

@author: Ferriss

pynams setup file
"""
from setuptools import setup
setup(
    name = "pynams",
    version = "0.2.0",
    package_dir = {'pynams' : 'pynams'},
    packages = ['pynams', 'pynams.uncertainties']
)
