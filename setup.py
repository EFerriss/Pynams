# -*- coding: utf-8 -*-
"""
Created on Wed Nov 04 17:50:54 2015

@author: Ferriss

pynams setup file
"""
from setuptools import setup


setup(
    description = "python for water in nominally anhydrous minerals",
    author = 'Elizabeth Ferriss',
    author_email = 'elizabeth.ferriss@gmail.com',
    url = 'https://github.com/EFerriss/Pynams',
    version = "0.2.4",
    packages = ['pynams', 'pynams.diffusion', 'pynams.example_FTIR_spectra'],
    name = "pynams",
    license = "MIT",
    install_requires = ['lmfit==0.8.3', 'uncertainties>=3.1.3'],
    include_package_data = True
    )
