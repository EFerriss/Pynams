# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 09:24:43 2015

__init__ file for pynams main package

@author: Ferriss
"""
from __future__ import print_function
from .samples import Sample
from .spectra import Spectrum
from .profiles import Profile
from .blocks import Block
from . import styles
from .diffusion import literaturevalues as dlib
from .diffusion.diffusivities import Diffusivities
from . import example_FTIR_spectra
import os
#ftirpath = os.path.dirname(example_FTIR_spectra.__file__)
#example_FTIR_file_location = ''.join((ftirpath,'\\'))
#thisfolder = ''.join((ftirpath,'\\..\\'))