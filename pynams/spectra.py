"""
Code for processing and plotting FTIR spectra.

The central concept is an object class called Spectrum, which contains
relevant attributes and methods.

A few defaults are set up that are geared toward H in nominally anhydrous
minerals (NAMs) such a plotting wavenumber range from 3000 to 4000 /cm and
creating baselines between 3200 and 3700 /cm. 

Copyright (c) 2015 Elizabeth Ferriss

This software was developed using Python 2.7.

"""
from __future__ import print_function, division, absolute_import
from future.builtins import range

import numpy as np
import os

#import styles
#import diffusion
#import gc
#import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
#from .uncertainties import ufloat
#from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
#import matplotlib.gridspec as gridspec
#from matplotlib.ticker import MultipleLocator
#import string as string
#import lmfit
#import uncertainties
#from matplotlib.backends.backend_pdf import PdfPages
#import xlsxwriter
#import json
#from scipy import signal as scipysignal
#import scipy.interpolate as interp


class Spectrum():
    """
    For plotting and manipulating Fourier transform infrared (FTIR) spectra.
    
    Requires an fname, which is the first part of the FTIR filename, e.g.,
    if your raw FTIR spectrum is saved in C:my_folder as my_file.CSV,
    the fname = 'my_file', the filetype='.CSV', and the folder should be 
    set to 'C:my_folder'. The fname is critical both for locating your file
    and importing your data and it is used as the default heading and label
    for a lot of the plotting functions, so it is required. 
    
    The simplest call would look like this:
    my_spectrum = Spectrum('test_spectrum')
    
    The filetype defaults to .CSV, and that's the filetype this code handles
    best, although it's done ok with .txt. 
    
    The folder defaults to the present working directory, but you'll 
    probably need to set the folder explicitly using a string like this one:
    'C:\\Users\\Ferriss\\Documents\\Code\\olivine\\'
    
    Once you make your Spectrum, it will automatically open the file and read 
    in the wavenumbers and absorbances. There is a test file test.CSV in 
    the pynams folder.
    
    The sample thickness can be input in the following ways: 
        1. By specifying a sample, a pynams Sample object, 
           that has thickness information input as thickness_thinslab_microns.
        2. By specifying a sample, a pynams Sample object, that has 
           length information in at least one direction AND specifying the 
           raypath as 'a', 'b', or 'c'. The thickness is then set the as 
           the length of the Sample parallel to the raypath.
        3. You can set thickness_microns directly
        4. You can have it guess the thickness based on the Si-O overtones
           by first making the Spectrum object and then calling the 
           thickness_from_SiO() method.
    
    The baseline is by default a line between wavenumbers 3200 cm-1 
    (base_low_wn) and 3700 (base_high_wn), and a quadratic baseline would
    deviate from that line relative to base_mid_wn in the middle. Those 
    values can be most easily set under make_baseline for an individual
    spectra, but for handling whole profiles at a time, it's convenient
    to be able to set these attributes and then be able to refer to them
    later on.
    """
    def __init__(self, fname, filetype='.CSV', folder='', sample=None,
                 thickness_microns=None, raypath=None, polar=None, 
                 base_low_wn=3200, base_high_wn=3700,
                 base_mid_wn=3550,
                 ):
        self.fname = fname
        self.folder = folder
        self.sample = sample
        self.thickness_microns = thickness_microns
        self.base_high_wn = base_high_wn
        self.base_low_wn = base_low_wn
        self.base_mid_wn = base_mid_wn
        self.polar = polar 
        self.filetype = filetype
        self.raypath = raypath        
        self.filename = self.folder + self.fname + self.filetype

        ### Get the wavenumber and absorbance data from the FTIR file
        if os.path.isfile(self.filename):
            # read in the signal
            if self.filetype == '.CSV':
                signal = np.loadtxt(self.filename, delimiter=',')
            elif self.filetype == '.txt':
                try:
                    signal = np.loadtxt(self.filename, delimiter='\t', 
                                        dtype=None) 
                except ValueError:
                    print('\nProblem reading this file format. Try .CSV')
            else:
                print('For now only CSV and txt files work.')
            # sort signal by wavenumber
            signal = signal[signal[:,0].argsort()]
            self.wn_full = signal[:, 0]
            self.abs_raw = signal[:, 1]    
            print('Yay. The test worked.')
            print('Try printing your_spectrum.wn_full for wavenumbers')
            print('and your_spectrum.abs_raw for the absorbances.')
        else:
            print('There is a problem finding the file.')
            print('filename =', self.filename)
            print('current working directory = ', os.getcwd())
            print('Maybe check the folder name')

#    # full range of measured wavenumber and absorbances
#    wn_full = None
#    abs_raw = None
#    abs_full_cm = None
#
#    # other metadata with a few defaults
#    position_microns_a = None
#    position_microns_b = None
#    position_microns_c = None
#    instrument = None
#    spot_size_x_microns = 100
#    spot_size_y_microns = 100
#    resolution_cm = 4
#    nscans = 100
#    date = None
#    other_name = None
#    # default baseline range
#    base_wn = None
#    base_abs = None
#    # need these for quadratic baseline
#    base_mid_yshift = 0.04
#    base_w_small = 0.02
#    base_w_large = 0.02
#    # calculated after baseline set in subtract_baseline
#    abs_nobase_cm = None
#    area = None
#    # peak fit information
#    peakpos = None
#    numPeaks = None
#    peak_heights = None
#    peak_widths = None
#    peak_areas = None
#    peakshape = None        
#    # Experimental information if applicable
#    temperature_celsius = None
#    time_seconds = None
#    D_area_wb = None
#   
#    def set_thick(self):
#        """Set spectrum thickness_microns from raypath and sample.thickness_microns.
#        Similar utility exists for Profiles"""
#        if self.sample is None and self.thickness_microns is None:
#            print('Estimating sample thickness from SiO overtones')
#            self.thickness_microns = self.thickness_from_SiO()
#        else:
#            s = self.sample
#        
#        if s.thickness_microns is None:
#            s.thickness_microns = get_3thick(s)
#
#        if isinstance(self.sample.thickness_microns, float) is True:
#            th = s.thickness_microns
#        elif isinstance(self.sample.thickness_microns, list) is True:
#            if len(s.thickness_microns) == 3:
#                if self.raypath == 'a':
#                    th = s.thickness_microns[0]
#                elif self.raypath == 'b':
#                    th = s.thickness_microns[1]
#                elif self.raypath == 'c':
#                    th = s.thickness_microns[2]
#                else:
#                    print('Need spectrum raypath a, b, or c')
#                    return False
#            else:
#                print('sample.thickness_microns must have 1 or 3 values')
#                return False
#        self.thickness_microns = th
#        return th
#
#    def thickness_from_SiO(self, show_plot=False, printout=False,
#                           size_inches=(3, 2.5)):
#        """Estimates the sample thickness based on area of the Si-O overtones
#           Using Eq 1 of Matveev and Stachel 2007"""
#        if self.wn_full is None:
#            check = self.get_data()
#            if check is False:
#                return False
#
#        # make temporary baseline under Si-O overtone peaks
#        if self.base_abs is not None:
#            self.save_baseline(baseline_ending='baseline-temp.CSV')
#            retrieveBaseline = True
#        else:
#            retrieveBaseline = False           
#
#        self.make_baseline(wn_low=1625, wn_high=2150, show_plot=False)
#        SiOarea = self.area_under_curve(show_plot=False, printout=printout)
#        thickness_microns = SiOarea / 0.6366
#
#        if show_plot is True:    
#            fig, ax = self.plot_showbaseline(size_inches=(6, 6), #pad_top=10.,
#                                             wn_xlim_left=2200., 
#                                             wn_xlim_right=1200)
#                                                         
#            ax.set_title(''.join((self.fname, '\n',
#                                    '{:.0f}'.format(thickness_microns),
#                                    ' $\mu$m thick')))
#
#        if retrieveBaseline is True:
#            self.get_baseline(baseline_ending='baseline-temp.CSV')
#        else:
#            self.base_abs = None
#
#        # This has to come after the plotting or else the baseline will be offset
#        self.thickness_microns = thickness_microns
#
#        if show_plot is True:
#            return fig, ax
#        return thickness_microns
#
#    def find_lowest_wn_over_given_range(self, wn_mid_range_high=3500., 
#                                        wn_mid_range_low=3300.,
#                                        relative=True):
#        """Take a spectrum and wavenumber range (default 3300-3500 cm-1)
#        and returns the wavenumber with the lowest absorbance within that range."""
#        if self.abs_full_cm is None:
#            self.start_at_zero()
#    
#        ## find minimum relative to a linear baseline
#        if relative is True:            
#            self.make_baseline(linetype='line', show_plot=False)
#            idx_mid_high = (np.abs(self.base_wn-wn_mid_range_high)).argmin()
#            idx_mid_low = (np.abs(self.base_wn-wn_mid_range_low)).argmin()
#            if idx_mid_high == idx_mid_low:
#                print('basenumber range not established. Check wavenumber range')
#                return False
#        
#            abs_nobase = self.subtract_baseline()
#            mid_abs_range = abs_nobase[idx_mid_low:idx_mid_high]
#            
#            mid_wn_range = self.base_wn[idx_mid_low:idx_mid_high]
#            idx_abs_mid = mid_abs_range.argmin()
#            WN_MID = mid_wn_range[idx_abs_mid]
#
#        else:        
#            ## finds absolute minimum over range
#            idx_mid_high = (np.abs(self.wn_full-wn_mid_range_high)).argmin()
#            idx_mid_low = (np.abs(self.wn_full-wn_mid_range_low)).argmin()
#            mid_abs_range = self.abs_full_cm[idx_mid_low:idx_mid_high]
#            mid_wn_range = self.wn_full[idx_mid_low:idx_mid_high]
#            idx_abs_mid = mid_abs_range.argmin()
#            WN_MID = mid_wn_range[idx_abs_mid]
#
#        return WN_MID
#    
#
#    def find_peaks(self, widths=np.arange(1,40), 
#                   linetype='line', shiftline=0.):
#        """Locate wavenumbers of npeaks number of most prominent peaks
#        in wavenumber range of interest wn_high to wn_low"""
#        if self.abs_nobase_cm is None:
#            self.make_baseline(linetype=linetype, 
#                               shiftline=shiftline, show_plot=False)
#        else:
#            print("Using previously fit baseline-subtracted absorbance")
#
#        self.abs_nobase_cm = self.subtract_baseline()               
#        peak_idx = scipysignal.find_peaks_cwt(self.abs_nobase_cm, widths)
#        print('change widths=np.arange(1,40) to alter sensitivity')
#        print('peaks found at the following wavenumbers:')
#        print(self.base_wn[peak_idx])
#
#        fig, ax = self.plot_subtractbaseline()
#        for idx in peak_idx:
#            wn = self.base_wn[idx]
#            ax.plot([wn, wn], ax.get_ylim(), '-r', linewidth=1.5)
#        
#
#    def make_average_spectra(self, spectra_list, folder=None):
#        """Takes list of spectra and returns average absorbance (/cm)
#        to the new spectrum (self)"""       
#        list_abs_to_average = []
#        list_wn = []
#        for sp in spectra_list:
#            if sp.abs_full_cm is None:
#                if folder is not None:
#                    sp.filename = ''.join((folder, sp.fname, '.CSV'))
#                else:
#                    sp.filename = ''.join((sp.fname, '.CSV'))
#                if os.path.isfile(sp.filename) is False:
#                    print(sp.filename)
#                    print('File not found')
#                    print()
#                    return False
#                sp.get_data()
#                check = sp.divide_by_thickness()
#                if check is False:
#                    return False
#            abs_to_append = sp.abs_full_cm
#            list_abs_to_average.append(abs_to_append)
#            list_wn.append(sp.wn_full)
#        ave_abs = np.mean(list_abs_to_average, axis=0)
#        waven = np.mean(list_wn, axis=0)
#        self.thickness_microns = 1. # dummy placeholder to show this is not raw data
#        self.wn_full = waven
#        self.abs_full_cm = ave_abs
#        self.start_at_zero()
#    
#    def divide_by_thickness(self):
#        """Divide raw absorbance by thickness"""
#        if self.thickness_microns is None:
#            check = self.set_thick()
#            if check is False:
#                return False
#                
#        if self.fname is None:
#            print('Need .fname to know what to call saved file')
#            return False
#            
#        # Get the data from the file
#        if self.wn_full is None or self.abs_raw is None:
#            self.get_data()
#
#        # Convert from numpy.float64 to regular python float
#        # or else element-wise division doesn't work.
#        if isinstance(self.thickness_microns, float) is True:
#            th = self.thickness_microns
#        elif isinstance(self.thickness_microns, int) is True:
#            th = float(self.thickness_microns)
#        else:
#            th = np.asscalar(self.thickness_microns)
#            
#        self.abs_full_cm = self.abs_raw * 1e4 / th
#        return self.abs_full_cm
#
#    def start_at_zero(self, wn_xlim_left=4000., wn_xlim_right=3000.):
#        """Divide raw absorbance by thickness and
#        shift minimum to 0 within specified wavenumber range specified
#        by wn_xlim_left and _right"""
#        if self.abs_full_cm is None:
#            check = self.divide_by_thickness()
#            if check is False:
#                return False
#
#        index_lo = (np.abs(self.wn_full-wn_xlim_right)).argmin()
#        index_hi = (np.abs(self.wn_full-wn_xlim_left)).argmin()
#
#        indices = list(range(index_lo, index_hi, 1))
#
#        try: 
#            self.abs_full_cm = self.abs_full_cm - min(self.abs_full_cm[indices])
#        except ValueError:
#            print('There was a problem.')
#            print('index_lo at wn', wn_xlim_right, ':', index_lo)
#            print('index_hi at wn', wn_xlim_left, ':', index_hi)
#            return False        
#        return self.abs_full_cm
#
#    def make_baseline(self, linetype='line', shiftline=0.02, 
#                      show_fit_values=False, show_plot=False,
#                      size_inches=(3., 2.5), wn_low=None, wn_high=None,
#                      wn_mid=None, 
#                      splinetype='quadratic', abs_high=None, abs_low=None,
#                      abs_smear_high=0, abs_smear_low=0):
#        """Make baseline that is a spline (linetype='spline'; default),
#        a line (linetype='line'), 
#        or quadratic baseline (linetype='quadratic' with 
#        argmument shiftline determining extent of curvature) 
#        and return baseline absorption curve. Shiftline value determines
#        how much quadratic deviates from linearity"""
#        if self.thickness_microns is not None:
#            if self.abs_full_cm is None:
#                check = self.start_at_zero()
#                if check is False:
#                    return False
#            absorbance = self.abs_full_cm
#        else: 
#            if self.abs_raw is None:
#                self.get_data()
#                check = self.start_at_zero()
#                if check is False:
#                    return False
#            absorbance = self.abs_raw
#
#        if wn_low is not None:
#            self.base_low_wn = wn_low
#        
#        if wn_high is not None:
#            self.base_high_wn = wn_high
#        
#        index_lo = (np.abs(self.wn_full-self.base_low_wn)).argmin()
#        index_hi = (np.abs(self.wn_full-self.base_high_wn)).argmin()        
#        self.base_wn = self.wn_full[index_lo:index_hi]
#
#        # Smearing start and stop over a range of wavenumbers
#        abs_smear_high=int(abs_smear_high)
#        abs_smear_low=int(abs_smear_low)
#        
#        if abs_high is None:
#            if abs_smear_high > 0:
#                gulp_hi_x = list(range(index_hi-abs_smear_high, 
#                                  index_hi+abs_smear_high))
#                gulp_hi_y = absorbance[gulp_hi_x]
#                yhigh = np.mean(gulp_hi_y)
#            else:
#                yhigh = absorbance[index_hi]
#        else:
#            yhigh = abs_high
#            
#        if abs_low is None:
#            if abs_smear_low > 0:
#                gulp_lo_x = list(range(index_lo-abs_smear_low, index_lo+abs_smear_low))
#                gulp_lo_y = absorbance[gulp_lo_x]
#                ylow = np.mean(gulp_lo_y)
#            else:
#                ylow = absorbance[index_lo]
#        else: 
#            ylow = abs_low
#
#        x = np.array([self.wn_full[index_hi], self.wn_full[index_lo]])
#        y = np.array([yhigh, ylow])
#        p = np.polyfit(x, y, 1)
#        
#        if linetype == 'line':
#            base_abs = np.polyval(p, self.base_wn)
#
#        elif linetype == 'quadratic':            
#            # add in a point to fit curve to
#            if wn_mid is not None:      
#                self.base_mid_wn = wn_mid
#                index_mid = (np.abs(self.wn_full-self.base_mid_wn)).argmin()
#                abs_at_wn_mid = absorbance[index_mid]
#                yadd = abs_at_wn_mid
#            elif shiftline is not None:
#                yshift = shiftline
#                yadd = np.polyval(p, self.base_mid_wn) - yshift
#            else:
#                yshift = self.base_mid_yshift
#                yadd = np.polyval(p, self.base_mid_wn) - yshift
#                
#            x = np.insert(x, 1, self.base_mid_wn)
#            y = np.insert(y, 1, yadd)
#            p2 = np.polyfit(x, y, 2)
#            base_abs = np.polyval(p2, self.base_wn)
#
#            if show_fit_values == True:
#                print('fitting x values:', x)
#                print('fitting y values:', y)
#
#        elif linetype == 'spline':
#            idx_max = self.wn_full.argmax()
#            xinterp = np.concatenate((self.wn_full[0:index_lo], 
#                                      self.wn_full[index_hi:idx_max]))                                      
#            yinterp = np.concatenate((absorbance[0:index_lo], 
#                                      absorbance[index_hi:idx_max]))
#            f = interp.interp1d(xinterp, yinterp, kind=splinetype)
#            base_abs = f(self.base_wn)
#        else:
#            print("linetype must be 'line', 'quadratic', or 'spline'")
#            return
#            
#        self.base_abs = base_abs
#        
#        if show_plot is True:
#            fig, ax = self.plot_showbaseline(size_inches=size_inches)
#            if show_fit_values is True:
#                if abs_smear_high > 0:
#                    ax.plot(self.wn_full[gulp_hi_x], 
#                            absorbance[gulp_hi_x], '-y', alpha=0.9)
#                    
#                if abs_smear_low > 0:
#                    ax.plot(self.wn_full[gulp_lo_x], 
#                            absorbance[gulp_lo_x], '-y', alpha=0.9)
#                    
#                if linetype == 'spline':
#                    ax.plot(self.wn_full[0:index_lo], 
#                            absorbance[0:index_lo], '-r')
#                    ax.plot(self.wn_full[index_hi:idx_max], 
#                            absorbance[index_hi:idx_max], '-r')
#                else:
#                    ax.plot(x, y, 'ro', alpha=0.4)
#                    
#            return fig, ax
#        else:
#            return base_abs
#
#    def subtract_baseline(self, linetype='line', bline=None, 
#                          show_plot=False, shiftline=0.02, yhigh=0.1,
#                          wn_low=None, wn_high=None,
#                          baseline_ending='-baseline.CSV'):
#        """Make baseline and return baseline-subtracted absorbance.
#        Preferred baseline used is (1) input bline, (2) self.base_abs, 
#        (3) created based on linetype information"""
#        if self.wn_full is None:
#            self.get_data()
#        
#        if wn_low is not None:
#            self.base_low_wn = wn_low
#        
#        if wn_high is not None:
#            self.base_high_wn = wn_high
#
#        if bline is not None:
#            base_abs = bline
#        elif self.base_abs is not None:
##            print 'Subtracting baseline stored in attribute base_abs'
#            base_abs = self.base_abs
#        else:
#            print('Making', linetype, 'baseline')
#            base_abs = self.make_baseline(linetype=linetype, shiftline=shiftline)
#        
#        index_lo = (np.abs(self.wn_full-self.base_low_wn)).argmin()
#        index_hi = (np.abs(self.wn_full-self.base_high_wn)).argmin()
##            
#        absorbance = self.absorbance_picker()
#        humps = absorbance[index_lo:index_hi]
#
#        # length check
#        if len(humps) > len(base_abs):
#            ndif = len(humps) - len(base_abs)
#            humps = humps[0:-ndif]
#        elif len(humps) < len(base_abs):
#            ndif = len(base_abs) - len(humps)
#            for n in range(ndif):
#                humps = np.append(humps, humps[-1])
#        
#        if len(humps) != len(base_abs):
#            print()
#            print('baseline subtraction failed; returning false')
#            print('len humps', len(humps))
#            print('len base_abs', len(base_abs))
#            return False
#            
#        abs_nobase_cm = humps - base_abs
#
#        if show_plot is True:
#            self.plot_showbaseline()
#
#        return abs_nobase_cm
#
#    def area_under_curve(self, linetype='line', show_plot=True, shiftline=0.02,
#                         printout=True, numformat='{:.1f}', 
#                         require_saved_baseline=False):
#        """Returns area under the curve in cm^2"""
#        if require_saved_baseline is True:
#            self.get_baseline()
#
#        if (self.base_high_wn is None) or (self.base_low_wn is None):
#            print('Need to specify spectrum baseline wavenumber range')
#            return False        
#
#        # Get or make absorbance and wavenumber range for baseline
#        if (self.base_abs is None) or (self.base_wn is None):
#            print('Making baseline...')
#            abs_baseline = self.make_baseline(linetype, shiftline)
#        else:
##            print ('Using previously fit baseline. ' + 
##                   'Use make_baseline to change it.')
#            abs_baseline = self.base_abs
#
#        abs_nobase_cm = self.subtract_baseline(bline=abs_baseline,
#                                               show_plot=show_plot)
#
#        dx = self.base_high_wn - self.base_low_wn
#        dy = np.mean(abs_nobase_cm)
#        area = dx * dy
#        if printout is True:
#            print(self.fname)
#            print('area:', numformat.format(area), '/cm^2')
#            
#        self.area = area
#        return area
#
#    def water(self, phase_name='cpx', calibration='Bell', numformat='{:.1f}',
#              show_plot=True, show_water=True, printout_area=False,
#              shiftline=None, linetype='line'):
#        """Produce water estimate without error for a single FTIR spectrum."""
#        area = self.area_under_curve(linetype, show_plot, shiftline, 
#                                     printout_area)
#        w = area2water(area, phase=phase_name, calibration=calibration)
#                        
#        # output for each spectrum
#        if show_water is True:
#            print(self.fname)
#            print('area:', numformat.format(area), '/cm^2')
#            print(''.join(('water: ', numformat.format(w), ' ppm H2O; *3 = ', 
#                    numformat.format(w*3.), ' ppm H2O')))
#        return area, w*3
#
#    def save_spectrum(self, delim='\t', file_ending='-per-cm.txt', 
#                      folder=None):
#        """Save entire spectrum divided by thickness to file with headings.
#        1st column: wavenumber /cm; 2nd column: absorbance /cm.
#        Formatted for upload to PULI"""
#        if self.fname is None:
#            print('Need .fname to know what to call saved file')
#            return
#        if self.abs_full_cm is None:
#            self.divide_by_thickness(folder=folder)
#
#        if folder is None:
#            folder = self.folder
#            
#        abs_filename = folder + self.fname + file_ending
#        print('Saving:')
#        print(abs_filename)
##        t = ['wavenumber (/cm)', 'absorbance (/cm)']
##        with open(abs_filename, 'w') as abs_file:
##            for item in t:
##                abs_file.write(item+delim)
##            abs_file.write('\n')
#        a = np.transpose(np.vstack((self.wn_full, self.abs_full_cm)))
#        with open(abs_filename, 'w') as abs_file:
#            np.savetxt(abs_file, a, delimiter=delim)
#        
#    def save_baseline(self, folder=None, 
#                 baseline_ending='-baseline.CSV', linetype='line', 
#                 shiftline=0.02, basel=None, delim=',', showplots=False):
#        """Save baseline with baseline-subtracted spectrum 
#        (-baseline) to file. These can be retrieved using get_baseline(). 
#        Use save_spectrum() to save full wavenumber and -per-cm absorbances.
#        """
#        if folder is None:
#            folder = self.folder
#            
#        if self.fname is None:
#            print('Need .fname to know what to call saved file')
#            return
#            
#        if self.base_abs is None:
#            print('making baseline...')
#            basel = self.make_baseline(linetype, shiftline)
#        else:
#            print('using existing baseline')
#            basel = self.base_abs
#
#        if showplots is True:
#            self.plot_showbaseline(baseline=basel)
#        self.abs_nobase_cm = self.subtract_baseline(linetype=linetype, 
#                                                    bline=basel)
#        
#        base_filename = folder + self.fname + baseline_ending
#        print('Saving', self.fname + baseline_ending)
#        
#        t = ['wavenumber (/cm)', 'baseline value (/cm)', 
#             'baseline-subtracted absorbance (/cm)']
#        with open(base_filename, 'w') as base_file:
#            for item in t:
#                base_file.write(item+delim)
#            base_file.write('\n')
#        a = np.transpose(np.vstack((self.base_wn, self.base_abs,
#                                    self.abs_nobase_cm)))
#        with open(base_filename, 'a') as base_file:
#            np.savetxt(base_file, a, delimiter=delim)
#
#    def get_baseline(self, folder=None, delim=',', 
#                     baseline_ending='-baseline.CSV'):
#        """Get baseline and -subtracted spectrum saved using save_baseline().
#        Same data that goes into MATLAB FTIR_peakfit_loop.m
#        Returns baseline absorbances and baseline-subtracted absorbances. 
#        For wavenumber range, it sets spectrum's base_wn attribute"""
#        if folder is None:
#            folder = self.folder
#        filename = ''.join((folder, self.fname, baseline_ending))
#        if os.path.isfile(filename) is False:
##            print ' '
##            print self.fname            
##            print 'use save_baseline() to make -baseline.CSV'
#            return
#        data = np.genfromtxt(filename, delimiter=',', dtype='float', 
#                             skip_header=1)
#
#        print('Baseline information taken from file', baseline_ending)
#        self.base_wn = data[:, 0]
#        self.base_abs = data[:, 1]
#        self.abs_nobase_cm = data[:, 2]
#        self.base_high_wn = max(data[:, 0])
#        self.base_low_wn = min(data[:, 0])
#        return self.base_abs, self.abs_nobase_cm
#        
#    def get_3baselines(self, folder=None, delim=',', 
#                baseline_ending='-3baselines.CSV'):
#        """Returns block of baseline data saved by water_from_spectra()
#        d[:,0] = baseline wavenumbers, d[:,1:4] = baselines, 
#        d[:,4:7] = baseline-subtracted absorbances"""
#        filename = folder + self.fname + baseline_ending
#        if os.path.isfile(filename) is False:
#            print(' ')
#            print(self.fname)            
#            print('Run water_from_spectra() with savebaselines=True')
#            return
#        data = np.genfromtxt(filename, delimiter=',', dtype='float', 
#                             skip_header=1)
#        return data                            
#                
#    def get_peakfit(self, folder=None, delim=',', 
#                    peak_ending='-peakfit.CSV', 
#                    baseline_ending='-baseline.CSV'):
#        """Get individual peaks for spectrum from fname-peakfit.CSV generated
#        in MATLAB using peakfit.m looping script and savefit.m
#        Return curves (assumes Gaussians) and summation"""
#        if folder is None:
#            folder = self.folder
#        filename = folder + self.fname + peak_ending
#        if os.path.isfile(filename) is True:
#            previous_fit = np.loadtxt(filename, delimiter=delim)
#            
#            # sort previous fit so wavenumbers are ascending
#            def getKey(item):
#                return item[0]
#            previous_fit= np.array(sorted(previous_fit, key=getKey))
#
#            peakpos = previous_fit[:, 0]
#            self.numPeaks = len(peakpos)
#            peak_heights = previous_fit[:, 1]
#            peak_widths = previous_fit[:, 2]
#            peak_areas = previous_fit[:, 3]
#            
#            self.peakpos = peakpos
#            self.peak_heights = peak_heights
#            self.peak_widths = peak_widths
#            self.peak_areas = peak_areas
#        else:
#            print(' ')            
#            print(self.fname)
#            print('Remember to savefit.m after fitting in matlab')
#            return False, False
#
#        if self.base_wn is None:
#            self.get_baseline(baseline_ending=baseline_ending)
#        
#        
#    def get_peakcurves(self, folder=None, delim=',', 
#                       peak_ending='-peakfit.CSV', 
#                       baseline_ending='-baseline.CSV'):
#        """Generates Gaussian curves from peakfitting info. 
#        Requires peak_ending and baseline_ending to find the 
#        right files.
#        Returns (1) peakfitcurves, all of the individual gaussians,
#        and (2) summed_spectrum, the sum of all the peakfitcurvesb"""
#        if self.peakpos is None:
#            self.get_peakfit(peak_ending=peak_ending,
#                             baseline_ending=baseline_ending)
#                             
#        peakfitcurves = np.ones([len(self.peakpos), len(self.base_wn)])
#        summed_spectrum = np.zeros_like(self.base_wn)
#        for k in range(len(self.peakpos)):
#            peakfitcurves[k] = make_gaussian(x=self.base_wn, 
#                                            pos=self.peakpos[k],
#                                            h=self.peak_heights[k], 
#                                            w=self.peak_widths[k])
#
#            summed_spectrum += peakfitcurves[k]
#        return peakfitcurves, summed_spectrum
#
#    def get_peakareas(self):
#        """Assign peak area based on Gaussian curves from current peak
#        width and height"""
#        peakfitcurves, summed_spectrum = self.get_peakcurves()
#        for k in range(len(self.peakpos)):
#            dx = self.base_high_wn - self.base_low_wn
#            dy = np.mean(peakfitcurves[k])
#            area = dx * dy
#            self.peak_areas[k] = area
#    
#    def ylim_picker(self, wn_xlim_left=4000, wn_xlim_right=3000, pad_top=0.1, 
#                    pad_bot=0., raw_data=False):
#        """Returns reasonable min and max values for y-axis of
#        plots based on the absorbance values for the specified wavenumber
#        range and padded top and bottom with pad variable"""
#        if self.wn_full is None:
#            self.get_data()
#            
#        if self.thickness_microns is None:
#            absorbance = self.abs_raw           
#        else:
#            self.start_at_zero(wn_xlim_left=wn_xlim_left,
#                               wn_xlim_right=wn_xlim_right)
#            absorbance = self.abs_full_cm
#            
#        idx_lo = (np.abs(self.wn_full-wn_xlim_right)).argmin()
#        idx_hi = (np.abs(self.wn_full-wn_xlim_left)).argmin()
#        
#        y = absorbance[idx_lo:idx_hi]
#
#        bottom = min(y) 
#        top = max(y)
#        ylow = bottom - pad_bot
#        yhigh = top + pad_top
#
#        return ylow, yhigh
#          
#    def plot_spectrum_outline(self, size_inches=(3., 2.5), shrinker=0.15,
#                              figaxis=None, wn_xlim_left=4000., 
#                              wn_xlim_right=3000., pad_top=0.1, 
#                              pad_bot=0., raw_data=False):
#        """Make standard figure outline for plotting FTIR spectra"""
#        if figaxis is None:
#            f, ax = plt.subplots(figsize=size_inches)
#        else:
#            ax = figaxis
#        ax.set_xlabel('Wavenumber (cm$^{-1})$')
#        ax.set_ylabel('Absorbance (cm$^{-1})$')
#        ax.set_title(self.fname)
#        ax.set_xlim(wn_xlim_left, wn_xlim_right)
#        ax.grid()
#        
#        if self.thickness_microns is None:
#            raw_data = True
#
#        ylow, yhigh = self.ylim_picker(wn_xlim_left=wn_xlim_left,
#                                       wn_xlim_right=wn_xlim_right, 
#                                       pad_top=pad_top, pad_bot=pad_bot,
#                                       raw_data=raw_data)
#        ax.set_ylim(ylow, yhigh)
#        
#        box = ax.get_position()
#        ax.set_position([box.x0 + box.width*shrinker, 
#                         box.y0 + box.height*shrinker, 
#                         box.width*(1.0-shrinker), 
#                         box.height*(1.0-shrinker)])
#
#        plt.setp(ax.get_xticklabels(), rotation=45)
#        
#        if figaxis is None:
#            return f, ax
#
#    def plot_spectrum(self, style=styles.style_1, figaxis=None, wn=None,
#                      size_inches=(3., 2.5), label=None, offset=0., 
#                      wn_xlim_left=4000., wn_xlim_right=3000., 
#                      pad_top=0.1, pad_bot=0., plot_raw=False):
#        """Plot the raw spectrum divided by thickness"""
#        if self.wn_full is None:
#            self.get_data()
#
#        if self.thickness_microns is not None:
#            check = self.start_at_zero()
#            if check is False:
#                return
#            absorbance = self.abs_full_cm
#        else:
#            absorbance = self.abs_raw
#            plot_raw = True
#                
#        if figaxis is None:
#            fig, ax = self.plot_spectrum_outline(size_inches=size_inches,
#                                                 wn_xlim_left=wn_xlim_left,
#                                                 wn_xlim_right=wn_xlim_right)
#        else:
#            fig = None
#            ax = figaxis
#
#        if plot_raw is True:
#            ax.set_ylabel('raw absorbance')
#
#        if label is None:
#            label = self.fname
#            
#        style_to_use = style.copy()
#        style_to_use.update({'label' : label})
#        
#        ax.plot(self.wn_full, absorbance + offset, **style_to_use)
#        
#        ax.set_xlim(wn_xlim_left, wn_xlim_right)
#        
#        ylow, yhigh = self.ylim_picker(wn_xlim_left=wn_xlim_left,
#                                       wn_xlim_right=wn_xlim_right, 
#                                       pad_top=pad_top, pad_bot=pad_bot)
#        ax.set_ylim(ylow, yhigh)
#
#        return fig, ax
#
#    def orientation(self, top=None):
#        """guess orientation based on Si-O overtones in FTIR spectrum"""
#        print('\nOrientations for olivine only. See Lemaire et al. 2004 Figure 1.')
#        fig, ax = self.plot_spectrum(pad_top=0.4, wn_xlim_left=2200., 
#                                     wn_xlim_right=1200)
#        fig.set_size_inches(6, 6)
#                                                    
#        if top is not None:
#            ax.set_ylim(0, top)
#        else:
#            ctop = ax.get_ylim()[1]
#            ax.set_ylim(0, ctop + 0.5*ctop)
#
#        ytext = ax.get_ylim()[1] - 0.1*ax.get_ylim()[1]
#        labels = ['E || a', 'E || b', 'E || c']
#        for idx, wn in enumerate([2035, 1670, 1785,]):
#            ax.plot([wn, wn], ax.get_ylim(), '-r')
#            ax.text(wn, ytext, labels[idx], rotation=90, backgroundcolor='w',
#                    va='center', ha='center', fontsize=12)
#        
#            
#        return fig, ax
#
#    def make_composite_peak(self, peak_idx_list):
#        """Make a new 'peak', e.g., for [Ti] in olivine, by summing up 
#        other peaks given by their indexes in peak_idx_list"""
#        peakpos_new = 0.
#        peak_height_new = 0.
#        peak_width_new = 0.
#        peak_area_new = 0.
#        for idx in peak_idx_list:
#            peakpos_new = peakpos_new + self.peakpos[idx]
#            peak_height_new = peak_height_new + self.peak_heights[idx]
#            peak_width_new = peak_width_new + self.peak_widths[idx]
#            peak_area_new = peak_area_new + self.peak_areas[idx]
#    
#
#        self.peakpos = np.append(self.peakpos, peakpos_new)
#        self.peak_heights = np.append(self.peak_heights, peak_height_new)
#        self.peak_widths = np.append(self.peak_widths, peak_width_new)
#        self.peak_areas = np.append(self.peak_areas, peak_area_new)
#        
#    def plot_peakfit(self, peak_ending='-peakfit.CSV', 
#                     baseline_ending='-baseline.CSV',
#                     style=styles.style_spectrum, 
#                     stylesum=styles.style_summed, 
#                     fig_ax = None, legend=True,
#                     stylepeaks=styles.style_fitpeak, top=None, legloc=1):
#        """Single spectrum: Plot peaks fit in MATLAB using peakfit.m
#        REQUIRES the peak_ending and baseline_ending so that it can
#        locate the file spectrum.fname+peak_ending and +baseline_ending.
#        If those aren't present or the file doesn't exist, it will screw 
#        up."""
#        # Take baseline-subtracted spectrum from saved file every time 
#        # to avoid any possible funny business from playing with baselines
#        reload(styles)
#        gaussian, summed_spectrum = self.get_peakcurves(peak_ending=peak_ending, 
#                                                        baseline_ending=baseline_ending)        
#        if gaussian is False:
#            return
#        
#        wn, observed = self.get_baseline(baseline_ending=baseline_ending)
#
#        if fig_ax is None:
#            fig, ax = self.plot_subtractbaseline(style=style, label='') 
#        else: 
#            ax = fig_ax
#            
#        ax.plot(wn, observed, label='observed', **style)
#        ax.plot(wn, gaussian[0], label='fit bands', **stylepeaks)
#        ax.plot(self.base_wn, summed_spectrum, label='peak sum', **stylesum)
#        
#        if fig_ax is None:
#            if legend is True:
#                ax.legend(loc=legloc)
#            ax.set_ylim(0., top)        
#        
#        if top is None:
#            topnat = ax.get_ylim()[1]
#            top = topnat + topnat*0.75
#            
#
#        for k in range(len(self.peakpos)):
#            ax.plot(self.base_wn, gaussian[k], **stylepeaks)
#
#        if fig_ax is None:
#            return fig, ax
#
#    def abs_at_given_wn(self, wn, absorbance='thickness normalized'):
#        """Input wavenumber, output absorbance that is thickness normalized
#        (absorbance='normalized' default), 'raw', or 'baseline-subtracted'"""
#
#        # check you know what the wavenumbers are
#        if self.wn_full is None:
#            check = self.get_data()
#            if check is False:
#                return False
#
#        idx = (np.abs(self.wn_full-wn)).argmin()
#                
#        if absorbance == 'thickness normalized':
#            if self.abs_full_cm is None:
#                self.start_at_zero()
#            abs_at_wn = self.abs_full_cm[idx]
#            return abs_at_wn
#            
#        elif absorbance == 'raw':
#            print('Sorry, raw not programmed in yet')
##            if self.abs_raw is None:
##                self.get_data()
#            return
#        elif absorbance == 'baseline-subtracted':
#            print('Sorry, baseline-subtracted not programmed in yet')
#            return
#        else: 
#            print('absorbance options are:')
#            print('thickness normalized (default)')
#            print('raw')
#            print('baseline-subtracted')
#
#    def absorbance_picker(self):
#        """Is this raw or thickness normalized absorbance you're after?"""
#        if self.thickness_microns is None:
#            if self.abs_raw is None:
#                self.get_data()
#            absorbance = self.abs_raw
#        else:
#            if self.abs_full_cm is None:
#                self.start_at_zero()
#            absorbance = self.abs_full_cm
#        return absorbance            
#        
#    def plot_showbaseline(self, linetype='line', shiftline=None,
#                          wn_baseline=None, abs_baseline=None, 
#                          style=styles.style_spectrum, 
#                          style_base=styles.style_baseline,
#                          figaxis=None, label=None, size_inches=(3., 2.5),
#                          offset=0.0, wn_xlim_left=4000, wn_xlim_right=3000.):
#        """Plot FTIR spectrum and show baseline. 
#        You can pass in your own baseline"""
#        if figaxis is None:
#            fig, ax = self.plot_spectrum_outline(size_inches=size_inches,
#                                                 wn_xlim_left=wn_xlim_left,
#                                                 wn_xlim_right=wn_xlim_right)
#        else:
#            fig = None
#            ax = figaxis
#        
#        # Get or make absorbance for baseline
#        if abs_baseline is None:          
#            if self.base_abs is None:
#                print('Making baseline...')
#                abs_baseline = self.make_baseline(linetype, shiftline)
#            else:
#                abs_baseline = self.base_abs
#
#        # get baseline wavenumber range
#        if wn_baseline is None:
#            if self.base_wn is not None:
#                wn_baseline = self.base_wn
#            else:
#                print(('Need to pass in baseline wavenumber range too ' +
#                       'either here as wn_baseline or as spectrum.base_wn'))
#                return                
#
#        if label is None:
#            label = self.fname
#            
#        style_to_use = style.copy()
#        style_to_use.update({'label' : label})
#        style_base['label'] = 'baseline'
#
#        absorbance = self.absorbance_picker()
#
#        ax.plot(self.wn_full, absorbance + offset, **style_to_use)
#        ax.plot(wn_baseline, abs_baseline + offset, **style_base)
#        
#        if self.thickness_microns is None:
#            ax.set_ylabel('Raw absorbance')
#        return fig, ax
#
#    def plot_subtractbaseline(self, polyorder=1, style=styles.style_spectrum, 
#                              figaxis=None, label=None, offset=0.,
#                              size_inches=(3, 2.5), 
#                              baseline_ending='-baseline.CSV'):
#        """Make and plot baseline-subtracted spectrum"""
#        abs_nobase_cm = self.subtract_baseline(polyorder)
#        
#        if figaxis is None:
#            fig, ax = self.plot_spectrum_outline(size_inches=size_inches)
#        else:
#            fig = None
#            ax = figaxis
#
#        if label is None:
#            label = self.fname
#            
#        style_to_use = style.copy()
#        style_to_use.update({'label' : label})
#            
#        pad = max(abs_nobase_cm)
#        yhigh = pad + 0.1*pad
#        ax.set_ylim(0, yhigh)
#        ax.set_xlim(self.base_high_wn, self.base_low_wn)
#        ax.plot(self.base_wn, abs_nobase_cm+offset, **style_to_use)
#        return fig, ax
#
#    def plot_peakfit_and_baseline(self, style=styles.style_spectrum, 
#                                  stylesum=styles.style_summed, 
#                                  stylepeaks=styles.style_fitpeak, 
#                                  style_base=styles.style_baseline,
#                                  top=None, bottom=0., legloc=1, 
#                                  label_spectrum='observed', 
#                                  peak_ending='-peakfit.CSV',
#                                  baseline_ending='-baseline.CSV'):
#        """Plot spectrum with baseline and peakfit information together"""
#        # Take baseline-subtracted spectrum from saved file every time 
#        # to avoid any possible funny business from playing with baselines
#        if self.peakpos is None:
#            print('need to run get_peakfit first')
#            return
#        if self.base_wn is None:
#            print('need to run get_baseline() first')
#            return
#        
#        gaussian, summed_spectrum = self.get_peakcurves(peak_ending=peak_ending, 
#                                              baseline_ending=baseline_ending)
#            
#        wn = self.base_wn
#
#        stylepeaks['zorder'] = 1 # otherwise peaks show up above baseline
#        
#        fig, ax = self.plot_spectrum(style=style, label=label_spectrum)
#        ax.plot(self.base_wn, self.base_abs, label='baseline', **style_base) 
#        ax.plot(wn, gaussian[0]+self.base_abs, label='fit bands', **stylepeaks)
#        ax.plot(self.base_wn, summed_spectrum+self.base_abs, label='peak sum', 
#                **stylesum)
#                
#        ax.legend(loc=legloc)
#        
#        if top is None:
#            topnat = ax.get_ylim()[1]
#            top = topnat + topnat*0.75
#            
#        ax.set_ylim(bottom, top)
#
#        for k in range(len(self.peakpos)):
#            ax.plot(self.base_wn, gaussian[k]+self.base_abs, **stylepeaks)
#
#        return fig, ax
#
#           
#
#def make_filenames(folder, classname=Spectrum, file_ending='.CSV'):
#    """ Set filename attribute based on folder and fname attribute
#    for all Spectra() with fname but no filename."""
#    for obj in gc.get_objects():
#        if isinstance(obj, classname):
#            try:
#                if obj.fname is not None and obj.filename is None:
#                    obj.filename = ''.join((folder, obj.fname, file_ending))
#            except AttributeError:
#                print('just chill, ok?')
#
#def water_from_spectra(list3, folder, phase='cpx', 
#                       proper3=False, numformat='{:.0f}',
#                       savebaselines=False, show_plots=True, 
#                       baseline_ending='-3baselines.CSV', delim=',', 
#                       calibration='Bell', main_yshift=None,
#                       window_large=None, window_small=None, yhigh=1.):
#    """Produce water estimate from list of FTIR spectra; 
#    Default calibration is Bell et al. 1995, ideally using 
#    3 spectra that are polarized in orthogonal directions (proper3=True)
#    Default (proper3=False) is to estimate area and total water from each spectrum.
#    Then average and std the results."""
#    if calibration == 'Bell':
#        print('Bell calibration')
#        pass
#    elif calibration == 'Paterson':
#        print('Paterson calibration not available quite yet...')
#        return
#    else:
#        print('Only Bell (default) or Paterson calibration so far')
#        return
#        
#    if proper3 == True and len(list3) != 3:
#        print(' ')
#        print('For proper3=True, list should contain only 3 spectra')
#        proper3 = False
#    
#    uarea_list = []
#    uwater_list = []
#    
#    for spec in list3:
#        if show_plots is True:
#            print('Showing spectrum for', spec.fname)
#            fig, ax = spec.plot_spectrum_outline()
#            plt.plot(spec.wn_full, spec.abs_full_cm, **styles.style_spectrum) 
#
#        # how much to shift each quadratic line
#        if main_yshift is None:        
#            main_yshift = spec.base_mid_yshift        
#        if window_large is None:
#            window_large = spec.base_w_large
#        if window_small is None:
#            window_small = spec.base_w_small
#
#        if spec.abs_full_cm is None:
#            spec.start_at_zero()
#        # Generate list of 3 possible areas under the curve
#        area_list = np.array([])
#        spec.make_baseline()
#        baseline_list = np.ones([3, len(spec.base_wn)])
#        pek_list = np.ones([3, len(spec.base_wn)])
##        
#        k = 0
#        for bam in [window_small, -window_small, -window_large]:
#            # Setting up 3 different baseline shapes
#            main_yshift = main_yshift - bam
#            base_abs = spec.make_baseline('quadratic', shiftline=main_yshift)
#            abs_nobase_cm = spec.subtract_baseline(bline=base_abs)
#            if show_plots is True:
#                ax.plot(spec.base_wn, base_abs, **styles.style_baseline)
#
#            area = spec.area_under_curve(area_plot=False, printout=False)
#            area_list = np.append(area_list, area)
#            baseline_list[k,:] = base_abs
#            pek_list[k,:] = abs_nobase_cm
#            k+=1
#
#        # list areas and water with uncertainties
#        uarea = ufloat(np.mean(area_list), np.std(area_list))
#        uarea_list.append(uarea)
#        uwater = area2water(uarea, phase=phase)
#
#        # Water estimate depends on if you are adding 3 areas 
#        # or just multiplying by three and later averaging
#        if proper3 is False:
#            uwater = 3*uwater
#        uwater_list.append(uwater)
#
#        # I show() the figures explicitly in this first loop because
#        # otherwise they show up inline *after* all the printed output
#        if show_plots is True:
#            if min(base_abs) > 0:
#                ylow = 0
#            else:
#                ylow = min(base_abs) + 0.1*min(base_abs)
#            plt.ylim(ylow, yhigh)
#            plt.show(fig)
#
#        if savebaselines is True:
#            base_filename = folder + spec.fname + baseline_ending
#            t = ['wavenumber (/cm)', 'upper baseline (/cm) (smaller area)', 
#                 'main baseline (/cm)', 'lower baseline (/cm) (larger area)',
#                 'absorbance - upper baseline (/cm)', 
#                 'absorbance - main baseline (/cm)',
#                 'absorbance - lower baseline (/cm)']
#            with open(base_filename, 'w') as base_file:
#                for item in t:
#                    base_file.write(item+delim)
#                base_file.write('\n')
#            a = np.transpose(np.vstack((spec.base_wn, baseline_list, pek_list)))
#            with open(base_filename, 'a') as base_file:
#                np.savetxt(base_file, a, delimiter=delim)
#            
#    # output for each spectrum
#    for x in range(len(list3)):
#        print(list3[x].fname)
#        print("area:", numformat.format(uarea_list[x]), "/cm^2")
#        print("water:", numformat.format(uwater_list[x]), "ppm H2O")
#        if savebaselines is True:
#            print('Saved baselines to', list3[x].fname+baseline_ending)
#        print(' ')
#    # Final output
#    if proper3 == True:
#        a = np.sum(uarea_list)
#        w = np.sum(uwater_list)
#        print('Sum of individual contributions')
#    else:
#        a = np.mean(uarea_list)
#        w = np.mean(uwater_list)
#        print('Averaged individual (*3) estimates')
#    print('area:', numformat.format(a), '/cm^2')
#    print('water:', numformat.format(w), 'ppm H2O')
#    print(' ')
#    return a, w  
#
#def list_with_attribute(classname, attributename, attributevalue):
#    """Gather all instances in specified class with a particular attribute"""
#    my_list = []
#    for obj in gc.get_objects():
#        if isinstance(obj, classname):
#            if getattr(obj, attributename) == attributevalue:
#                my_list.append(obj.fname)
#    return my_list
#
