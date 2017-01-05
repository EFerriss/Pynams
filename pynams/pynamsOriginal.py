"""
FTIR class definitions and functions for processing and plotting FTIR spectra
that are saved in .CSV format either on their own, as profiles, and/or
 in groups of three othogonal profiles ('whole-block' data).
Also includes some diffusion modeling (pynams.diffusion) and 
using my style preferences (pynams.styles).

The main unit is a class called Spectrum, which includes attributes 
such as sample thickness and lists of wavenumbers and absorbances. 
A few defaults are set up that are geared toward H in nominally anhydrous
minerals (NAMs) such a plotting wavenumber range from 3000 to 4000 /cm and
creating baselines between 3200 and 3700 /cm. 

Copyright (c) 2015 Elizabeth Ferriss

This software was developed using Python 2.7.

"""
#import styles
#import diffusion
#import gc
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.lines as mlines

from .uncertainties import ufloat

print '\nIf you see this message and a +/- number, my test worked!'
y = ufloat(3, 0.5)
print y

#%%
#import os.path
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


#def make_gaussian(pos, h, w, x=np.linspace(3000, 4000, 150)):
#    """Make and return Gaussian curve over specified range"""
#    y = h * np.e**(-((x-pos) / (0.6005615*w))**2)
#    return y
#
#def absorption_coefficients(phase, calibration, peak_idx=None):
#    """Input phase, calibration, and peak 
#    Returns absorption coefficient"""
#    x = ufloat(5, 0.2)
#    print x
#    if (calibration == 'Bell') and (phase=='cpx'):
#        absorption_coeff = 1.0 / ufloat(7.09, 0.32)
#
#    # Bell et al. 2003        
#    elif (calibration == 'Bell') and (phase=='olivine'):
#        absorption_coeff = ufloat(0.188, 0.012)
#        
#    # Withers et al. 2012
#    elif (calibration == 'Withers') and (phase=='olivine'):
#        absorption_coeff = ufloat(0.119, 0.006)
#        
#    else:
#        print 'Calibrations supported so far:',
#        print '   Bell et al. 1995 for cpx'
#        print '   Bell et al. 2003 for olivine'
#        print '   Withers et al. 2012 for olivine'
#        return
#
#    return absorption_coeff        
#    
#def area2water(area_cm2, phase='cpx', calibration='Bell', peak_idx=None):
#    """Takes area in cm-2, multiplies by absorption coefficient, 
#    return water concentration in ppm H2O"""
#    absorption_coeff = absorption_coefficients(phase=phase, 
#                                               calibration=calibration,
#                                               peak_idx=peak_idx)
#    w = absorption_coeff * area_cm2
#    return w
#        
##%% Define classes and functions related to FTIR spectra        
#class Spectrum():
#    def __init__(self, fname, folder='', sample=None, filename=None,
#                 thick_microns=None, raypath=None,
#                 base_low_wn=3200, base_high_wn=3700,
#                 base_mid_wn=3550, polar=None, other_name=None,
#                 filetype='.CSV'):
#        self.fname = fname
#        self.folder = folder
#        self.sample = sample
#        self.filename = filename
#        self.thick_microns = thick_microns
#        self.base_high_wn = base_high_wn
#        self.base_low_wn = base_low_wn
#        self.base_mid_wn = base_mid_wn
#        self.polar = polar 
#        self.other_name = other_name
#        self.folder = folder
#        self.filetype = filetype
#        self.raypath = raypath
#        
#        if self.filename is None and self.fname is not None:
#            self.filename = self.folder + self.fname + self.filetype
#
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
#        """Set spectrum thick_microns from raypath and sample.thick_microns.
#        Similar utility exists for Profiles"""
#        if self.sample is None and self.thick_microns is None:
#            print 'Estimating sample thickness from SiO overtones'
#            self.thick_microns = self.thickness_from_SiO()
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
#                    print 'Need spectrum raypath a, b, or c'
#                    return False
#            else:
#                print 'sample.thickness_microns must have 1 or 3 values'
#                return False
#        self.thick_microns = th
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
#        self.thick_microns = thickness_microns
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
#                print 'basenumber range not established. Check wavenumber range'
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
#            print "Using previously fit baseline-subtracted absorbance"
#
#        self.abs_nobase_cm = self.subtract_baseline()               
#        peak_idx = scipysignal.find_peaks_cwt(self.abs_nobase_cm, widths)
#        print 'change widths=np.arange(1,40) to alter sensitivity'
#        print 'peaks found at the following wavenumbers:'
#        print self.base_wn[peak_idx]
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
#                    print sp.filename
#                    print 'File not found'
#                    print
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
#        self.thick_microns = 1. # dummy placeholder to show this is not raw data
#        self.wn_full = waven
#        self.abs_full_cm = ave_abs
#        self.start_at_zero()
#    
#    def get_data(self):
#        """Get the data from file"""
#        if self.filename is None:
#            self.filename = self.folder + self.fname + self.filetype
#
#        if os.path.isfile(self.filename):
#            if self.filetype == '.CSV':
#                signal = np.loadtxt(self.filename, delimiter=',')
#            elif self.filetype == '.txt':
#                try:
#                    signal = np.loadtxt(self.filename, delimiter='\t', 
#                                        dtype=None) 
#                except ValueError:
#                    print '\nProblem reading this file format'
#                    return False
#            else:
#                print 'check filetype'
#        else:
#            print 'There is a problem finding the file.'
##            print 'You may need to run pynams.make_filenames(folder=)'
#            print 'filename =', self.filename
#            return False
#        
#        # sort signal by wavenumber
#        signal = signal[signal[:,0].argsort()]
#        self.wn_full = signal[:, 0]
#        self.abs_raw = signal[:, 1]
#        return True
#    
#    def divide_by_thickness(self):
#        """Divide raw absorbance by thickness"""
#        if self.thick_microns is None:
#            check = self.set_thick()
#            if check is False:
#                return False
#                
#        if self.fname is None:
#            print 'Need .fname to know what to call saved file'
#            return False
#            
#        # Get the data from the file
#        if self.wn_full is None or self.abs_raw is None:
#            self.get_data()
#
#        # Convert from numpy.float64 to regular python float
#        # or else element-wise division doesn't work.
#        if isinstance(self.thick_microns, float) is True:
#            th = self.thick_microns
#        elif isinstance(self.thick_microns, int) is True:
#            th = float(self.thick_microns)
#        else:
#            th = np.asscalar(self.thick_microns)
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
#        indices = range(index_lo, index_hi, 1)
#
#        try: 
#            self.abs_full_cm = self.abs_full_cm - min(self.abs_full_cm[indices])
#        except ValueError:
#            print 'There was a problem.'
#            print 'index_lo at wn', wn_xlim_right, ':', index_lo
#            print 'index_hi at wn', wn_xlim_left, ':', index_hi
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
#        if self.thick_microns is not None:
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
#                gulp_hi_x = range(index_hi-abs_smear_high, 
#                                  index_hi+abs_smear_high)
#                gulp_hi_y = absorbance[gulp_hi_x]
#                yhigh = np.mean(gulp_hi_y)
#            else:
#                yhigh = absorbance[index_hi]
#        else:
#            yhigh = abs_high
#            
#        if abs_low is None:
#            if abs_smear_low > 0:
#                gulp_lo_x = range(index_lo-abs_smear_low, index_lo+abs_smear_low)
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
#                print 'fitting x values:', x
#                print 'fitting y values:', y
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
#            print "linetype must be 'line', 'quadratic', or 'spline'"
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
#            print 'Making', linetype, 'baseline'
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
#            print
#            print 'baseline subtraction failed; returning false'
#            print 'len humps', len(humps)
#            print 'len base_abs', len(base_abs)
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
#            print 'Need to specify spectrum baseline wavenumber range'
#            return False        
#
#        # Get or make absorbance and wavenumber range for baseline
#        if (self.base_abs is None) or (self.base_wn is None):
#            print 'Making baseline...'
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
#            print self.fname
#            print 'area:', numformat.format(area), '/cm^2'
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
#            print self.fname
#            print 'area:', numformat.format(area), '/cm^2'
#            print ''.join(('water: ', numformat.format(w), ' ppm H2O; *3 = ', 
#                    numformat.format(w*3.), ' ppm H2O'))
#        return area, w*3
#
#    def save_spectrum(self, delim='\t', file_ending='-per-cm.txt', 
#                      folder=None):
#        """Save entire spectrum divided by thickness to file with headings.
#        1st column: wavenumber /cm; 2nd column: absorbance /cm.
#        Formatted for upload to PULI"""
#        if self.fname is None:
#            print 'Need .fname to know what to call saved file'
#            return
#        if self.abs_full_cm is None:
#            self.divide_by_thickness(folder=folder)
#
#        if folder is None:
#            folder = self.folder
#            
#        abs_filename = folder + self.fname + file_ending
#        print 'Saving:'
#        print abs_filename
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
#            print 'Need .fname to know what to call saved file'
#            return
#            
#        if self.base_abs is None:
#            print 'making baseline...'
#            basel = self.make_baseline(linetype, shiftline)
#        else:
#            print 'using existing baseline'
#            basel = self.base_abs
#
#        if showplots is True:
#            self.plot_showbaseline(baseline=basel)
#        self.abs_nobase_cm = self.subtract_baseline(linetype=linetype, 
#                                                    bline=basel)
#        
#        base_filename = folder + self.fname + baseline_ending
#        print 'Saving', self.fname + baseline_ending
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
#        print 'Baseline information taken from file', baseline_ending
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
#            print ' '
#            print self.fname            
#            print 'Run water_from_spectra() with savebaselines=True'
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
#            print ' '            
#            print self.fname
#            print 'Remember to savefit.m after fitting in matlab'
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
#        if self.thick_microns is None:
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
#        if self.thick_microns is None:
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
#        if self.thick_microns is not None:
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
#        print '\nOrientations for olivine only. See Lemaire et al. 2004 Figure 1.'
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
#            print 'Sorry, raw not programmed in yet'
##            if self.abs_raw is None:
##                self.get_data()
#            return
#        elif absorbance == 'baseline-subtracted':
#            print 'Sorry, baseline-subtracted not programmed in yet'
#            return
#        else: 
#            print 'absorbance options are:'
#            print 'thickness normalized (default)'
#            print 'raw'
#            print 'baseline-subtracted'
#
#    def absorbance_picker(self):
#        """Is this raw or thickness normalized absorbance you're after?"""
#        if self.thick_microns is None:
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
#                print 'Making baseline...'
#                abs_baseline = self.make_baseline(linetype, shiftline)
#            else:
#                abs_baseline = self.base_abs
#
#        # get baseline wavenumber range
#        if wn_baseline is None:
#            if self.base_wn is not None:
#                wn_baseline = self.base_wn
#            else:
#                print ('Need to pass in baseline wavenumber range too ' +
#                       'either here as wn_baseline or as spectrum.base_wn')
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
#        if self.thick_microns is None:
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
#            print 'need to run get_peakfit first'
#            return
#        if self.base_wn is None:
#            print 'need to run get_baseline() first'
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
#                print 'just chill, ok?'
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
#        print 'Bell calibration'
#        pass
#    elif calibration == 'Paterson':
#        print 'Paterson calibration not available quite yet...'
#        return
#    else:
#        print 'Only Bell (default) or Paterson calibration so far'
#        return
#        
#    if proper3 == True and len(list3) != 3:
#        print ' '
#        print 'For proper3=True, list should contain only 3 spectra'
#        proper3 = False
#    
#    uarea_list = []
#    uwater_list = []
#    
#    for spec in list3:
#        if show_plots is True:
#            print 'Showing spectrum for', spec.fname
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
#        print list3[x].fname
#        print "area:", numformat.format(uarea_list[x]), "/cm^2"
#        print "water:", numformat.format(uwater_list[x]), "ppm H2O"
#        if savebaselines is True:
#            print 'Saved baselines to', list3[x].fname+baseline_ending
#        print ' '
#    # Final output
#    if proper3 == True:
#        a = np.sum(uarea_list)
#        w = np.sum(uwater_list)
#        print 'Sum of individual contributions'
#    else:
#        a = np.mean(uarea_list)
#        w = np.mean(uwater_list)
#        print 'Averaged individual (*3) estimates'
#    print 'area:', numformat.format(a), '/cm^2'
#    print 'water:', numformat.format(w), 'ppm H2O'
#    print ' '
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
##
## Define classes and functions for working with profiles of spectra
##
#
#class Profile():
#    def __init__(self, profile_name=None, time_seconds=None, folder='',
#                 fname_list=[], positions_microns = np.array([]),
#                 sample=None, direction=None, raypath=None, short_name=None,
#                 spectra_list=[], set_thickness=False,
#                 initial_profile=None, base_low_wn=None, base_high_wn=None,
#                 diffusivity_log10m2s=None, diff_error=None, length_microns=None,
#                 peak_diffusivities=[], peak_diff_error=[], thick_microns=None):
#        """fname_list = list of spectra filenames without the .CSV extension.
#        Raypath and direction expressed as 'a', 'b', 'c' with thickness/length
#        info contained in sample's length_a_microns, length_b_microns, and length_c_microns.
#        base_low_wn and base_high_wn can be used to set the wavenumber
#        range of the baseline for the spectra.
#        
#        """
#        self.profile_name = profile_name
#        self.spectra_list = spectra_list
#        self.folder = folder
#        self.fname_list = fname_list
#        self.positions_microns = positions_microns
#        self.sample = sample
#        self.length_microns = length_microns
#        self.direction = direction
#        self.raypath = raypath
#        self.initial_profile = initial_profile
#        self.short_name = short_name
#        self.time_seconds = time_seconds
#        self.diffusivity_log10m2s = diffusivity_log10m2s
#        self.diff_error = diff_error
#        self.peak_diffusivities = peak_diffusivities
#        self.peak_diff_error = peak_diff_error
#        self.thick_microns = thick_microns
#        
##        if (self.fname_list is not None) and (self.sample is not None):
#        if base_low_wn is not None:
#            for spectrum in self.spectra_list:
#                spectrum.base_low_wn = base_low_wn
#
#        if base_high_wn is not None:
#            for spectrum in self.spectra_list:
#                spectrum.base_low_wn = base_high_wn
#
#        self.make_spectra_list(set_thickness=set_thickness)
#
#    short_name = None # short string for saving diffusivities, etc.
#    thick_microns_list = None
#
#    
#    # for constructing whole-block profiles
#
#    # The actual list of spectra is made automatically by make_spectra_list.
#    # Set spectrum_class_name to use non-default spectra classes
#    # e.g., to set different baseline limits consistently
#    spectra_list = []
#    spectrum_class_name = None 
#    avespec = None # averaged spectra made by self.average_spectra()
#    iavespec = None # initial averaged spectra
#    length_microns = None # length, but I usually use set_len() directly each time
#    waters_list = []
#    waters_errors = []
#    areas_list = None
#    bestfitline_areas = None
#    # Bulk whole-block 3D-WB information and diffusivities if applicable
#    wb_areas = None 
#    wb_water = None
#    # for plotting 
#    style_base = styles.style_profile
#    style_x_marker = None
#    style_y_marker = None
#    style_z_marker = None
#    style_x_line = None
#    style_y_line = None
#    style_z_line = None
#    # peak fitting. i = initial
#    peakpos = None
#    peak_heights = None
#    peak_widths = None
#    peak_areas = None
#    peak_iheights = None
#    peak_iwidths = None
#    peak_iareas = None
#    peak_wb_areas = None
#    peak_wb_heights = None
#    peak_wb_widths = None
#    peak_initial_wb_ratio = 1.0    
#    # Maximum values to scale up to with diffusion curves
#    maximum_area = None
#    maximum_wb_area = 1.0
#    peak_maximum_areas = None
#    peak_maximum_heights = None
#    peak_maximum_areas_wb = None
#    peak_maximum_heights_wb = None
#    
#    # Diffusivities (D) and associated errors all in log10 m2/s
#    # D's can be derived from any or all of the following
#    # - absolute area bulk H: D_area
#    # - absolute area peak-specific: D_peakarea
#    # - absolute height peak-specific: D_height
#    # - whole-block area bulk: D_area_wb 
#    # - whole-block area peak-specific: D_peakarea_wb
#    # - whole-block height peak-specific: D_height_wb
#
#    # absolute value-derived diffusivities
#    D_area = 0.
#    D_peakarea = None
#    D_height = None
#
#    D_area_error = 0.
#    D_peakarea_error = None
#    D_height_error = None
#
#    # whole-block-derived diffusivities
#    D_area_wb = 0.
#    D_peakarea_wb = None
#    D_height_wb = None
#
#    D_area_wb_error = 0.
#    D_peakarea_wb_error = None
#    D_height_wb_error = None
#    
#    def set_all_thicknesses_from_SiO(self):
#        """Individually set thicknesses for all spectra based on the area
#        under their Si-O overtone peaks"""
#        self.make_spectra_list()
#
#    def set_len(self):
#        """Set profile.length_microns from profile.direction and 
#        profile.sample.thick_microns""" 
#
#        if self.sample is None:
#            print '\n', self.profile_name
#            print 'Need to specify profile sample\n'
#            return False
#        else:
#            s = self.sample
#
#        if s.thickness_microns is None:
#            s.thickness_microns = get_3thick(s)                
#        
#        if self.direction == 'a':
#           self.length_microns = s.thickness_microns[0]
#        elif self.direction == 'b':
#            self.length_microns = s.thickness_microns[1]
#        elif self.direction == 'c':
#            self.length_microns = s.thickness_microns[2]
#        else:
#            print 'Set direction of profile to a, b, or c'
#
#        return self.length_microns
#
#    def set_thick(self):
#        """Set profile.thick_microns from profile.raypath and
#        profile.sample.thick_microns"""
#        if self.sample is None:
#            print 'Need to specify profile sample or thickness'
#            return False
#        else:
#            s = self.sample
#        
#        if s.thickness_microns is None:
#            s.thickness_microns = get_3thick(s)
#
#        if self.raypath == 'a':
#           self.thick_microns = s.thickness_microns[0]
#        elif self.raypath == 'b':
#            self.thick_microns = s.thickness_microns[1]
#        elif self.raypath == 'c':
#            self.thick_microns = s.thickness_microns[2]
#        else:
#            print 'Need raypath'
#            return False
#            
#        return self.thick_microns
#
#    def plot_thicknesses(self, figaxis=None):
#        """Plot thickness across profile"""
#        if self.length_microns is None:
#            self.length_microns = max(self.positions_microns)+1.
#        if figaxis is not None:
#            ax = figaxis
#        else:
#            fig, ax, ax_right = self.plot_area_profile_outline(centered=False)
#        ax.plot(self.positions_microns, self.thick_microns_list, 'o')
#        ax.set_ylabel('thickness ($\mu$m)')
#        ax.set_title(self.profile_name)
#        ax.set_ylim(min(self.thick_microns_list)-0.05*min(self.thick_microns_list), 
#                    max(self.thick_microns_list)+0.05*max(self.thick_microns_list))
#        return fig, ax            
#
#    def make_spectra_list(self, set_thickness=True):
#        """Set profile length and generate spectra_list 
#        with key attributes"""
#        try:
#            if (self.raypath is not None) and (self.raypath == self.direction):
#                print "raypath cannot be the same as profile direction"
#                return False
#        except AttributeError:
#            self.raypath = None
#
#        # construct each spectrum from fnames
#        if len(self.spectra_list) == 0:
#            if len(self.fname_list) == 0:
#                print 'Need fnames'
#                return False                
#            fspectra_list = []
#            for x in self.fname_list:
#                newspec = Spectrum(fname=x, folder=self.folder)
#                newspec.fname = x
#                newspec.thick_microns = self.thick_microns
#                fspectra_list.append(newspec)
#            self.spectra_list = fspectra_list
#
#        # set sample, raypath for all
#        for spec in self.spectra_list:
#            spec.sample = self.sample
#            spec.raypath = self.raypath
#
#        if set_thickness is True:
#            self.set_thicknesses()
#        return
#
#    def set_thicknesses(self):
#        """Sets thickness for each spectrum and makes list of thickness for profile"""
#        if self.thick_microns_list is None:
#            self.thick_microns_list = []
#        for spec in self.spectra_list:
#           if len(self.thick_microns_list) < len(self.spectra_list):
#                if self.thick_microns is not None:      
#                    spec.thick_microns = self.thick_microns
#                elif spec.thick_microns is None:
#                    spec.thickness_from_SiO()
#                self.thick_microns_list.append(spec.thick_microns)   
#
#    def set_length(self):
#        if self.length_microns is None:
#            self.length_microns = max(self.thick_microns_list) + 50.
#
#    def average_spectra(self):
#        """Creates and returns averaged spectrum and stores it in
#        attribute avespec"""
#        spec_list = self.spectra_list
#        avespec = Spectrum(folder=None, fname=None)
#        avespec.make_average_spectra(spec_list, folder=self.folder)
#        
#        if self.profile_name is not None:
#            avespec.fname = (self.profile_name + '\naverage profile')
#        else:
#            avespec.fname = 'average profile'
#            
#        avespec.base_high_wn = spec_list[0].base_high_wn
#        avespec.base_low_wn = spec_list[0].base_low_wn
#        avespec.start_at_zero()
#        self.avespec = avespec
#        return avespec
#
#    def plot_spectra(self, show_baseline=True, show_initial_ave=True, 
#                     show_final_ave=True, plot_all=False, 
#                     initial_and_final_together=False, style=styles.style_spectrum, 
#                     stylei=styles.style_initial, wn=None,
#                     figsize=(3.2, 3.2)):
#        """Plot averaged spectrum across profile. Returns figure and axis."""
#        for spec in self.spectra_list:
#            spec.plot_spectrum()
##        if self.spectra_list is None or len(self.spectra_list) < 1:
##            self.make_spectra_list()
##
##        if show_initial_ave is True or initial_and_final_together is True:
##            if self.initial_profile is None:
##                print 'No initial_profile attribute specified'
##                show_initial_ave = False
##                initial_and_final_together = False
##        
##        f = None
##        ax = None
##        
##        if plot_all is True:
##            for spec in self.spectra_list:
##                if show_baseline is True:
##                    spec.plot_showbaseline(style=style)
##                else:
##                    spec.plot_spectrum(style=style, wn=wn)
##            if show_initial_ave is True:
##                for spec in self.initial_profile.spectra_list:
##                    if show_baseline is True:
##                        spec.plot_showbaseline(style=stylei)
##                    else:
##                        spec.plot_spectrum(style=stylei, wn=wn)
##                        
##        if show_final_ave is True or initial_and_final_together is True:
##            avespec = self.average_spectra()
##        
##        if show_initial_ave is True or initial_and_final_together is True:
##            initspec = self.initial_profile.average_spectra()
##            
##        # Plot initial average
##        if show_initial_ave is True:
##            if show_baseline is True:
##                initspec.plot_showbaseline(style=stylei)
##            else:
##                initspec.plot_spectrum(style=stylei, wn=wn)
##
##        # Plot final average
##        if show_final_ave is True:
##            if show_baseline is True:
##                f, ax = avespec.plot_showbaseline(style=style)
##            else:            
##                f, ax = avespec.plot_spectrum(style=style, wn=wn)
##        
##        # Plot average spectra together
##        if initial_and_final_together is True:
##            f, ax = avespec.plot_spectrum_outline()
##            ax.plot(avespec.wn_full, avespec.abs_full_cm, label='Final', **style)
##            ax.plot(initspec.wn_full, initspec.abs_full_cm, label='Initial', **stylei)            
##            ax.legend()
##            tit = self.profile_name + '\nAverage profiles'
##            ax.set_title(tit)
##            if wn is not None:
##                ax.plot([wn, wn], [ax.get_ylim()[0], ax.get_ylim()[1]], 
##                        color='r')
##        return f, ax
#
#    def change_baseline(self, highwn=3800, lowwn=3000, shift=None):
#        """Change baseline parameters for all spectra, final and initial"""
#        for spectrum in self.spectra_list + self.initial_profile.spectra_list:
#            spectrum.base_high_wn = highwn
#            spectrum.base_low_wn = lowwn
#            if shift is not None:
#                spectrum.base_mid_yshift = shift
#
#    def make_composite_peak(self,peak_idx_list):
#        """Make composite peaks for all spectra in profile"""
#        for spec in self.spectra_list:
#            spec.make_composite_peak(peak_idx_list)
#        self.get_peak_info()
#
#    def make_baselines(self, linetype='line', shiftline=None, 
#                       wn_high=3700., wn_low=3200., wn_mid=None,
#                       show_fit_values=False, show_plot=False,
#                       size_inches=(3., 2.5), abs_high=None):
#        """Make baselines for all final and initial spectra in profile"""
#        if len(self.spectra_list) < 1:
#            self.make_spectra_list()
#        for spectrum in self.spectra_list:
#            spectrum.base_high_wn = wn_high
#            spectrum.base_low_wn = wn_low
#            spectrum.make_baseline(linetype=linetype, shiftline=shiftline,
#                                    show_fit_values=show_fit_values,
#                                    show_plot=show_plot, wn_mid=wn_mid,
#                                    size_inches=size_inches)
#
#    def get_baselines(self, initial_too=False, folder=None, delim=',', 
#                      baseline_ending='-baseline.CSV'):
#        """Get previously saved baselines for all spectra in profile"""
#        for spectrum in self.spectra_list:
#            spectrum.get_baseline(baseline_ending=baseline_ending,
#                                  folder=folder, delim=delim)
#            
#        if initial_too is True:
#            for spectrum in self.initial_profile.spectra_list:
#                spectrum.get_baseline()
#
#    def matlab(self):
#        """Print out spectra fnames for FTIR_peakfit_loop.m"""
#        string = "{"
#        for spec in self.spectra_list:
#            stringname = spec.fname
#            string = string + "'" + stringname + "' "
#        string = string + "};"
#        print '\nfilenames ready for FTIR_peakfit_loop.m:'
#        print string
#                        
#    def save_baselines(self, printnames=True):
#        """Save all baselines in profile"""
#        for spectrum in self.spectra_list:
#            spectrum.save_baseline()
#        
#        if printnames is True:
#            self.print_names4matlab()
#            
#    def make_area_list(self, polyorder=1, show_plot=False, set_class=None,
#                       shiftline=None, printout_area=False, peak=None):
#        """Make list of areas (no error) under the curve for an FTIR profile.
#        Default is bulk area. Set peak=wavenumber for peak-specific profile"""
#        # You need the list of spectra for the profile
#        if len(self.spectra_list) < 1:
#            check = self.make_spectra_list(class_from_module=set_class)
#            if check is False:
#                return False
#              
#        areas = []
#        if peak is None:
#            if self.areas_list is not None:
#                print 'Overwriting previous areas list'
#                            
#            if self.spectra_list[0].area is None:                
#                print 'generating bulk areas under curves...'
#                for spec in self.spectra_list:
#                    a = spec.area_under_curve(polyorder, show_plot, shiftline, 
#                                           printout_area, 
#                                           require_saved_baseline=False)
#                    areas.append(a)            
#                self.areas_list = np.array(areas)
#            else:
#                for spec in self.spectra_list:
#                    areas.append(spec.area)
#                
#            self.areas_list = areas
#            
#        else:
#            peaklist = list(self.spectra_list[0].peakpos)
#            print 'peak at', peak
#
#            if peak in peaklist:
#                idx = peaklist.index(peak)
#                for x in self.spectra_list:
#                    a = x.peak_areas[idx]
#                    areas.append(a)
#            else:
#                print 'No peak at wavenumber', peak
#                print 'peak positions:', peaklist
#                return False
#        return areas
#        
#    def get_peakfit(self, peak_ending='-peakfit.CSV',
#                    baseline_ending='-baseline.CSV'):
#        """Get fit peaks from MATLAB for all spectra in profile, including 
#        in the initial profile. The resulting numpy arrays are in dimensions
#        (number of peaks, number of spectra in profile) and stored in
#        the profiles's attributes peak_heights, peak_widths, and peak_areas"""
#        for spectrum in self.spectra_list:
#            spectrum.get_peakfit(peak_ending=peak_ending,
#                                 baseline_ending=baseline_ending)
#        self.get_peak_info()
#            
#    def get_peak_info(self):
#        """Pull peak info from individual spectra into a single profile
#        attribute"""
#        if len(self.spectra_list) < 1:
#            self.make_spectra_list()
#        peakpos = self.spectra_list[0].peakpos
#
#        if peakpos is None:
#            return False
#        hbig = []
#        wbig = []
#        abig = []
#        
#        for p in range(len(peakpos)):
#            h = []
#            w = []
#            a = []        
#            for spectrum in self.spectra_list:
#                h.append(spectrum.peak_heights[p])
#                w.append(spectrum.peak_widths[p])
#                a.append(spectrum.peak_areas[p])
#            hbig.append(h)
#            wbig.append(w)
#            abig.append(a)
#                
#        self.peak_heights = np.array(hbig)
#        self.peak_widths = np.array(wbig)
#        self.peak_areas = np.array(abig)
#        self.peakpos = peakpos
#        
#        self.peak_wb_areas = np.zeros_like(peakpos)
#        self.D_peakarea = np.zeros_like(peakpos)
#        self.D_peakarea_error = np.zeros_like(peakpos)
#        self.D_peakarea_wb = np.zeros_like(peakpos)
#        self.D_peakarea_wb_error = np.zeros_like(peakpos)
#        self.peak_maximum_areas = np.zeros_like(peakpos)
#        self.peak_maximum_areas_wb = np.ones_like(peakpos)        
#        
#        self.D_height = np.zeros_like(peakpos)
#        self.D_height_error = np.zeros_like(peakpos)
#        self.D_height_wb = np.zeros_like(peakpos)
#        self.D_height_wb_error = np.zeros_like(peakpos)
#        self.peak_maximum_heights = np.zeros_like(peakpos)
#        self.peak_maximum_heights_wb = np.ones_like(peakpos)
#
#    def get_area_total(self):
#        """Sum up total area of all peaks"""
#        total_peak_area = []
#        total_height = []
#        
#        for spec in self.spectra_list:
#            spec.get_peakareas()
#            total_peak_area.append(sum(spec.peak_areas))
#            total_height.append(sum(spec.peak_heights))
#        total_area = np.mean(total_peak_area)
#        total_h = np.mean(total_height)
#        self.total_peak_area = total_area
#        self.total_peak_height = total_h
#        return total_area, total_h
#            
#    def print_peakfits(self):
#        """Print out peakfit information for all spectra in profile"""
#        print '\n', self.profile_name
#
#        poscounter = 0
#        for spectrum in self.spectra_list:
#            print '\n', spectrum.fname, \
#                    self.positions_microns[poscounter], 'microns'
#            poscounter += 1
#            
#            if spectrum.peakpos is None:
#                check = spectrum.get_peakfit()
#                if check is False:
#                    print 'trouble with getting peakfit'
#            
#            if spectrum.peakpos is not None:
#                for k in range(len(spectrum.peakpos)):
#                    print spectrum.peakpos[k], spectrum.peak_heights[k], \
#                          spectrum.peak_widths[k], spectrum.peak_areas[k]        
#
#
#    def print_peakfits_ave(self, printout=True):
#        """Computes, prints, and returns average values in profile for
#        peak areas, peak heights, and the summed average peak areas"""
#        # get peak information for all spectra
#        for spectrum in self.spectra_list:
#            if spectrum.peak_areas is None:
#                self.get_peakfit()
#
#        average_peakareas = np.average(self.peak_areas, axis=1)
#        average_peakheights = np.average(self.peak_heights, axis=1)
#        total_area = np.sum(average_peakareas)
#        
#        max_peakareas = np.max(self.peak_areas, axis=1)
#        max_peakheights = np.max(self.peak_heights, axis=1)
#
#        if printout is True:
#            print '\n', self.profile_name
#            print 'peak positions (cm-1)\n', self.spectra_list[0].peakpos
#            print 'average peak areas (cm-2)\n', average_peakareas
#            print 'average peak heights (cm-1)\n', average_peakheights
#            print 'summed areas (cm-2)', total_area
#            print 'max peak areas (cm-2)\n', max_peakareas
#            print 'max peak heights (cm-1)\n', max_peakheights
#            
#        return average_peakareas, average_peakheights, total_area, \
#                max_peakareas, max_peakheights
#
#    def plot_peakfits(self, initial_too=False, legloc=1, top=1.2):
#        """Show fit peaks for all spectra in profile"""
#        if len(self.spectra_list) < 1:
#            self.make_spectra_list()
#            
#        for spectrum in self.spectra_list:
#            spectrum.plot_peakfit(legloc=legloc, top=top)
#        
#        if initial_too is True and self.initial_profile is not None:
#            for spectrum in self.initial_profile.spectra_list:
#                spectrum.plot_peakfit(legloc=legloc)
#                
#    def make_wholeblock(self, peakfit=True, show_plot=False, bulk=True):
#        """Take initial and final profiles and make self.wb_areas.
#        If peakfit=True, then get_peakfit and make peak-specific ratios
#        peak_wb_areas, peak_wb_heights, peak_wb_widths"""
#
#        if self.initial_profile is None:
#            print '\n profile name:', self.profile_name
#            print 'Need to specify an initial profile'
#            return False
#
#        # Bulk H whole-block
#        if bulk is True:
#            # get both initial and final area lists
#            init = self.initial_profile
#            fin = self
#            for prof in [init, fin]:
#                if prof.areas_list is None:
#                    check = prof.make_area_list(1, show_plot=False)
#                    if check is False:
#                        return False
#
#                if ((len(init.areas_list) == 2) and 
#                    (init.positions_microns[0] == init.positions_microns[1])):
#                    p = (0, np.mean(init.areas_list))
#                elif len(init.areas_list) > 1:
#                    p = np.polyfit(init.positions_microns, init.areas_list, 1)
#                else:
#                    p = (0, init.areas_list[0])
#            
#            init_line = np.polyval(p, self.positions_microns)
#            area_ratio = self.areas_list / init_line
#            self.wb_areas = area_ratio
#        if peakfit is False:
#            return area_ratio
#    
#        # Peak-specific whole-block values
#        if self.peak_areas is None:
#            self.get_peakfit()
#            if self.peak_areas is None:
#                print '\nCould not get peak areas'
#                return False
#            
#        if self.initial_profile.peak_areas is None:
#            self.initial_profile.get_peakfit()
#            if self.initial_profile.peak_areas is None:
#                print '\nCould not get initial peak areas'
#                return False
#            
#        ipos = self.initial_profile.positions_microns
#        iareas = self.initial_profile.peak_areas
#        iheights = self.initial_profile.peak_heights
#        iwidths = self.initial_profile.peak_widths
#        pos = self.positions_microns       
#        areas = self.peak_areas
#        heights = self.peak_heights
#        widths = self.peak_widths        
#        npeaks = len(iareas)
#        
#        wb_areas = []
#        wb_heights = []
#        wb_widths = []
#        for k in range(npeaks):           
#            pa = np.polyfit(ipos, iareas[k], 1)
#            ph = np.polyfit(ipos, iheights[k], 1)
#            pw = np.polyfit(ipos, iwidths[k], 1)
#            anormalizeto = np.polyval(pa, pos)
#            hnormalizeto = np.polyval(ph, pos)
#            wnormalizeto = np.polyval(pw, pos)
#            
#            wb_areas.append(areas[k] / anormalizeto)                
#            wb_heights.append(heights[k] / hnormalizeto)
#            wb_widths.append(widths[k] / wnormalizeto)
#
##            print 'peak #', k
##            print 'normalizing to areas', anormalizeto
##            print 'areas', areas[k]
##            print 'whole-block areas', areas[k] / anormalizeto
#
#        self.peak_wb_areas = wb_areas
#        self.peak_wb_heights = wb_heights
#        self.peak_wb_widths = wb_widths
#        return wb_areas
#
#    def get_peak_wb_areas(self, peak_idx=0, peakwn=None, 
#                          heights_instead=False):
#        """Returns peak-specific whole-block areas for the profile
#        AND peak wavenumber in cm-1 because usually the peak index is easier
#        to pass.
#        Defaults to the first peak in the peak position list"""
#        # Add check that argument is actually a profile
#        
#        if self.peakpos is None: 
#            print 'Getting peaks fit in matlab for', self.profile_name
#            self.get_peakfit()
#
#        # Getting peak from wavenumber if index not given
#        # Or peak wavenumber from index
#        if peak_idx is None:
#            if peakwn not in self.peakpos:
#                print 'There is no peak at', peakwn
#                print 'Peaks are at', self.peakpos
#                return
#            peak_idx = np.where(self.peakpos==peakwn)[0][0]
#            print 'peak at', peakwn, 'is index', peak_idx
#        else:
#            peakwn = self.peakpos[peak_idx]
#
#        # Make peak-specific whole-block areas
#        self.make_wholeblock(peakfit=True, show_plot=False, 
#                                 bulk=False)
#               
#        if heights_instead is False:
#            returnvals = self.peak_wb_areas[peak_idx]
#        else:
#            returnvals = self.peak_wb_heights[peak_idx]
#    
#        return returnvals, peakwn
#
#    def make_style_subtypes(self):
#        """Make direction-specific marker and line style dictionaries for 
#        profile based on profile's self.style_base"""
#        if self.style_base is None:
#            print 'Set base style (style_base) for profile first'
#            return False
#        self.style_x_marker =  dict(self.style_base.items() + styles.style_Dx.items())
#        self.style_y_marker =  dict(self.style_base.items() + styles.style_Dy.items())
#        self.style_z_marker =  dict(self.style_base.items() + styles.style_Dz.items())
#        self.style_x_line = make_line_style('x', self.style_base)
#        self.style_y_line = make_line_style('y', self.style_base)
#        self.style_z_line = make_line_style('z', self.style_base)
#        
#        if self.raypath == 'a':
#            self.style_x_marker.update({'marker' : styles.style_Rx['marker']})
#            self.style_y_marker.update({'marker' : styles.style_Rx['marker']})
#            self.style_z_marker.update({'marker' : styles.style_Rx['marker']})
#        elif self.raypath == 'b':
#            self.style_x_marker.update({'marker' : styles.style_Ry['marker']})
#            self.style_y_marker.update({'marker' : styles.style_Ry['marker']})
#            self.style_z_marker.update({'marker' : styles.style_Ry['marker']})
#        elif self.raypath == 'c':
#            self.style_x_marker.update({'marker' : styles.style_Rz['marker']})        
#            self.style_y_marker.update({'marker' : styles.style_Rz['marker']})
#            self.style_z_marker.update({'marker' : styles.style_Rz['marker']})
#        return
#
#    def choose_line_style(self):
#        """Returns line style with direction information"""
#        if self.style_base is None:
#            print 'Using default styles. Set profile style_base to change'
#            self.style_base = styles.style_profile            
#        if self.style_x_line is None:
#            self.make_style_subtypes()
#        if self.direction == 'a':
#            style_bestfitline = self.style_x_line
#        elif self.direction == 'b':
#            style_bestfitline = self.style_y_line
#        elif self.direction == 'c':
#            style_bestfitline = self.style_z_line
#        else:
#            style_bestfitline = {'linestyle' : '-'}
#        return style_bestfitline
#    
#    def choose_marker_style(self):
#        """Returns marker style with direction and ray path information"""
#        if self.style_base is None:
#            print 'Using default styles. Set profile style_base to change'
#            self.style_base = styles.style_profile
#
#        if self.style_x_marker is None:
#            self.make_style_subtypes()
#
#        if self.direction == 'a':
#            style = self.style_x_marker
#        elif self.direction == 'b':
#            style = self.style_y_marker
#        elif self.direction == 'c':
#            style = self.style_z_marker
#        else:
#            style = self.style_base
#
#        return style
#        
#    def plot_area_profile_outline(self, centered=True, peakwn=None,
#                                  set_size=(6.5, 4), top=1.2, 
#                                  wholeblock=False, heights_instead=False,
#                                  show_water_ppm=True):
#        """Set up area profile outline and style defaults. 
#        Default is for 0 to be the middle of the profile (centered=True)."""
#        if self.style_base is None:
#            self.style_base = styles.style_profile
#        self.make_style_subtypes()
#        
#        if self.length_microns is None:
#            leng = self.set_len()
#        else:
#            leng = self.length_microns
#
#        fig = plt.figure(figsize=set_size)
#        ax = SubplotHost(fig, 1,1,1)
#        fig.add_subplot(ax)
#
#        ax_ppm = ax.twin()
#        ax_ppm.axis["top"].major_ticklabels.set_visible(False)
#        
#        if show_water_ppm is True:
#            pass
#        else:
#            ax_ppm.axis["right"].major_ticklabels.set_visible(False)    
#        
#        ax.set_xlabel('Position ($\mu$m)')
#        
#        # Set y-label
#        if wholeblock is True:
#            if heights_instead is False:
#                ax.set_ylabel('Area/Area$_0$')
#            else:
#                ax.set_ylabel('Height/Height$_0$')            
#        else:
#            if heights_instead is False:
#                ax.set_ylabel('Area (cm$^{-2}$)')
#            else:
#                ax.set_ylabel('Height (cm$^{-1}$)')
#
#        if top is None:
#            if len(self.areas_list) > 1:
#                top = max(self.areas_list)+0.2*max(self.areas_list)
#            else:
#                top = 1.
#            
#        ax.set_ylim(0, top)
#
#        if centered is True:
#            ax.set_xlim(-leng/2.0, leng/2.0)
#        else:
#            ax.set_xlim(0, leng)
#                    
#        ax.grid()
#        return fig, ax, ax_ppm
#
#    def plot_area_profile(self, polyorder=1, centered=True, 
#                          heights_instead=False, figaxis=None, 
#                          bestfitline=False, style_bestfitline=None,
#                          show_FTIR=False, show_water_ppm=True,
#                          show_values=False, set_class=None,
#                          peakwn=None, peak_idx=None,
#                          style=styles.style_points, show_initial_areas=False,
#                          error_percent=0, wholeblock=False,
#                          label=None, initial_style=None,
#                          initial_label=None, phase='olivine',
#                          calibration='Bell', orientation_factor=3.):
#        """Plot area profile. Centered=True puts 0 in the middle of the x-axis.
#        figaxis sets whether to create a new figure or plot on an existing axis.
#        bestfitline=True draws a best-fit line through the data.
#        Set peak=wavenumber for peak-specific profile"""
#        if wholeblock is True and self.initial_profile is None:
#            print 'Need to specify an initial profile'
#            return False, False
#        
#        if len(self.spectra_list) < 1:
#            self.make_spectra_list()
#        
#        # Check for or create positions and areas
#        if len(self.positions_microns) < 1:
#            print 'Need positions_microns for profile'
#            return
#
#        # Get list of areas
#        if wholeblock is False:
#            if peakwn is None and peak_idx is None:
#                # bulk hydrogen
#                self.get_baselines()
#                areas = self.make_area_list(polyorder, show_FTIR, 
#                                        set_class, peak=peakwn)
#            else:
#                # peak-specific
#                if self.peak_areas is None:
#                    self.get_peakfit()
#    
#                if peak_idx is None:
#                    peak_idx = np.where(self.peakpos==peakwn)[0][0]
#                    print 'peak at', peakwn, 'is index', peak_idx
#                else:
#                    peakwn = self.peakpos[peak_idx]
#
#            if heights_instead is True and peak_idx is not None:
#                areas = self.peak_heights[peak_idx]
#            elif peak_idx is not None:
#                areas = self.peak_areas[peak_idx]
#
#        # whole-block
#        else:
#            # bulk hydrogen
#            if peak_idx is None:
#                if self.wb_areas is None:
#                    self.make_wholeblock(peakfit=False, bulk=True)
#                areas = self.wb_areas
#            # peak-specific
#            else:
#                if self.peak_wb_areas is None:
#                    self.make_wholeblock(peakfit=True, bulk=False)
#    
#                if heights_instead is False:
#                    areas = self.peak_wb_areas[peak_idx]
#                else:
#                    areas = self.peak_wb_heights[peak_idx]
#                    
#        if areas is False:
#            return            
#
#        if np.shape(areas) != np.shape(self.positions_microns):
#            print 'Area and positions lists are not the same size!'
#            print 'area:', np.shape(self.areas_list)
#            print 'positions:', np.shape(self.positions_microns)
#            return
#
#        # Use new or old figure axes
#        if figaxis is None:
#            f, ax, ax_ppm = self.plot_area_profile_outline(centered, peakwn,
#                           heights_instead=heights_instead, wholeblock=wholeblock,
#                           show_water_ppm=show_water_ppm)
#        else:
#            ax = figaxis
#            show_water_ppm = False
#
#        # Set length
#        if self.length_microns is None:
#            leng = self.set_len()
#            if self.length_microns is None:
#                print 'Need some information about profile length'
#                return
#        else:
#            leng = self.length_microns
#
#        # Set up plotting styles
#        if style_bestfitline is None:
#            style_bestfitline = self.choose_line_style()
#        if style is None:
#            style = self.choose_marker_style()
#        if label is None:
#            style['label'] = self.profile_name
#        else:
#            style['label'] = label
#
#        # Plot best fit line beneath data points
#        if bestfitline is True:
#            if ((len(areas) == 2) and 
#                (self.positions_microns[0] == self.positions_microns[1])):
#                p = (0, np.mean(areas))
#            elif len(areas) > 1:
#                p = np.polyfit(self.positions_microns, areas, 1)
#            else:
#                p = (0, areas[0])
#
#            self.bestfitline_areas = p
#            x = np.linspace(0, leng, 100)
#            y = np.polyval(p, x)
#            if centered is True:
#                ax.plot(x-leng/2.0, y, **style_bestfitline)
#            else:
#                ax.plot(x, y, **style_bestfitline)
#
#        # Plot data
#        if centered is True:
#            x = np.array(self.positions_microns) - leng/2.0
#        else:
#            x = self.positions_microns            
#        
#        yerror = np.array(areas)*error_percent/100.
#        
#        if error_percent == 0:
#            ax.plot(x, areas, **style)
#        else:
#            ax.errorbar(x, areas, yerr=yerror, **style)
#            
#        # Plot initial profile areas
#        if show_initial_areas is True:
#            if initial_style is None:
#                initial_style = styles.style_initial_point
#            if initial_label is None:
#                initial_style['label'] = 'initial'
#            else:
#                initial_style['label'] = label
#            self.initial_profile.plot_area_profile(style=initial_style,
#                                                   figaxis=ax)            
#        ax.set_ylim(0, max(areas)+0.2*(max(areas)))
#
#        # Title
#        if peak_idx is None:
#            tit = 'Bulk hydrogen'
#        else:
#            peakwn = self.peakpos[peak_idx]
#            tit = ' '.join(('Peak at', str(peakwn) ,'/cm'))
#        ax.set_title(tit)
#
#        if show_water_ppm is True:
#            ax_ppm.set_ylabel(''.join(('ppm H2O in ', phase, ', ', calibration, 
#                                       ' calibration *', str(orientation_factor))))
#            parasite_tick_locations = np.linspace(ax.get_ylim()[0],
#                                                  ax.get_ylim()[1], 5)
#            abs_coeff = absorption_coefficients(phase=phase, 
#                                                calibration=calibration, 
#                                                peak_idx=peak_idx)
##            ppm_labels = parasite_tick_locations * abs_coeff * orientation_factor
#            ppm_labels = parasite_tick_locations * abs_coeff.n * orientation_factor
#
#            ax_ppm.set_yticks(parasite_tick_locations)
##            ax_ppm.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
#            np.set_printoptions(precision=1)
#            ax_ppm.set_yticklabels(['{:.1f}'.format(i) for i in ppm_labels])
##    parasite_tick_locations = 1e4/(celsius_labels + 273.15)
##    ax_celsius.set_xticks(parasite_tick_locations)
##    ax_celsius.set_xticklabels(celsius_labels)
##    fig.add_subplot(ax)
##    ax.axis["bottom"].set_label("10$^4$/Temperature (K$^{-1}$)")
##    ax.axis["left"].set_label("log$_{10}$diffusivity (m$^{2}$/s)")
##    ax_celsius.axis["top"].set_label("Temperature ($\degree$C)")
##    ax_celsius.axis["top"].label.set_visible(True)
##    ax_celsius.axis["right"].major_ticklabels.set_visible(False)
#
##            max_water = area2water(ax.get_ylim()[1]*3, phase=phase, 
##                                   calibration=calibration)
##            ax_ppm.axis["right"].major_ticklabels.set_visible(False) 
##            ax_ppm.axis["right"].set_label(''.join(("max ppm H$_2$O ~ ", 
##                                                    str(max_water))))
#         
#        if figaxis is None:
#            return f, ax
#        else:
#            return
#
#
##    def plot_wb_water(self, polyorder=1, centered=True, style=styles.style_profile):
##        """Plot area ratios and scale up with initial water"""
##        if self.sample.initial_water is None:
##            print 'Set sample initial water content please.'
##            return
##        a, w = self.make_3DWB_water_list(polyorder=1)
##        
##        fig, ax, leng = self.plot_area_profile_outline(centered)
##        top = max(a) + 0.2*max(a)
##        ax.set_ylim(0, top)
##        ax.set_ylabel('Final area / Initial area')
##        
##        if centered is True:
##            ax.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
##            ax.plot(self.positions_microns - leng/2.0, a, **style)           
##        else:
##            ax.plot([0, leng], [1, 1], **styles.style_1)
##            ax.plot(self.positions_microns, a, **style)
##        return fig, ax, leng
#
#    def y_data_picker(self, wholeblock, heights_instead, peak_idx=None):
#        """Pick out and return peak area or height data of interest."""
#        if wholeblock is True:
#            if peak_idx is None:
#                if self.wb_areas is None:
#                    self.make_wholeblock(peakfit=False, bulk=True)
#                y = self.wb_areas
#            # Peak-specific 
#            else:
#                if self.peak_wb_areas is None:
#                    self.make_wholeblock(peakfit=True, bulk=False)
#                    
#                if heights_instead is False:
#                    y = self.peak_wb_areas[peak_idx]
#                else:
#                    y = self.peak_wb_heights[peak_idx]
#        else:
#            if self.areas_list is None:
#                self.make_area_list()
#            y = self.areas_list
#        return y
#
#    def D_picker(self, wholeblock=True, heights_instead=False, peak_idx=None):
#        """Returns attribute name and diffusivity of interest. 
#        Consider using get_diffusivity() first"""
#        log10D_m2s = None
#        
#        if peak_idx is not None:
#            if heights_instead is True:
#                if wholeblock is False and self.D_height[peak_idx] != 0.0:
#                    log10D_m2s = self.D_height[peak_idx]
#                elif wholeblock is True and self.D_height_wb[peak_idx] != 0.0:
#                    log10D_m2s = self.D_height_wb[peak_idx]
#            else:
#                if wholeblock is False and self.D_peakarea[peak_idx] != 0.0:
#                    log10D_m2s = self.D_peakarea[peak_idx]
#                elif wholeblock is True and self.D_peakarea_wb[peak_idx] != 0.0:
#                    log10D_m2s = self.D_peakarea_wb[peak_idx]
#        elif wholeblock is False and self.D_area != 0.0:
#            log10D_m2s = self.D_area
#        elif self.D_area_wb != 0.0:
#            log10D_m2s = self.D_area_wb
#            
#        if log10D_m2s is None or log10D_m2s == 0.:
#            print '\nNeed to set profile bulk or peak diffusivity'
#            return None
#            
#        return log10D_m2s
#
#    def D_saver(self, D, error, wholeblock=True, heights_instead=False, 
#                peak_idx=None):
#        """Saves a diffusivity in the profile attribute of interest.
#        Consider using get_diffusivity() first"""
#        if peak_idx is not None:
#            if heights_instead is True:
#                if wholeblock is False:
#                    self.D_height[peak_idx] = D
#                    self.D_height_error[peak_idx] = error
#                elif wholeblock is True:
#                    self.D_height_wb[peak_idx] = D
#                    self.D_height_wb_error[peak_idx] = error
#            else:
#                if wholeblock is False:
#                    self.D_peakarea[peak_idx] = D
#                    self.D_peakarea_error[peak_idx] = error
#                elif wholeblock is True:
#                    self.D_peakarea_wb[peak_idx] = D
#                    self.D_peakarea_wb_error[peak_idx] = error
#        elif wholeblock is False:
#            self.D_area = D
#            self.D_area_error = error
#        else:
#            self.D_area_wb = D
#            self.D_area_wb_error = error
#
#    def scaling_factor_picker(self, maximum_value=None, wholeblock=True,
#                              heights_instead=False, peak_idx=None):
#        """Returns value to scale diffusion to"""
#        if maximum_value is not None:
#            scaling_factor = maximum_value
#            
#        elif wholeblock is True:
#            if peak_idx is None:
#                scaling_factor = self.maximum_wb_area
#            else:
#                if heights_instead is False:
#                    scaling_factor = self.peak_maximum_areas_wb[peak_idx]
#                else:
#                    scaling_factor = self.peak_maximum_heights_wb[peak_idx]
#        
#        else: 
#            # bulk water                
#            if peak_idx is None:
#                if len(self.areas_list) < 1:
#                    self.make_area_list()
#                if self.maximum_area is not None:
#                    scaling_factor = self.maximum_area
#                else:
#                    scaling_factor = np.max(self.areas_list)
#                    self.maximum_area = scaling_factor
#        
#            # peak-specific heights                
#            elif heights_instead is True: 
#                if self.peak_maximum_heights[peak_idx] != 0.0:
#                    scaling_factor = self.peak_maximum_heights[peak_idx]
#                else:                
#                    scaling_factor = np.max(self.peak_heights[peak_idx])
#                    self.peak_maximum_heights[peak_idx] = scaling_factor
#            # peak-specific areas
#            else:
#                if self.peak_maximum_areas[peak_idx] != 0.0:
#                    scaling_factor = self.peak_maximum_areas[peak_idx]
#                else:
#                    scaling_factor = np.max(self.peak_areas[peak_idx])
#                    self.peak_maximum_areas[peak_idx] = scaling_factor
#    
#        return scaling_factor
#
#    def plot_diffusion(self, log10D_m2s=None, time_seconds=None, 
#                       peak_idx=None, top=1.2, wholeblock=False,
#                       maximum_value=None, style=None, centered=False,
#                       heights_instead=False, points=200., symmetric=True,
#                       labelD=True, labelDx=None, labelDy=None,
#                       erf_or_sum='erf', label4legend=None,
#                       initial_unit_value=1., final_unit_value=0.):
#        """Plot diffusion curve with profile data."""
#        if wholeblock is True and self.initial_profile is None:
#            print 'Need to specify an initial profile'
#            return False, False
#        
#        if symmetric is False:
#            centered = False
#            
#        fig, ax = self.plot_area_profile(wholeblock=wholeblock, 
#                                         peak_idx=peak_idx, centered=centered,
#                                         heights_instead=heights_instead)
#                                                       
#        if self.length_microns is not None:
#            microns = self.length_microns
#        elif self.sample.thickness_microns is not None:
#            microns = self.set_len()
#        else:
#            print '\nNeed to set profile length directly or with a sample'
#            return False, False
#
#        # Get time
#        if time_seconds is not None:
#           pass
#        elif self.time_seconds is not None:
#            time_seconds = self.time_seconds
#        else:
#            print '\nNeed to set profile time_seconds attribute'
#            return False, False, False
#
#        # Get diffusivity - input directly or stored in an attribute
#        if log10D_m2s is None:
#           log10D_m2s = self.D_picker(wholeblock, heights_instead, 
#                                      peak_idx)   
#        if log10D_m2s is None:
#            print '\nNeed to set profile bulk or peak diffusivity'
#            return False, False, False
#
#        if peak_idx is not None and self.peakpos is None:
#            self.get_peakfit()
#
#        # Get scaling factor for diffusion curve
#        scaling_factor = self.scaling_factor_picker(maximum_value, 
#                                        wholeblock, heights_instead, peak_idx)
#                        
#        ax.plot(ax.get_xlim(), [scaling_factor, scaling_factor], '--k')
#        print '\nScaling to {:.2f} maximum value'.format(scaling_factor)
##        print 'You can set maximum_value to scale diffusion curve'
#       
#        # Setup and plot diffusion curves
#        if symmetric is True:
#            params = diffusion.params_setup1D(microns, log10D_m2s, time_seconds,
#                                              init=initial_unit_value, 
#                                              fin=final_unit_value)
#            x_diffusion, y_diffusion = diffusion.diffusion1D_params(params, 
#                                                                    points=points)        
#            if centered is False:
#                x_diffusion = x_diffusion + (self.length_microns/2.)
#        else:
#            params = diffusion.params_setup1D(microns*2, log10D_m2s, 
#                                              time_seconds,
#                                              init=initial_unit_value, 
#                                              fin=final_unit_value)
#            x_diffusion, y_diffusion = diffusion.diffusion1D_params(params, 
#                                                                points=points)        
#            x_diffusion = x_diffusion[int(points/2):]
#            y_diffusion = y_diffusion[int(points/2):]
#
#        ### FINALLY PLOTTING
#        ax.plot(x_diffusion, y_diffusion*scaling_factor, label=label4legend)
#            
#            
#        # label diffusivity on plot
#        if labelD is True:
#            if labelDx is None:
#                if centered is True:
#                    labelDx = -microns/2.
#                else:
#                    labelDx = self.length_microns * 0.05
#            if labelDy is None:
#                top = ax.get_ylim()[1]
#                labelDy = top-0.15*top
#            strD = "{:.2f}".format(log10D_m2s) 
#            
#            ax.text(labelDx, labelDy,
#                    ''.join(('  D=10$^{', strD, '}$ m$^2$/s')), 
#                    ha='left', fontsize=16, backgroundcolor='w') 
#                  
#        return fig, ax
#
#    def diffusion_residuals(self, time_seconds, log10D_m2s, wholeblock=True,
#                            heights_instead=False, peak_idx=None,
#                            initial_unit_value=1., final_unit_value=0.,
#                            show_plot=True, top=1.2, 
#                            maximum_value=None):
#        """Compare 1D diffusion curve with profile data.
#        Returns vector containing the residuals and 
#        and the variance = sqrt(sum of squares of the residuals)"""
#        if self.positions_microns is None:
#            print 'Need to set profile positions'
#            return
#
#        if peak_idx is not None:
#            if peak_idx is not None and self.peakpos is None:
#                self.get_peakfit()
#
#        y = self.y_data_picker(wholeblock, heights_instead, peak_idx)
#        scaling_factor = self.scaling_factor_picker(maximum_value, 
#                                        wholeblock, heights_instead, peak_idx)
#        
#        x = self.positions_microns
#        L = self.set_len()
#        t = time_seconds
#
#        # Need to work on this        
#        init = scaling_factor
#            
#        fin = final_unit_value
#        
#        params = diffusion.params_setup1D(L, log10D_m2s, t, init, fin, 
#                                          False, False, False)
#
#        xdif, model = diffusion.diffusion1D_params(params)
#        resid = diffusion.diffusion1D_params(params, x, y)
#        RSS = np.sum(resid**2)
##        plt.plot(x-L/2., y, '+k')
##        plt.plot(xdif, model, '-r')
#        
#        if show_plot is True:
#            f, ax = self.plot_diffusion(log10D_m2s, t, peak_idx, top,
#                                        wholeblock=wholeblock, 
#                                        heights_instead=heights_instead,
#                                        maximum_value=maximum_value)
#        return resid, RSS
#            
#    def fitD(self, time_seconds=None, points=200, 
#             initial_unit_value=1., vary_initial=True,
#             varyD=True, guess=-13., peak_idx=None, top=1.2, 
#             peakwn=None, wholeblock=True, centered=False,
#             show_plot=True, polyorder=1, heights_instead=False,
#             final_unit_value=0., vary_final=False, 
#             maximum_value=None, min_water=0., symmetric=False):
#        """Fits 1D diffusion curve to profile data."""
#        if self.time_seconds is None and time_seconds is None:
#            print 'Need time_seconds'
#            return
#        elif time_seconds is None:
#            time_seconds = self.time_seconds            
#
#        if self.length_microns is None:
#            print 'Need to set profile attribute length_microns'
#            return
#            
#        if self.positions_microns is None:
#            print 'Need to set profile positions'
#            return
#
#        if self.areas_list is None:
#            self.make_area_list()
#        
#        ### Choose y data to fit to ###
#        y = self.y_data_picker(wholeblock, heights_instead, peak_idx)
#        scaling_factor = max(y)
#        y = y / scaling_factor # scale down to unit, so y range between 0 and 1
#        
#        # Set up x data and other parameters
#        x = self.positions_microns
#        L = self.length_microns
#        t = time_seconds
#        D0 = guess
#        init = initial_unit_value
#        fin = final_unit_value
#        params = diffusion.params_setup1D(L, D0, t, init, fin, vD=varyD, 
#                                          vinit=vary_initial, vfin=vary_final)
#
#        dict_fitting = {'points' : points, 
##                        'symmetric' : symmetric,
##                        'centered' : centered
#                        }
#
##        # minimization
#        lmfit.minimize(diffusion.diffusion1D_params, params, args=(x, y), 
#                       kws=dict_fitting)
#        best_D = ufloat(params['log10D_m2s'].value, 
#                        params['log10D_m2s'].stderr)
#        best_init = ufloat(params['initial_unit_value'].value, 
#                         params['initial_unit_value'].stderr)   
#        print 'best-fit log10D m2/s', best_D
#        print 'best-fit initial    ', best_init*scaling_factor
#
##        # save results as attributes
##        # Use save_diffusivities to save to a file
##        if wholeblock is True:
##            if peak_idx is not None:
##                if heights_instead is False:
##                    self.D_peakarea_wb[peak_idx] = best_D.n
##                    self.D_peakarea_wb_error[peak_idx] = best_D.s
##                    self.peak_maximum_areas_wb[peak_idx] = best_init.n
##                else:
##                    self.D_height_wb[peak_idx] = best_D.n
##                    self.D_height_wb_error[peak_idx] = best_D.s
##                    self.peak_maximum_heights_wb[peak_idx] = best_init.n
##
##            else:
##                self.D_area_wb = best_D.n
##                self.D_area_wb_error = best_D.s
##                self.maximum_wb_area = best_init.n
##        
##        resid, RSS = self.diffusion_residuals(best_D.n, wholeblock, 
##                          heights_instead, peak_idx, best_init.n, 
##                          final_unit_value, show_plot=False, 
##                          maximum_value=best_init.n)
#
#        resid, RSS = self.diffusion_residuals(time_seconds=time_seconds, 
#                                              log10D_m2s=best_D.n,
#                                              wholeblock=wholeblock,
#                                              heights_instead=heights_instead,
#                                              peak_idx=peak_idx,
#                                              initial_unit_value=best_init.n,
#                                              final_unit_value=final_unit_value,
#                                              show_plot=False,
#                                              maximum_value=best_init.n)
#        # report results
#        print '\ntime in hours:', params['time_seconds'].value / 3600.
#        print 'initial unit value:', '{:.2f}'.format(best_init*scaling_factor)
#        print 'bestfit log10D in m2/s:', '{:.2f}'.format(best_D)
#        print 'residual sum of squares', '{:.2f}'.format(RSS)
#        if show_plot is True:
#            fig, ax = self.plot_diffusion(time_seconds=time_seconds,
#                                          log10D_m2s=best_D.n,
#                                          peak_idx=peak_idx, top=top, 
#                                          wholeblock=wholeblock,
#                                          centered=centered,
#                                          symmetric=symmetric,
#                                          heights_instead=heights_instead, 
#                                          maximum_value=best_init.n*scaling_factor,
#                                          final_unit_value=final_unit_value
#                                          )
#            return fig, ax
##        else:
##            return 1, 2
###        return best_D, best_init, RSS
#
#    def save_diffusivities(self, folder=None, 
#                           file_ending='-diffusivities.txt'):
#        """Save diffusivities for profile to a file"""
#        if folder is None:
#            folder = self.folder
#            
#        if self.short_name is None:
#            print 'Need profile short_name attribute to label file'
#            return
#            
#        if self.peakpos is None:
#            self.get_peakfit()
#
#        # Same format as output by print_diffusivities:
#        # Diffuvities from WB areas - area - WB heights - heights 
#        # Diffusivity - error - maximum value to scale up to for each
#
#        # Saving diffusivities in first column, errors in 2nd, 
#        # Bulk H in first row, then each peak-specific value
#        a = []
#        a.append([self.D_area_wb, self.D_area_wb_error, self.maximum_wb_area])
#        a.append([self.D_area,    self.D_area_error,    self.maximum_area])
#        
#        for k in range(len(self.peakpos)):
#            a.append([self.D_peakarea_wb[k], self.D_peakarea_wb_error[k],
#                     self.peak_maximum_areas_wb[k]])
#            a.append([self.D_peakarea[k], 
#                      self.D_peakarea_error[k], 
#                      self.peak_maximum_areas[k]])
#            a.append([self.D_height_wb[k], self.D_height_wb_error[k],
#                     self.peak_maximum_heights_wb[k]])
#            a.append([self.D_height[k],
#                     self.D_height_error[k], self.peak_maximum_heights[k]])
#        
#        workfile = ''.join((folder, self.short_name, file_ending))
#        with open(workfile, 'w') as diff_file:
#            diff_file.write(json.dumps(a))
#
#
#    def get_diffusivities(self, folder=None, 
#                           file_ending='-diffusivities.txt'):
#        """Get saved diffusivities for profile from a file"""
#        if folder is None:
#            folder = self.folder
#        
#        # Look into json format        
#        if self.short_name is None:
#            print 'Need profile short_name attribute to label file'
#            return
#            
#        if self.peakpos is None:
#            self.get_peakfit()
#
#        workfile = ''.join((folder, self.short_name, file_ending))
#        if os.path.isfile(workfile) is False:
#            print ' '
#            print self.fname      
#            print 'use save_diffusivities() to make -diffusivities.txt'
#            return
#
#        with open(workfile, 'r') as diff_file:
#            diffusivities_string = diff_file.read()
#
#        diffusivities = json.loads(diffusivities_string)
#
#        self.D_area_wb = diffusivities[0][0]
#        self.D_area_wb_error = diffusivities[0][1]
#        self.maximum_wb_area = diffusivities[0][2]
#        
#        self.D_area = diffusivities[1][0]
#        self.D_area_error = diffusivities[1][1]
#        self.maximum_area = diffusivities[1][2]
#        
#        if self.peakpos is None:
#            self.get_peakfit()
#        if self.peakpos is None:
#            print 'Having trouble getting peakfit info to grab diffusivities'
#            return            
#                    
#        npeaks = len(self.peakpos)
#        for k in range(npeaks):
#            self.D_peakarea_wb[k] = diffusivities[2+4*k][0]
#            self.D_peakarea_wb_error[k] = diffusivities[2+4*k][1]
#            self.peak_maximum_areas_wb[k] = diffusivities[2+4*k][2]
#            self.D_peakarea[k] = diffusivities[3+4*k][0]
#            self.D_peakarea_error[k] = diffusivities[3+4*k][1]
#            self.peak_maximum_areas[k] = diffusivities[3+4*k][2]
#            self.D_height_wb[k] = diffusivities[4+4*k][0]
#            self.D_height_wb_error[k] = diffusivities[4+4*k][1]
#            self.peak_maximum_heights_wb[k] = diffusivities[4+4*k][2]
#            self.D_height[k] = diffusivities[5+4*k][0]
#            self.D_height_error[k] = diffusivities[5+4*k][1]
#            self.peak_maximum_heights[k] = diffusivities[5+4*k][2]
#        
#        return diffusivities
#        
#    def print_diffusivities(self, show_on_screen=True):
#        """Print out all diffusivities, including peak-specific, in profile"""
#        if self.peakpos is None:
#            self.get_peakfit()
#
#        D_area_wb = ufloat(self.D_area_wb, self.D_area_wb_error)   
#        wbmax = '{:.2f}'.format(self.maximum_wb_area)
#
#        print '\n', self.profile_name
#        print ' Diffusivities as log10(D in m2/s), errors, then max value A/A0'
#        print 'bulkH', D_area_wb, wbmax
#
#        # peak-specific
#        for k in range(len(self.peakpos)):
#            D_area_wb = ufloat(self.D_peakarea_wb[k], 
#                               self.D_peakarea_wb_error[k])   
#            awb = '{:.2f}'.format(self.peak_maximum_areas_wb[k])            
#            Da_wb = '{:.2f}'.format(D_area_wb)            
#            str1 = ' '.join((Da_wb, str(awb)))
#            string = ' '.join((str(self.peakpos[k]), str1))
#            print string
#
##        ### Showing all ways to generate ###
##        # bulk H diffusivities and initial values
##        D_area_wb = ufloat(self.D_area_wb, self.D_area_wb_error)   
##        D_area = ufloat(self.D_area, self.D_area_error)        
##        wb = '{:.2f}'.format(D_area_wb)
##        a = '{:.2f}'.format(D_area)
##        wbmax = '{:.2f}'.format(self.maximum_wb_area)
##        if self.maximum_area is not None:
##            amax = '{:.1f}'.format(self.maximum_area)  
##        else:
##            amax = 'n/a'
##        st1 = ''.join((wb, '(', wbmax, ')'))
##        st2 = ''.join((a, '(', amax, ')'))
##        bulkstring = ''.join(('bulk H : ', st1, ';  ', st2, 
##                              ';  n/a;         n/a'))
##
##        if show_on_screen is True:
##            print '\n', self.profile_name
##            print ' Diffusivities as log10(D in m2/s) (max value)'
##            print '         wb areas;         areas;       wb heights;     heights'
##            print bulkstring
##
##        # peak-specific
##        peakstrings = []
##        for k in range(len(self.peakpos)):
##            D_area_wb = ufloat(self.D_peakarea_wb[k], 
##                               self.D_peakarea_wb_error[k])   
##            D_area = ufloat(self.D_peakarea[k], self.D_peakarea_error[k])   
##            D_h_wb = ufloat(self.D_height_wb[k], self.D_height_wb_error[k])   
##            D_h = ufloat(self.D_height[k], self.D_height_error[k])
##            a = '{:.2f}'.format(self.peak_maximum_areas[k])
##            awb = '{:.2f}'.format(self.peak_maximum_areas_wb[k])
##            h = '{:.2f}'.format(self.peak_maximum_heights[k])
##            hwb = '{:.2f}'.format(self.peak_maximum_heights_wb[k])
##            
##            Da = '{:.2f}'.format(D_area)
##            Da_wb = '{:.2f}'.format(D_area_wb)
##            Dh = '{:.2f}'.format(D_h)
##            Dh_wb = '{:.2f}'.format(D_h_wb)   
##            
##            st1 = ''.join((Da_wb, '(', str(awb), ')'))
##            st2 = ''.join((Da, '(', str(a), ')'))
##            st3 = ''.join((Dh_wb, '(', str(hwb), ')'))
##            st4 = ''.join((Dh, '(', str(h), ')'))
##            
##            string0 = str(self.peakpos[k])
##            stringD = ';  '.join((st1, st2, st3, st4))
##            string = ' '.join((string0, ':', stringD))
##            if show_on_screen is True:
##                print string
##            peakstrings.append(string)
##
#    def start_at_arbitrary(self, wn_matchup=3000, offset=0.):
#        """For each spectrum in profile, divide raw absorbance by thickness 
#        and set spectra abs_full_cm such that they overlap at the specified 
#        wavenumber wn_matchup with specified offset up from zero"""
#        for x in self.spectra_list:
#            # Divide by thickness if not done already
#            if x.abs_full_cm is None:
#                check = x.divide_by_thickness(folder=self.folder)
#                if check is False:
#                    return False
#            # print 'Setting to zero at wn_matchup'
#            index = (np.abs(x.wn_full - wn_matchup)).argmin()
#            abs_matched = (x.abs_full_cm - x.abs_full_cm[index]) + offset
#            x.abs_full_cm = abs_matched
#        return
#
##
## Other functions
##
#
#class TimeSeries(Profile):
#    def __init__(self, sample=None, fname_list=[], time_hours=[], folder='',
#                 thick_microns=None, style_base=styles.style_points):
#        self.sample = sample
#        self.thick_microns = thick_microns
#        self.fname_list = fname_list
#        self.times_hours = time_hours
#        self.style_base = style_base        
#        self.folder = folder
#        self.make_spectra_list()
#    
#    def plot_timeseries(self, y=None, peak_idx=None, tit=None, D_list=[], 
#                        thickness_microns=None, max_hours=None,
#                        style=None, idx_C0=0, figaxis=None): 
#        """Plot and return figure and axis of time-series data and
#        all diffusivities in D_list in m2/s
#        """
#        if max_hours is None:
#            max_hours = max(self.times_hours)
#            
#        if thickness_microns is None:
#            thickness_microns = self.thick_microns
#            
#        if style is None:
#            style = self.style_base
#    
#        if figaxis is None:
#            fig = plt.figure()
#            fig.set_size_inches(3.2, 3.2)
#            ax = fig.add_subplot(111)
#            fig.autofmt_xdate()
#        else:
#            ax = figaxis
#            
#        ax.set_xlabel('Time (hours)', fontsize=12)
#        ax.set_ylabel('Concentration/\nMaximum Concentration', fontsize=12)
#        ax.set_xlim(0, max_hours)
#
#        # curves for diffusivities in D_list    
#        for D in D_list:
#             t, cc = diffusion.diffusionThinSlab(log10D_m2s=D, 
#                                        thickness_microns=thickness_microns, 
#                                        max_time_hours=max_hours)
#             ax.plot(t, cc, '-k', linewidth=1)
#    
#        # plot area data
#        x = self.times_hours
#        if y is None:
#            if peak_idx is None:
#                if len(self.areas_list) < 1:
#                    self.make_area_list()
#                C = self.areas_list
#            else:
#                print 'not ready for peak-specific yet'
#            C0 = C[idx_C0]
#            y = np.array(C) / C0
#        
#        ax.plot(x, y, **style)
#        
#        if figaxis is None:
#            return fig, ax
#
#
#
#def subtract_2spectra(list2, wn_high=4000, wn_low=3000):
#    """Subtract spectrum 1 from spectrum 0 input as list between given
#    wavenumbers (defaults to 4000 and 3000 cm-1). Spectra do not need 
#    to be the same length or have the exact same wavenumbers."""
#    # initial checks
#    if len(list2) != 2:
#        print 'Takes a list of exactly 2 spectra'
#        print 'length:', len(list2)
#        return
#    for x in list2:
#        # Check they are spectra
#        if isinstance(x, Spectrum) is False:
#            print x, 'is not a Spectrum'
#            return
#        # Check for and if necessary make absorbance and wavenumber full range
#        if x.abs_full_cm is None:
#            check = x.start_at_zero()
#            if check is False:
#                return False
#        if x.wn_full is None:
#            x.make_baseline()
#
#    # Set up wavenumber list for just the main spectrum 0
#    x = list2[0]
#    index_lo = (np.abs(x.wn_full-wn_low)).argmin()
#    index_hi = (np.abs(x.wn_full-wn_high)).argmin()
#    
#    wn_upper_spectrum = x.wn_full[index_lo:index_hi]
#    abs_upper_spectrum = x.abs_full_cm[index_lo:index_hi]
#    abs_difference = np.zeros_like(wn_upper_spectrum)
#
#    idx = 0
#    for wn in wn_upper_spectrum:
#        # find index of nearest wavenumber in spectrum to be subtracted off
#        idx_subtract = (np.abs(list2[1].wn_full-wn)).argmin()
#        # subtract
##        abs_difference[idx] = (x.abs_full_cm[idx_upper_spectrum] - 
##                           list2[1].abs_full_cm[idx_subtract])
#        abs_difference[idx] = (abs_upper_spectrum[idx] - 
#                                list2[1].abs_full_cm[idx_subtract])
#        idx += 1
#
#    dx = wn_upper_spectrum[-1] - wn_upper_spectrum[0]
#    dy = np.mean(abs_difference)
#    area = dx * dy
#    return area
#
#    if len(list2[0].base_wn) != len(list2[1].base_wn):
#        print 'Length problem in subtract_spectra. To be dealt with.'
#        return
#    
#    # subtract
##    difference = np.zeros_like(list2[0].base_wn)
##    for inx in range(len(list2[0].base_wn)):
##        difference[inx] = list2[0].abs - list2[1].abs
##    return difference
#    
#
#def make_all_specta_lists(classname=Profile):
#    for obj in gc.get_objects():
#        if isinstance(obj, classname):
#            obj.make_spectra_list()
#    
#
##
## General plot setup
##
#
#def plotsetup_3x3minus2(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
#                        xtickgrid=250, ytickgrid=0.5):
#    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
#    fig = plt.figure()
#    fig.set_size_inches(6.5, 6.5)
##    fig.set_size_inches(3, 3)
#    gs = gridspec.GridSpec(3,3)
#    ax1 = plt.subplot(gs[0, 0])
#    ax2 = plt.subplot(gs[0, 1])
#    ax3 = plt.subplot(gs[0, 2])
#    ax5 = plt.subplot(gs[1, 1])
#    ax6 = plt.subplot(gs[1, 2])
#    ax8 = plt.subplot(gs[2, 1])
#    ax9 = plt.subplot(gs[2, 2])
#
#    axis_list = [ax1, ax2, ax3, ax5, ax6, ax8, ax9]
#    gs.update(bottom=0.14, left=0.15, right=0.95, top=0.95)    
#    xmajorLocator = MultipleLocator(xtickgrid)
#    ymajorLocator = MultipleLocator(ytickgrid)
#    for ax in axis_list:
#        ax.grid()
#        ax.set_xlim([xhi, xlo])
#        ax.set_ylim([ylo, yhi])
#        ax.xaxis.set_major_locator(xmajorLocator)
#        ax.yaxis.set_major_locator(ymajorLocator)
#        plt.setp(ax.get_xticklabels(), visible=False)
#        plt.setp(ax.get_yticklabels(), visible=False)
#        ax.get_xaxis().get_major_formatter().set_useOffset(False)
#    for k in range(len(axis_list)):
#        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
#                string.ascii_uppercase[k], fontweight='bold')
#    ax1.set_ylabel('Absorbance (cm$^{-1}$)\nPolarized\nBefore heating')
#    ax5.set_ylabel('Absorbance (cm$^{-1}$)\nUnpolarized\nBefore heating')
#    ax8.set_ylabel('Absorbance (cm$^{-1}$)\nUnpolarized\nAfter heating')
#    ax8.set_xlabel('Wavenumber (cm$^{-1}$)')
#    plt.setp(ax5.get_yticklabels(), visible=True)
#    plt.setp(ax9.get_xticklabels(), visible=True, rotation=45)
#    for ax in [ax1, ax8]:
#        plt.setp(ax.get_yticklabels(), visible=True)
#        plt.setp(ax.get_xticklabels(), visible=True, rotation=45)    
#    return axis_list
#
#def plotsetup_3x3(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
#                  xtickgrid=250, ytickgrid=0.5,
#                  fig_size_inches=(6.5, 6.5)):
#    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
#    fig = plt.figure()
#    fig.set_size_inches(fig_size_inches)
#    gs = gridspec.GridSpec(3,3)
#    ax1 = plt.subplot(gs[0, 0])
#    ax2 = plt.subplot(gs[0, 1])
#    ax3 = plt.subplot(gs[0, 2])
#    ax4 = plt.subplot(gs[1, 0])
#    ax5 = plt.subplot(gs[1, 1])
#    ax6 = plt.subplot(gs[1, 2])
#    ax7 = plt.subplot(gs[2, 0])
#    ax8 = plt.subplot(gs[2, 1])
#    ax9 = plt.subplot(gs[2, 2])
#    axis_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
#    gs.update(bottom=0.14, left=0.15, right=0.95, top=0.95)    
#    xmajorLocator = MultipleLocator(xtickgrid)
#    ymajorLocator = MultipleLocator(ytickgrid)
#    for ax in axis_list:
#        ax.grid()
#        ax.set_xlim([xhi, xlo])
#        ax.set_ylim([ylo, yhi])
#        ax.xaxis.set_major_locator(xmajorLocator)
#        ax.yaxis.set_major_locator(ymajorLocator)
#        plt.setp(ax.get_xticklabels(), visible=False)
#        plt.setp(ax.get_yticklabels(), visible=False)
#        ax.get_xaxis().get_major_formatter().set_useOffset(False)
#    for k in range(len(axis_list)):
#        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
#                string.ascii_uppercase[k], fontweight='bold')
#    ax8.set_xlabel('Wavenumber (cm$^{-2}$)')
#    for ax in [ax1, ax4, ax7]:
#        plt.setp(ax.get_yticklabels(), visible=True)
#    for ax in [ax7, ax8, ax9]:
#        plt.setp(ax.get_xticklabels(), visible=True, rotation=45)    
#
##    blue_line = mlines.Line2D([], [], label='Observed', 
##                              **sp.style_spectrum)
##    green_line = mlines.Line2D([], [], label='Fit peaks', **sp.style_fitpeak)
##    dashed_line = mlines.Line2D([], [], label='Sum of fit peaks', 
##                                **sp.style_summed)
##    leg = axis_list[7].legend(handles=[blue_line, green_line, dashed_line], ncol=3,
##            loc = 'upper center', bbox_to_anchor = (0.5, -0.5), fancybox=True)
#        
#    return axis_list
#
#def plotsetup_3stacked(yhi=2, ylo=0, xlo = 3200, xhi = 3800,
#                  xtickgrid=200, ytickgrid=0.5):
#    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
#    fig = plt.figure()
#    fig.set_size_inches(3, 6.5)
#    gs = gridspec.GridSpec(3, 1)
#    ax1 = plt.subplot(gs[0, 0])
#    ax2 = plt.subplot(gs[1, 0])
#    ax3 = plt.subplot(gs[2, 0])
#    axis_list = [ax1, ax2, ax3]
#    gs.update(bottom=0.14, left=0.35, right=0.95, top=0.95)    
#    xmajorLocator = MultipleLocator(xtickgrid)
#    ymajorLocator = MultipleLocator(ytickgrid)
#
#    for ax in axis_list:
##        ax.grid()
#        ax.set_xlim([xhi, xlo])
#        ax.set_ylim([ylo, yhi])
#        ax.xaxis.set_major_locator(xmajorLocator)
#        ax.yaxis.set_major_locator(ymajorLocator)
#        ax.get_xaxis().get_major_formatter().set_useOffset(False)
#        plt.setp(ax.get_xticklabels(), rotation=45)    
#
#     # Place label A., B. C.
##    for k in range(len(axis_list)):
##        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
##                string.ascii_uppercase[k], fontweight='bold')
#
#    ax3.set_xlabel('Wavenumber (cm$^{-1}$)')
#    for ax in [ax1, ax2]:
#        plt.setp(ax.get_xticklabels(), visible=False)
#    return axis_list
#
##%% Generate 3D whole-block area and water profiles
#def make_3DWB_area_profile(final_profile, 
#                           initial_profile=None, 
#                           initial_area_list=None, 
#                           initial_area_positions_microns=None,
#                           show_plot=True, top=1.2, fig_ax=None,
#                           peakwn=None):
#    """Take final 3D-wholeblock FTIR profile and returns
#    a profile of the ratio of the two (A/Ao). A/A0 is also saved in the 
#    profile attribute wb_areas.
#    Requires information about initial state, so either an initial 
#    profile (best) as the profile's initial_profile attribute,
#    or else an initial area list passed in with its position. 
#    Defaults to making a plot with max y-value 'top'.
#    Note initial_area_positions_microns is assumed to start at 0 and then 
#    gets shifted for the fit.
#    """
#    fin = final_profile
#    leng = fin.set_len()
#
#    # If whole-block areas are already made, use them. 
#    # Otherwise make them.
#    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
#        wb_areas = fin.wb_areas
#    else:
#        print fin.wb_areas
#        # initial checks
#        if len(fin.positions_microns) == 0:
#            print 'Need position information'
#            return False
#        if fin.length_microns is None:
#            check = fin.set_len()
#            if check is False:
#                print 'Need more info to set profile length'
#                return False
#        if fin.length_microns is None:
#            fin.set_len()
#        if len(fin.areas_list) == 0:
#            print 'making area list for profile'
#            fin.make_area_list(show_plot=False)
#    
#        # What to normalize to? Priority given to self.wb_initial_profile, then
#        # initial_profile passed in here, then initial area list passed in here.
#        if fin.initial_profile is not None:
#            
#            init = fin.initial_profile
#            if init.length_microns != fin.length_microns:
#                print 'initial and final lengths must be the same!'
#                return False
#            # Make sure area lists are populated
#            for profile in [init, fin]:
#                if len(profile.areas_list) == 0:
#                    print profile.profile_name
#                    print 'making area list for profile'
#                    profile.make_area_list(show_plot=False)
#            A0 = init.areas_list
#            positions0 = init.positions_microns           
#    
#        elif initial_profile is not None:
#            init = initial_profile
#            if isinstance(init, Profile) is False:
#                print 'initial_profile argument must be a pynams Profile.'
#                print 'Consider using initial_area_list and positions instead'
#                return False
#            # Make sure area lists are populated
#            for profile in [init, fin]:
#                if len(profile.areas_list) == 0:
#                    print 'making area list for profile'
#                    profile.make_area_list(show_plot=False)
#            A0 = init.areas_list
#            positions0 = init.positions_microns
#    
#        elif initial_area_list is not None:
#            if initial_area_positions_microns is None:
#                print 'Need initial_area_positions_microns for initial_area_list'
#                return False
#            A0 = initial_area_list
#            positions0 = initial_area_positions_microns
#            if len(fin.areas_list) == 0:
#                print 'making area list for final profile'
#                fin.make_area_list(show_plot=False)
#        else:
#            print 'Need some information about initial state'
#            return False
#    
#        # More initial checks
#        if len(fin.areas_list) != len(fin.positions_microns):
#            print 'area and position lists do not match'
#            print 'length areas_list:', len(fin.areas_list)
#            print 'length positions list:', len(fin.positions_microns)
#            return False    
#        if len(A0) < 1:        
#            print 'Nothing in initial area list'
#            return False
#        if len(positions0) < 1:
#            print 'Nothing in initial positions list'
#            return False
#        if len(A0) == 1:
#            print 'Using single point to generate initial line'
#            A0.extend([A0[0], A0[0]])
#            positions0.extend([0, fin.length_microns])
#            
#        # Use best-fit line through initial values to normalize final data
#        p = np.polyfit(positions0-(leng/2.), A0, 1)
#        
#        normalizeto = np.polyval(p, fin.areas_list)
#        wb_areas = fin.areas_list / normalizeto
#         
#        # Save whole-block areas as part of profile
#        fin.wb_areas = wb_areas    
#    
#    if show_plot is True:
#        if fig_ax is None:
#            f, ax = final_profile.plot_area_profile_outline()
#        else:
#            ax = fig_ax
#        ax.set_ylim(0, top)
#        ylabelstring = 'Final area / Initial area'
#        if peakwn is not None:
#            print 'NOT READY FOR PEAKFITTING YET'
##            extrabit = '\n for peak at ' + str(peakwn) + ' cm$^{-1}$'
##            ylabelstring = ylabelstring + extrabit
#        ax.set_ylabel(ylabelstring)
#        style = fin.choose_marker_style()
#        ax.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
#        ax.plot(fin.positions_microns - (leng/2.0), wb_areas, **style)
#        if fig_ax is None:
#            return wb_areas, f, ax
#        else:
#            return wb_areas
#    else:
#        return wb_areas
#
#
#def make_3DWB_water_profile(final_profile, water_ppmH2O_initial=None,
#                            initial_profile=None, 
#                            initial_area_list=None, 
#                            initial_area_positions_microns=None,
#                            show_plot=True, top=1.2, fig_ax=None):
#    """Take a profile and initial water content.
#    Returns the whole-block water concentration profile based on
#    the profile's attribute wb_areas. If wb_areas have not been made, 
#    some initial profile information and various options are passed
#    to make_3DWB_area_profile().
#    Default makes a plot showing A/Ao and water on parasite y-axis
#    """
#    fin = final_profile
#    init = initial_profile
#
#    # Set initial water
#    if water_ppmH2O_initial is not None:
#        w0 = water_ppmH2O_initial
#    else:
#        if fin.sample is not None:
#            if fin.sample.initial_water is not None:
#                w0  = fin.sample.initial_water
#        elif init is not None:
#            if init.sample is not None:
#                if init.sample.initial_water is not None:
#                    w0 = init.sample.initial_water
#        else:
#            print 'Need initial water content.'
#            return False
#    
#    # Set whole-block areas
#    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
#        wb_areas = fin.wb_areas
#    else:  
#        wb_areas = make_3DWB_area_profile(fin, initial_profile, 
#                                          initial_area_list, 
#                                          initial_area_positions_microns)
#    water = wb_areas * w0
#    if show_plot is True:
#        # Use a parasite y-axis to show water content
#        fig = plt.figure()
#        ax_areas = SubplotHost(fig, 1,1,1)
#        fig.add_subplot(ax_areas)
#        area_tick_marks = np.arange(0, 100, 0.2)
#        ax_areas.set_yticks(area_tick_marks)
#        ax_water = ax_areas.twin()
#        ax_water.set_yticks(area_tick_marks)
#        if isinstance(w0, uncertainties.Variable):
#            ax_water.set_yticklabels(area_tick_marks*w0.n)
#        else:
#            ax_water.set_yticklabels(area_tick_marks*w0)
#        ax_areas.axis["bottom"].set_label('Position ($\mu$m)')
#        ax_areas.axis["left"].set_label('Final area / Initial area')
#        ax_water.axis["right"].set_label('ppm H$_2$O')
#        ax_water.axis["top"].major_ticklabels.set_visible(False)
#        ax_water.axis["right"].major_ticklabels.set_visible(True)
#        ax_areas.grid()
#        ax_areas.set_ylim(0, 1.2)
#        if fin.length_microns is not None:
#            leng = fin.length_microns
#        else:
#            leng = fin.set_len()
#        ax_areas.set_xlim(-leng/2.0, leng/2.0)
#            
#        style = fin.choose_marker_style()
#        ax_areas.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
#        ax_areas.plot(fin.positions_microns-leng/2.0, wb_areas, **style)
#        return water, fig, ax_areas
#    else:
#        return water
#
#                
##%% Group profiles together as whole-block unit
#class WholeBlock():
#    def __init__(self,profiles=[], name='', peakfit=False, folder='',
#                 make_wb_areas=False, time_seconds=None, worksheetname=None,
#                 style_base = None, temperature_celsius=None,
#                 diffusivities_log10_m2s=None, get_baselines=False,
#                 diffusivity_errors=None, sample=None, D_area_wb=[-12., -12., -12.]):
#        self.profiles = profiles
#        self.folder = folder
#        self.name = name
#        self.time_seconds = time_seconds
#        self.worksheetname = worksheetname
#        self.style_base = style_base
#        self.temperature_celsius = temperature_celsius
#        self.sample = sample        
#        self.D_area_wb = D_area_wb
#        self.D_area_wb_error = [0., 0., 0.]
#        self.peak_diffusivities = []
#        self.peak_diffusivities_errors = []
#        
#        if len(self.profiles) > 0:
#            self.setupWB(peakfit=peakfit, make_wb_areas=make_wb_areas,
#                         get_baselines=get_baselines)
#                
#    def setupWB(self, peakfit=False, make_wb_areas=False, get_baselines=False):
#        """Sets up and checks WholeBlock instance
#        - Check that profiles list contains a list of three (3) profiles
#        - Generate list of initial profiles
#        - Generate list of profile directions
#        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
#        - Generate list of ray paths
#        - Verify three ray path directions are compatible with directions list
#        """
#        if len(self.profiles) != 3:
#            print 'For now, only a list of 3 profiles is allowed'
#            return False
#        d = []
#        r = []
#        ip = []
#        L = []
#        
#        for prof in self.profiles:
##            if isinstance(prof, Profile) is False:
##                print 'Only profile objects allowed in profile list!'
##                return False
#            d.append(prof.direction)
#            r.append(prof.raypath)
#            L.append(prof.set_len())
#            prof.time_seconds = self.time_seconds
#
#            if prof.sample is None and self.sample is not None:
#                prof.sample = self.sample
#            elif prof.sample is not None and self.sample is None:
#                self.sample = prof.sample
#            elif prof.sample != self.sample:
#                print 'Warning: profile sample does not match wb sample'
#                
#            if prof.positions_microns is None:
#                print 'profile', prof.profile_name
#                print 'Needs positions_microns attribute'
#
#            if prof.initial_profile is None:
##                print 'Need to set initial_profile attribute for each profile.'
##                print 'Assuming these are initial profiles...'
#                prof.initial_profile = prof
#            ip.append(prof.initial_profile)
#            
#            if make_wb_areas is True:
#                check = prof.make_wholeblock(peakfit=peakfit)
#                if check is False:
#                    return False
#
#            if prof.spectra_list is None:
#                prof.make_spectra_list()
#                
#            for spec in prof.spectra_list:
#                if spec.thick_microns is None and self.sample is not None:
#                    if prof.raypath == 'a':
#                        spec.thick_microns = np.mean(self.sample.length_a_microns)
#                    elif prof.raypath == 'b':
#                        spec.thick_microns = np.mean(self.sample.length_b_microns)
#                    elif prof.raypath == 'c':
#                        spec.thick_microns = np.mean(self.sample.length_c_microns)
#                    else:
#                        print 'Need wb.direction to assign thickness'
#
#        self.directions = d
#        self.raypaths = r
#        self.initial_profiles = ip 
#        self.lengths = L
#
#        if peakfit is True:
#            for prof in self.profiles:
#                prof.get_peakfit()
#
#        if get_baselines is True:        
#            self.get_baselines()
#
#        return True
#
#    def get_peakfit(self, peak_ending='-peakfit.CSV', 
#                    baseline_ending='-baseline.CSV'):
#        """Get peakfit information for all profiles"""
#        for prof in self.profiles:
#            prof.get_peakfit(peak_ending=peak_ending,
#                             baseline_ending=baseline_ending)
#
#    def print_baseline_limits(self, initial_too=True):
#        """Print out baseline wavenumber range for each profile"""
#        for prof in self.profiles:
#            print '\n', prof.profile_name
#            print prof.spectra_list[0].base_low_wn, 
#            print prof.spectra_list[0].base_high_wn
#            if initial_too is True:
#                print prof.initial_profile.profile_name
#                print prof.initial_profile.spectra_list[0].base_low_wn, 
#                print prof.initial_profile.spectra_list[0].base_high_wn
#
#    def get_baselines(self, initial_too=False, folder=None, delim=',', 
#                      baseline_ending='-baseline.CSV'):
#        """Get baselines for all spectra in whole block"""
#        if self.initial_profiles is None:
#            self.setupWB()
#        for prof in self.profiles:
#            for spectrum in prof.spectra_list:
#                spectrum.get_baseline(baseline_ending=baseline_ending,
#                                      folder=folder, delim=delim)
#        if initial_too is True:
#            for prof in self.initial_profiles:
#                for spectrum in prof.spectra_list:
#                    spectrum.get_baseline(baseline_ending=baseline_ending,
#                                          folder=folder, delim=delim)
#            
#    def plot_showbaselines(self):
#        """Plot baselines for all spectra in the whole block"""
#        for prof in self.profiles:
#            for spec in prof.spectra_list:
#                spec.plot_showbaseline()
#
#
#    def make_area_lists(self, polyorder=1, show_plot=False, set_class=None,
#                       shiftline=None, printout_area=False, peak=None):
#            """Make list of areas from all profiles"""
#            self.areas = []
#            for prof in self.profiles:
#                a = prof.make_area_list(polyorder, show_plot, set_class,
#                                        shiftline, printout_area, peak)
#                self.areas.append(a)
#           
#
#    def plot_spectra(self, profile_idx=None, show_baseline=True, 
#                     show_initial_ave=True,
#                     show_final_ave=True, plot_all=False, 
#                     initial_and_final_together=False, style=styles.style_spectrum, 
#                     stylei=styles.style_initial, wn=None):
#        """Plot all spectra in all or specified profile in whole-block"""
#        if profile_idx is None:
#            proflist = ([self.profiles[0]] + 
#                        [self.profiles[1]] + 
#                        [self.profiles[2]])
#        else:
#            proflist = [self.profiles[profile_idx]]
#            
#        for prof in proflist:
#            print prof.profile_name
#            prof.plot_spectra(show_baseline=show_baseline, 
#                  show_initial_ave=show_initial_ave, plot_all=plot_all,
#                  show_final_ave=show_final_ave,
#                  initial_and_final_together=initial_and_final_together,
#                  style=style, stylei=stylei, wn=None)
#
#    def plot_3panels_ave_spectra(self, peak_idx=None, peakwn=None, 
#                                 top=1., high=4000, low=3000,
#                                style=styles.style_spectrum_red, 
#                                stylei=styles.style_initial, show_raypaths=False,
#                                figsize=(6., 4), show_initial=True,
#                                legloc=5, label='Final', figax3=None):
#        """Three suplots showing average initial and final spectra in each
#        direction"""
#        if self.initial_profiles is None:
#            self.setupWB()
#
#        if figax3 is None:
#            f, ax3 = plt.subplots(1,3)
#            f.set_size_inches(figsize)
#            for idx, ax in enumerate(ax3[:3]):
#                ax.set_xlim(high, low)
#                ax.set_ylim(0., top)
#                if show_raypaths is True:
#                    raypathstring = ''.join(('profile || ', self.profiles[idx].direction, 
#                                             '\nray path || ', self.profiles[idx].raypath))
#                    ax.text(3500, top-top*0.22, raypathstring, 
#                            backgroundcolor='w', horizontalalignment='center')
#            
#        else:
#            ax3 = figax3
#            
#        for k in range(3):
#            prof = self.profiles[k]
#            avespec = prof.average_spectra()
#            ax3[k].plot(avespec.wn_full, avespec.abs_full_cm, label=label, 
#                        **style)
#        
#            if peak_idx is not None:
#                if self.profiles[k].peakpos is None:
#                    self.profiles[k].get_peakfit()
#                peakpos = self.profiles[k].peakpos
#                peakwn = peakpos[peak_idx]
#                ax3[k].plot([peakwn, peakwn], [0, top], color='r')
#
#            if show_initial is True:
#                iprof = self.initial_profiles[k]
#                if iprof is not None:
#                    initspec = iprof.average_spectra()
#                    ax3[k].plot(initspec.wn_full, initspec.abs_full_cm, 
#                            **stylei)
#
#        plt.setp(ax3[1].get_yticklabels(), visible=False)
#        plt.setp(ax3[2].get_yticklabels(), visible=False)
#        ax3[1].set_xlabel('wavenumber (cm$^{-1}$)')
#        ax3[0].set_ylabel('absorbance (cm$^{-1}$)')
#        
#
#        if show_initial is True:
#            ax3[1].legend(loc=legloc)
#        if figax3 is None:
#            plt.tight_layout()
#            plt.gcf().autofmt_xdate()
#            tit = ' '.join(('Averaged profiles for', self.name))
#            if peak_idx is not None:
#                tit = str.join(tit, ', peak at ', str(peakpos[peak_idx]), '/cm')
#            ax3[1].set_title(tit, zorder=100) # must come after tight_layout
#
#        plt.subplots_adjust(top=0.85, bottom=0.25)
#        if figax3 is None:
#            return f, ax3
#
#    def xy_picker(self, peak_idx=None, wholeblock=True, heights_instead=False,
#                  centered=True, unit='microns'):
#        """Picks out and returns appropriate x and y-data for 3D"""
#        positions = []
#        y = []
#            
#        for prof in self.profiles:
#            positions.append(prof.positions_microns)
#
#            # Bulk hydrogen            
#            if peak_idx is None:
#                
#                # whole-block
#                if wholeblock is True:
#                    if prof.wb_areas is None:
#                        print '\nMaking whole block area ratios'
#                        check = prof.make_wholeblock(peakfit=False, bulk=True,
#                                                     show_plot=False)
#                        if check is False:
#                            return False        
#                    y_to_add = prof.wb_areas
#                
#                # absolute areas
#                else:
#                    if prof.areas_list is None:
#                        prof.make_area_list()
#                    y_to_add = prof.areas_list
#
#
#            # Peak-specific                
#            else:
#                for idx, spec in enumerate(prof.spectra_list):
#                    if spec.peak_areas is None:
#                        spec.get_peakareas()
#                    prof.peak_heights[peak_idx][idx]=spec.peak_heights[peak_idx]
#                    prof.peak_areas[peak_idx][idx]=spec.peak_areas[peak_idx]
#
#                if wholeblock is True:
#                    peak_wb, peakwn = prof.get_peak_wb_areas(peak_idx, 
#                                           heights_instead=heights_instead)
#                    y_to_add = peak_wb
#
#                else:
#                    if heights_instead is False:
#                        y_to_add = prof.peak_areas[peak_idx]
#                    else:
#                        y_to_add = prof.peak_heights[peak_idx]
#            y.append(y_to_add)
#            
#        
#        if centered is True:
#            a = np.mean(self.profiles[0].sample.length_a_microns) / 2.
#            b = np.mean(self.profiles[1].sample.length_b_microns) / 2.
#            c = np.mean(self.profiles[2].sample.length_c_microns) / 2.
#            halflengths = [a, b, c]
#            for idx in xrange(3):
#                positions[idx] = positions[idx] - halflengths[idx]
#            
#        return positions, y
#
#    def plot_areas_3panels(self, peak_idx=None, fig_ax3=None, centered=True,
#                           top=None, wn=None, figsize=(6.5, 2.5), 
#                           show_spectra=True, percent_error=0., 
#                           xerror=0., yerror=None, pie=True,
#                           label4legend=[None, None, None],
#                           styles3=[styles.style_points]*3,
#                           use_area_profile_styles=True, unit='microns',
#                           heights_instead=False, wholeblock=True,
#                           show_line_at_1=True, get_saved_baseline=True,
#                           show_errorbars=True, peak_group=None):
#        """Plot whole-block ratio of Area/Initial Area (default) 
#        OR just areas (set wholeblock=False) on three panels"""
##        if get_saved_baseline is True:
##            self.get_baselines()
##
#        if wholeblock is True:
#            if ((self.directions is None) or (self.raypaths is None) or
#                (self.initial_profiles is None) or (self.lengths is None)):
#                    check = self.setupWB(make_wb_areas=False, peakfit=False)
#                    if check is False:
#                        print 'Problem setting up whole block in setupWB'                    
#                        return False
#        
#        if peak_idx is not None:
#            for prof in self.profiles:
#                for spec in prof.spectra_list:
#                    if spec.peak_areas is None:
#                        spec.get_peakareas()    
#
#        peakpos = self.profiles[0].peakpos
#
#        # concatenate positions and areas across three profiles to 
#        # send in to the plotting function
#        if peak_group is not None:
#            positions, y_placeholder = self.xy_picker(peak_idx=peak_idx, 
#                                                      wholeblock=wholeblock,
#                                                      heights_instead=heights_instead, 
#                                                      centered=centered, 
#                                                      unit=unit)
#            y = np.zeros_like(y_placeholder)
#            tit = 'Sum of peaks'
#
#            for peak_group_idx in peak_group:
#                positions, y_add = self.xy_picker(peak_idx=peak_group_idx, 
#                                                  wholeblock=wholeblock,
#                                                  heights_instead=heights_instead, 
#                                                  centered=centered, unit=unit)
#                y = y + y_add
#                tit = ' '.join((tit, str(peak_group_idx)))
#            
#        else:
#            
#            positions, y = self.xy_picker(peak_idx=peak_idx, wholeblock=wholeblock,
#                                          heights_instead=heights_instead, 
#                                          centered=centered, unit=unit)
#            
#            # Change title if peak-specific rather than bulk
#            if peak_idx is not None:
#                tit = ' '.join(('Peak at', str(peakpos[peak_idx]), '/cm'))
#            else:
#                tit = 'Bulk hydrogen'
#
#        if top is None:
#            z = [max(y[0]), max(y[1]), max(y[2])]
#            top = max(z) + 0.1*max(z)
#
#            
#        if use_area_profile_styles is True and None in styles3:
#            for k in range(3):
#                styles3[k] = self.profiles[k].choose_marker_style()
#                styles3[k]['markersize'] = 10
#        
#        if unit == 'microns':
#            lengths = self.lengths
#        elif unit == 'mm':
#            lengths = np.array(self.lengths) / 1000.
#        else:
#            print 'unit must be microns (default) or mm'
#            return
#                
#        # Sent positions and areas to plotting command
#        if fig_ax3 is not None:
#            styles.plot_3panels(positions, y, lengths, figaxis3=fig_ax3,
#                                styles3=styles3, top=top, wholeblock=wholeblock,
#                                show_line_at_1=show_line_at_1,
#                                heights_instead=heights_instead,
#                                label4legend=label4legend,
#                                use_errorbar=show_errorbars,
#                                yerror=yerror, unit=unit,
#                                percent_error=percent_error,
#                                xerror=xerror, centered=centered)
#            fig_ax3[1].set_title(tit)                                
#        else:
#            fig, ax = styles.plot_3panels(positions, y, lengths,
#                                          styles3=styles3, top=top, 
#                                          wholeblock=wholeblock,
#                                          show_line_at_1=show_line_at_1,
#                                          label4legend=label4legend,
#                                          heights_instead=heights_instead,
#                                          use_errorbar=show_errorbars,
#                                          percent_error=percent_error,
#                                          yerror=yerror, unit=unit,
#                                          xerror=xerror, centered=centered)
#            ax[1].set_title(tit)
#            fig.set_size_inches(figsize)
#            fig.autofmt_xdate()
#
#            if pie is True:
#                # add pie chart showing % of total height or area
#                if peak_idx is None:
#                    pass
#                else:
#                    ax_pie = fig.add_subplot(339)
#                    atot, htot = prof.get_area_total()
#                    
#                    if heights_instead is False:
#                        a = spec.peak_areas[peak_idx] / atot
#                        size = [a, 1. - a]
#                        tit = '% total area'
#            
#                    else:
#                        h = spec.peak_heights[peak_idx] / htot
#                        size = [h, 1. - h]
#                        tit = '% sum of heights'
#                
#                    colors = ['k', 'w']
#                    plt.pie(size, colors=colors, 
#                            startangle=90,
#                            radius=0.25, center=(0, 0), frame=False)
#                    ax_pie.axis('equal')
#                    ax_pie.set_title(tit)
#                
#            return fig, ax
#
#    def make_composite_peak(self, peak_idx_list):
#        """Make composite peaks for all spectra in whole block"""
#        for prof in self.profiles:
#            for spec in prof.spectra_list:
#                spec.make_composite_peak(peak_idx_list)
#            prof.get_peak_info()
#
#    def print_spectra_names(self, show_initials=True):
#        """Print out fnames of all spectra associated with the whole-block 
#        instance."""
#        if show_initials is True:
#            if self.initial_profiles is None:
#                self.setupWB()
#                
#            if self.initial_profiles is not None:
#                print '--Initial profiles--'
#                for prof in self.initial_profiles:
#                    print prof.profile_name
#                    spec_list = []
#                    for spectrum in prof.spectra_list:
#                        spec_list.append(spectrum.fname)
#                    print spec_list
#                    print ' '
#            else:
#                print 'No initial profiles given'
#
#        if self.profiles is not None:
#
#            print '--Final profiles--'
#            for prof in self.profiles:
#                print prof.profile_name
#                spec_list = []
#                for spectrum in prof.spectra_list:
#                    spec_list.append(spectrum.fname)
#                print spec_list
#                print ' '
#
#    def make_profile_list(self, initial_too=False):
#        """Return False or a list of profiles"""
#        if initial_too is True:
#            if self.initial_profiles is None:
#                self.setupWB()
#            if self.initial_profiles is None:
#                print 'No initial profiles'
#                profile_list = self.profiles
#            else:
#                profile_list = self.initial_profiles + self.profiles
#        else:
#            profile_list = self.profiles
#        return profile_list
#
#    def make_baselines(self, initial_too=False, linetype='line', 
#                       wn_high=3700., wn_low=3200., shiftline=None, 
#                       show_fit_values=False, show_plot=False,
#                       wn_mid=None): 
#        """Make spectra baselines for all spectra in all profiles in 
#        whole-block"""        
#        profile_list = self.make_profile_list(initial_too)        
#        for prof in profile_list:
#            prof.make_baselines(linetype=linetype, shiftline=shiftline, 
#                                show_fit_values=show_fit_values, 
#                                show_plot=show_plot, wn_mid=wn_mid,
#                                wn_high=wn_high, wn_low=wn_low) 
#
#    def save_baselines(self, initial_too=True):
#        """Make and save spectra baselines for all spectra."""
#        profile_list = self.make_profile_list(initial_too)
#        for prof in profile_list:
#            for spectrum in prof.spectra_list:
#                spectrum.save_baseline()
#
#    def matlab(self, initial_too=False):
#        """Print out a list of spectra names in a matlab-friendly way"""        
#        print '\nFor use in FTIR_peakfit_loop.m\n'        
#        
#        if initial_too is True:
#            if self.initial_profiles is None:
#                self.setupWB()
#            string = "{"
#            for prof in self.initial_profiles:
#                for spec in prof.spectra_list:
#                    stringname = spec.fname
#                    string = string + "'" + stringname + "' "
#            string = string + "};"
#            print string, '\n'
#        
#        string = "{"
#        for prof in self.profiles:
#            for spec in prof.spectra_list:
#                stringname = spec.fname
#                string = string + "'" + stringname + "' "
#        string = string + "}"
#        print string
#
#    def plot_areas(self, profile_index=None, peak_idx=None, peakwn=None, 
#                   show_initials=False, show_finals=True, show_legend=True,
#                   legloc=1, frame=False, top=None, bestfitlines=False,
#                   heights_instead=False, together=False):
#        """Plots profiles on one figure.
#        Set initial_instead_of_final to True to see initial areas.
#        Need to add legend and checks is wavenumber not in list."""
#
#        # Which profiles to plot
#        if show_initials is True:
#            if self.initial_profiles is None:
#                self.setupWB()
#            if self.initial_profiles is None:
#                print 'Need initial profile'
#                return
#            if self.initial_profiles is None:
#                self.setupWB(False, False)
#
#        # get wavenumber if only peak_idx is givin
#        if peak_idx is not None:
#            if profile_index is not None:
#                idx = profile_index
#            else:
#                idx = 0            
#            prof = self.profiles[idx]
#            if prof.peakpos is None:
#                prof.get_peakfit()
#            if peak_idx is None:
#                peak_idx = np.where(prof.peakpos==peakwn)[0][0]
#
#        f, ax = self.profiles[0].plot_area_profile_outline(peakwn=peakwn)
#
#        if profile_index is None:
#            if show_finals is True:
#                ai = self.profiles[0]
#                bi = self.profiles[1]
#                ci = self.profiles[2]
#                for prof in [ai, bi, ci]:
#                    prof.plot_area_profile(figaxis=ax, peakwn=peakwn, 
#                                           peak_idx=peak_idx,
#                                           bestfitline=bestfitlines,
#                                           heights_instead=heights_instead)
#
#            if show_initials is True:
#                ai = self.initial_profiles[0]
#                bi = self.initial_profiles[1]
#                ci = self.initial_profiles[2]            
#                ai.plot_area_profile(figaxis=ax, peakwn=peakwn, 
#                                     peak_idx=peak_idx,
#                                     heights_instead=heights_instead)
#                bi.plot_area_profile(figaxis=ax, peakwn=peakwn, 
#                                     peak_idx=peak_idx,
#                                     heights_instead=heights_instead)
#                ci.plot_area_profile(figaxis=ax, peakwn=peakwn, 
#                                     peak_idx=peak_idx,
#                                     heights_instead=heights_instead)
#        else:
#            if show_finals is True:
#                self.profiles[profile_index].plot_area_profile(
#                        figaxis=ax, peakwn=peakwn, peak_idx=peak_idx,
#                        heights_instead=heights_instead)
#            if show_initials is True:
#                self.initial_profiles[profile_index].plot_area_profile(
#                        figaxis=ax, peakwn=peakwn, peak_idx=peak_idx,
#                        heights_instead=heights_instead)
#
#        if show_legend is True:
#            leg_handle_list = []
#            descript = ['profile || a*', 'raypath || a*', 'profile || b', 
#                        'raypath || b',  'profile || c', 'raypath || c']
#            
#            bstylelinebase = {'marker' : 's', 'color' : 'black', 'alpha' : 0.5,
#                         'markersize' : 10, 'linestyle': 'none'}
#            bstyleline = [None, None, None, None, None, None]
#            bstyleline[0] = dict(bstylelinebase.items() + styles.style_Dx.items())
#            bstyleline[1] = dict(bstylelinebase.items() + styles.style_Rx.items())
#            bstyleline[2] = dict(bstylelinebase.items() + styles.style_Dy.items())
#            bstyleline[3] = dict(bstylelinebase.items() + styles.style_Ry.items())
#            bstyleline[4] = dict(bstylelinebase.items() + styles.style_Dz.items())
#            bstyleline[5] = dict(bstylelinebase.items() + styles.style_Rz.items()) 
#            
#            for k in range(6):
#                add_marker = mlines.Line2D([], [], label=descript[k], 
#                                           **bstyleline[k])
#                leg_handle_list.append(add_marker)
#            
#            # Shrink current axis's height by 10% on the bottom
#            box = ax.get_position()
#            ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                             box.width, box.height * 0.9])
#            
#            # Put a legend below current axis
#            ax.legend(handles=leg_handle_list,
#                      loc='upper center', bbox_to_anchor=(0.5, -0.3),
#                      fancybox=True, ncol=3)
#
#        # move y-axis upper limit to accomodate all data
#        all_areas = np.array([])
#        if top is not None:
#            ax.set_ylim(0, top)
#        else:
#            top = ax.get_ylim()[1]
#            if peak_idx is None:
#                for prof in [ai, bi, ci]:
#                    all_areas = np.append(all_areas, prof.make_area_list())
#            else:
#                for prof in [ai, bi, ci]:
#                    if heights_instead is False:
#                        all_areas = np.append(all_areas, 
#                                              prof.peak_areas[peak_idx])
#                    else:
#                        all_areas = np.append(all_areas, 
#                                              prof.peak_heights[peak_idx])
#                                              
#            maxtop = np.max(all_areas)
#            if maxtop > top:
#                top = maxtop + 0.15*maxtop
#        ax.set_ylim(0, top)
#                    
##        if peak_idx is None:
##            tit = 'Bulk hydrogen'
##        else:
##            tit = ' '.join(('Peak at',
##                            str(self.profiles[idx].peakpos[peak_idx]),'/cm'))
##        ax.text(0, top + 0.05*top, tit, horizontalalignment='center')
#                   
#        f.set_size_inches((6.5, 3.5))
#        plt.subplots_adjust(top=0.9, bottom=0.35)
#                   
#        return f, ax
#
#    def plot_peakfits(self, initial_too=False, profile_idx=None, legloc=1):
#        """Whole block: Plot peakfits for all spectra in all profiles"""
#        if profile_idx is None:
#            for prof in self.profiles:
#                prof.plot_peakfits(initial_too, legloc=legloc)
#        else: 
#            self.profiles[profile_idx].plot_peakfits(initial_too, 
#                                                    legloc=legloc)
#
#
#    def excelify(self, filename=None, exceldpi=150, top=1.0, 
#                 bestfitlines_on_areas=True, include_peakfits=True):
#        """Save all peak fit info AND peak figures to an excel spreadsheet """
#        # file setup
#        if filename is None and self.worksheetname is None:
#            print 'Need filename information here'
#            print 'or in worksheetname attribute of whole block instance'
#            return
#        if filename is None:
#            filename = ''.join((self.worksheetname, '.xlsx'))            
#        if self.profiles[0].peak_wb_areas is None:
#            self.setupWB(peakfit=True, make_wb_areas=True)
#        workbook = xlsxwriter.Workbook(filename)
#        worksheetname = self.worksheetname        
#        worksheet = workbook.add_worksheet(worksheetname)
#
#        # Column locations
#        col1 = 0
#        col2 = 4
#        col2pt5 = 8
#        col3 = 12
#        col4 = 28
#        col_heights = 20
#        col_wb_heights = 36
#        col5 = 44
#
#        # formats        
#        worksheet.set_column(col1, col1, width=13)
#        worksheet.set_column(col5, col5+50, width=15)
#        worksheet.set_column(col1+1, col2, width=11)
#        
#        wraptext = workbook.add_format()
#        wraptext.set_text_wrap()
#        wraptext.set_align('center')
#        wraptext.set_bold()
#        
#        boldtext = workbook.add_format()
#        boldtext.set_bold()
#
#        italictext = workbook.add_format()
#        italictext.set_italic()
#
#        one_digit_after_decimal = workbook.add_format()
#        one_digit_after_decimal.set_num_format('0.0')
#        one_digit_after_decimal.set_italic()
#
#        two_digits_after_decimal = workbook.add_format()
#        two_digits_after_decimal.set_num_format('0.00')
#
#        highlight = workbook.add_format()
#        highlight.set_bg_color('yellow')
#
#        peakpos = self.profiles[0].spectra_list[0].peakpos
#        
#        worksheet.write(0, 0, self.name)
#        worksheet.write(1, 0, 'peak position (/cm)', wraptext)
#        worksheet.write(1, 1, 'peak height (/cm)', wraptext)
#        worksheet.write(1, 2, 'peak width (/cm)', wraptext)
#        worksheet.write(1, 3, 'peak area (/cm2)', wraptext)
#        worksheet.write(1, 5, 'Baselines', wraptext)        
#        worksheet.write(1, 9, 'Baseline-subtracted', wraptext)
#        worksheet.set_column(9, 9, width=11)
#        
#        if include_peakfits is True:
#            # Bulk area and peak fits in 1st column for individual fits
#            row = 2        
#            pic_counter = 0
#            for prof in self.profiles:
#                row = row + 1
#                worksheet.write(row, col1, prof.profile_name, boldtext)
#                row = row + 1
#    
#                pos_idx = 0
#                for spec in prof.spectra_list:                
#                    row = row + 1
#                    worksheet.write(row, 0, spec.fname)                               
#                    worksheet.write(row, 1, prof.positions_microns[pos_idx],
#                                    one_digit_after_decimal)
#                    worksheet.write(row, 2, 'microns from edge', italictext)                   
#                    
#                    # Show baselines in 2nd column
#                    spec.get_baseline()
#                    f, ax = spec.plot_showbaseline()
#                    newfig = 'baseline{:d}.png'.format(pic_counter)
#                    pic_counter = pic_counter + 1
#                    f.savefig(newfig, dpi=exceldpi, format='png')
#                    worksheet.insert_image(row, col2, newfig, 
#                                           {'x_scale' : 0.5, 'y_scale' : 0.5})
#    
#                    # Peak fit figures in 2nd column
#                    f, ax = spec.plot_peakfit()
#                    newfig = 'peakfit{:d}.png'.format(pic_counter)
#                    pic_counter = pic_counter + 1
#                    f.savefig(newfig, dpi=exceldpi, format='png')
#                    worksheet.insert_image(row, col2pt5, newfig, 
#                                           {'x_scale' : 0.5, 'y_scale' : 0.5})
#    
#                    pos_idx = pos_idx + 1
#                    row = row + 1
#    
#                    if spec.peakpos is None:
#                        check = spec.get_peakfit()
#                        if check is False:
#                            print 'trouble with getting peakfit'
#                    sumarea = 0                    
#                    for k in range(len(spec.peakpos)):
#                        worksheet.write(row, 0, spec.peakpos[k])
#                        worksheet.write(row, 1, spec.peak_heights[k])
#                        worksheet.write(row, 2, spec.peak_widths[k])
#                        worksheet.write(row, 3, spec.peak_areas[k])                   
#                        
#                        sumarea = sumarea + spec.peak_areas[k]
#                        row = row + 1
#    
#                    worksheet.write(row, 0, 'Sum of peaks areas')
#                    worksheet.write(row, 3, sumarea, two_digits_after_decimal)
#                    row = row + 1
#    
#                    worksheet.write(row, 0, 'Observed total area')
#                    worksheet.write(row, 3, spec.area, two_digits_after_decimal)
#                    row = row + 1
#
#        # 3 panel averaged spectra in 3rd column
#        f, ax = self.plot_3panels_ave_spectra(top=top)
#        f.savefig('panel3.png', dpi=exceldpi, format='png')
#        worksheet.insert_image(0, col3, 'panel3.png', {'x_scale' : 0.5, 
#                               'y_scale' : 0.5})
#            
#        # area profiles, height profiles, whole-block profiles for both
#        # heights are for peak-specific only
#        worksheet.write(11, col3, 'Peak areas profiles', boldtext)
#        worksheet.write(11, col4, 'Whole-block area profiles', boldtext)
#        worksheet.write(25, col_heights, 'Peak height profiles', boldtext)
#        worksheet.write(25, col_wb_heights, 'Whole-block height profiles', 
#                        boldtext)
#        
#        errorstring1 = ', '.join(('+/- 2% errors in area', 
#                                  '+/- 50 microns errors in position'))
#        errorstring2 = ', '.join(('+/- 3% errors in whole-block area ratio', 
#                                  '+/- 50 microns errors in position'))
#        errorstring3 = ', '.join(('assuming +/- 2% errors in heights', 
#                                  '+/- 50 microns errors in position'))
#        errorstring4 = ', '.join(('+/- 3% errors in whole-block height ratio', 
#                                  '+/- 50 microns errors in position'))
#        worksheet.write(12, col3, errorstring1)
#        worksheet.write(12, col4, errorstring2)
#        worksheet.write(26, col_heights, errorstring3)
#        worksheet.write(26, col_wb_heights, errorstring4)
#
#        worksheet.write(13, col3, 'Errors typically plot in or near symbols')
#        worksheet.write(13, col4, 'Errors typically plot in or near symbols')
#        worksheet.write(27, col_heights, 
#                                'Errors typically plot in or near symbols')
#        worksheet.write(27, col_wb_heights, 
#                                'Errors typically plot in or near symbols')
#        
#        f, ax = self.plot_areas_3panels(wholeblock=False)
#        f.savefig('bulk_areas.png', dpi=exceldpi, format='png')
#        worksheet.insert_image(15, col3, 'bulk_areas.png', {'x_scale' : 0.5, 
#                               'y_scale' : 0.5})
#        
#        f, ax = self.plot_areas_3panels(wholeblock=True)
#        f.savefig('wbbulk.png', dpi=exceldpi, format='png')
#        worksheet.insert_image(15, col4, 'wbbulk.png', {'x_scale' : 0.5, 
#                               'y_scale' : 0.5})
#                       
#        for peak_idx in range(len(self.profiles[0].spectra_list[0].peakpos)):
#            f, ax = self.plot_areas_3panels(wholeblock=False, 
##                                            bestfitlines=bestfitlines_on_areas, 
#                                            peak_idx=peak_idx)
#            newfig = 'area{:d}.png'.format(peak_idx)
#            f.savefig(newfig, dpi=exceldpi, format='png')
#            worksheet.insert_image(29+(14*peak_idx), col3, newfig, 
#                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            
#
#            f, ax = self.plot_areas_3panels(peak_idx=peak_idx, wholeblock=True)
#            newfig = 'wbareas{:d}.png'.format(peak_idx)
#            f.savefig(newfig, dpi=exceldpi, format='png')
#            worksheet.insert_image(29+(14*peak_idx), col4, newfig, 
#                                   {'x_scale' : 0.5, 'y_scale' : 0.5})
#                   
#            f, ax = self.plot_areas_3panels(wholeblock=False,
##                                            bestfitlines=bestfitlines_on_areas, 
#                                            peak_idx=peak_idx, 
#                                            heights_instead=True)
#            newfig = 'height{:d}.png'.format(peak_idx)
#            f.savefig(newfig, dpi=exceldpi, format='png')
#            worksheet.insert_image(29+(14*peak_idx), col_heights, newfig, 
#                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            
#
#            f, ax = self.plot_areas_3panels(peak_idx=peak_idx, 
#                                            heights_instead=True,
#                                            wholeblock=True)
#            newfig = 'wbheights{:d}.png'.format(peak_idx)
#            f.savefig(newfig, dpi=exceldpi, format='png')
#            worksheet.insert_image(29+(14*peak_idx), col_wb_heights, newfig, 
#                                   {'x_scale' : 0.5, 'y_scale' : 0.5})
#
#        # Initial profiles and baseline information in 4rd column
#        row = 1
#        worksheet.write(row, col4, 'Initial profiles', boldtext) 
#        worksheet.write(row+1, col4, self.initial_profiles[0].profile_name)
#        worksheet.write(row+2, col4, self.initial_profiles[1].profile_name)
#        worksheet.write(row+3, col4, self.initial_profiles[2].profile_name)   
#        
#        worksheet.write(row+5, col4, 'Baseline wavenumber ranges (/cm)', 
#                        boldtext)
#        for k in range(3):
#            prof = self.profiles[k]
#            iprof = self.initial_profiles[k]
#            spec = prof.spectra_list[0]
#            ispec = iprof.spectra_list[0]
#            worksheet.write(row+6+k, col4, 
#                            ''.join(('final || ', prof.direction)))
#            worksheet.write(row+6+k, col4+3, 
#                            ''.join(('initial || ',iprof.direction)))
#            
#            if spec.base_high_wn != ispec.base_high_wn:
#                worksheet.write(row+6+k, col4+1, spec.base_high_wn, highlight)
#                worksheet.write(row+6+k, col4+4, ispec.base_high_wn, highlight)
#            else:
#                worksheet.write(row+6+k, col4+1, spec.base_high_wn)
#                worksheet.write(row+6+k, col4+4, ispec.base_high_wn)
#            
#            if spec.base_low_wn != ispec.base_low_wn:
#                worksheet.write(row+6+k, col4+2, spec.base_low_wn, highlight)
#                worksheet.write(row+6+k, col4+5, ispec.base_low_wn, highlight)
#            else:
#                worksheet.write(row+6+k, col4+2, spec.base_low_wn)
#                worksheet.write(row+6+k, col4+5, ispec.base_low_wn)
#
#
#        ### Summary list at the end - positions, areas, whole-block areas
#        # labels
#        worksheet.write(1, col5, 'centered peak position (/cm)', wraptext)
#        worksheet.write(1, col5+1, 'bulk area (/cm2)', wraptext)
#        col = col5 + 2
#        for peak in peakpos:
#            label1 = ''.join((str(peak), ' area (/cm2)'))
#            worksheet.write(1, col, label1, wraptext)
#            col = col + 1
#        worksheet.write(1, col, 'bulk whole-block ratio (/cm2)', wraptext)
#        col = col + 1
#        for peak in peakpos:
#            label2 = ''.join((str(peak), ' whole-block area ratio'))
#            worksheet.write(1, col, label2, wraptext)
#            col = col + 1
#        for peak in peakpos:
#            label1 = ''.join((str(peak), ' height (/cm)'))
#            worksheet.write(1, col, label1, wraptext)
#            col = col + 1
#        for peak in peakpos:
#            label2 = ''.join((str(peak), ' whole-block height ratio'))
#            worksheet.write(1, col, label2, wraptext)
#            col = col + 1
#        worksheet.write(1, col, 'file label', wraptext)
#
#        # values filling in the rows
#        row = 2
#        for prof in self.profiles:
#            col = col5
#  
#            # label profile name with spaces on either side
#            worksheet.write(row, col, prof.profile_name)            
#            halflen = prof.set_len() / 2.            
#            pos_idx = 0
#            row = row + 1
#
#            for spec in prof.spectra_list:
#                col = col5 # Restart at first column each time
#                worksheet.write(row, col, 
#                                prof.positions_microns[pos_idx]-halflen,
#                                two_digits_after_decimal)
#                col = col + 1
#
#                # bulk and peak fit areas
#                worksheet.write(row, col, prof.areas_list[pos_idx],
#                                two_digits_after_decimal)
#
#                col = col + 1
#                for k in range(len(peakpos)):
#                    area = prof.peak_areas[k][pos_idx]
#                    worksheet.write(row, col, area, 
#                                    two_digits_after_decimal)
#                    col = col + 1
#
#                # whole-block bulk and peakfit ratios
#                worksheet.write(row, col, prof.wb_areas[pos_idx],
#                                two_digits_after_decimal)
#
#                col = col + 1
#                for k in range(len(peakpos)):
#                    wb = prof.peak_wb_areas[k][pos_idx]
#                    if np.isnan(wb):                        
#                        worksheet.write(row, col, 'nan')
#                    elif np.isinf(wb): 
#                        worksheet.write(row, col, 'inf')
#                    else:
#                        worksheet.write(row, col, wb, 
#                                    two_digits_after_decimal)
#                    col = col + 1
#
#                # peak heights
#                for k in range(len(peakpos)):
#                    height = prof.peak_heights[k][pos_idx]
#                    worksheet.write(row, col, height, 
#                                    two_digits_after_decimal)
#                    col = col + 1
#                # peak whole-block heights
#                for k in range(len(peakpos)):
#                    wbh = prof.peak_wb_heights[k][pos_idx]
#                    if np.isnan(wbh):  
#                        worksheet.write(row, col, 'nan')
#                    elif np.isinf(wbh): 
#                        worksheet.write(row, col, 'inf')
#                    else:
#                        worksheet.write(row, col, wbh, 
#                                    two_digits_after_decimal)
#                    col = col + 1
#
#                # filenames at the end of summary
#                worksheet.write(row, col, spec.fname)
#
#                row = row + 1
#                pos_idx = pos_idx + 1
#
#        workbook.close()        
#
#    def print_peakfits(self, initial_too=False, excelfriendly=True):
#        """Print out all peakfit information for each spectrum in 
#        each profile"""
#        if initial_too is True and self.initial_profiles is not None:
#            proflist = self.profiles + self.initial_profiles
#        else:
#            proflist = self.profiles
#
#        if excelfriendly is True:
#            print 'position height width area'
#            
#        for prof in proflist:
#            prof.print_peakfits()
#
#    def print_max_arearatio(self, peak_idx=None, heights_instead=False):
#        """ Prints out the maximum whole-block area ratio observed 
#        in any profile of the wholeblock for specified peak_idx"""
#        if peak_idx is None:
#            self.setupWB(False, True)
#        else:
#            self.get_peakfit()
#            self.setupWB(True, False)
#
#        a = []
#
#        for prof in self.profiles:
#            if peak_idx is None:
#                maxval = max(prof.wb_areas)
#            else:
#                prof.make_wholeblock(peakfit=True)
#                if heights_instead is False:
#                    maxval = max(prof.peak_wb_areas[peak_idx])
#                else:
#                    maxval = max(prof.peak_wb_heights[peak_idx])
#
##            if np.isnan(maxval) is False:
#            a.append(maxval)
##            else:
##                a.append(0)
#        
#        print '\n', self.name
#        print max(a)
#        
#        return max(a)
#
#    def print_peakfits_ave(self, printall=True, print_max=False, 
#                           print_aves=True):
#        """ Prints out and returns average peak areas, peak heights, and 
#        sum of the average peak areas summed over all profiles"""
#        areas = []
#        heights = []
#        totalareas = []
#        max_a = []
#        max_h = []
#        for prof in self.profiles:
#            a, h, totala, ma, mh = prof.print_peakfits_ave(printout=printall)
#            areas.append(a)
#            heights.append(h)
#            totalareas.append(totala)
#            max_a.append(ma)
#            max_h.append(mh)
#
#            if print_max is True:
##                print '\n', prof.profile_name
##                print 'max areas'
#                print ma
##                print 'max_heights'
#                print mh
#            
#        if print_aves is True:
#            asum = np.array(np.sum(areas, axis=0))
#            hsum = np.array(np.sum(heights, axis=0))
#            tasum = np.sum(totalareas)
#    
#            print '\n', self.name
#            print 'peak positions (cm-1)'
#            print self.profiles[0].spectra_list[0].peakpos
#            print '\naverage peak areas summed over all profiles (cm-2)'
#            print asum
#            print '\naverage peak heights summed over all profiles (cm-1)'
#            print hsum
#            print '\naverage total area summed over all profiles (cm-1)'
#            print tasum
#
#
##        return asum, hsum, tasum, masum, mhsum
#
#
#    def print_diffusivities(self, peak_idx=None, profile_idx=None,
#                            show_plot=False, top=1.5):
#        """Print diffusivities for each profile"""
#        if show_plot is True:
#            self.plot_3panels_ave_spectra(peak_idx=peak_idx, top=top)
#                     
#        if peak_idx is None and profile_idx is None:
#            for prof in self.profiles:
#                prof.get_diffusivities()
#                prof.print_diffusivities()
#            return
#        
#        if profile_idx is None:
#            for prof in self.profiles:
#                k = peak_idx
#                print '\n', prof.profile_name
#                print 'peak position and log10(diffusivity in m2/s)'
#                print 'bulk H :', prof.D_area_wb, '+/-', \
#                    prof.D_area_wb_error
#                print prof.peakpos[k], ':', prof.D_peakarea_wb[k],\
#                    '+/-', prof.peak_D_area_wb_error[k]
#            return
#            
#        prof = self.profiles[profile_idx]
#        print '\n', prof.profile_name
#        print 'peak positions, log10(D in m2/s), D errors'
#        print 'bulk H :  ', prof.D_area_wb, '+/-', \
#                prof.D_area_wb_error
#        print prof.peakpos
#        print prof.D_peakarea_wb
#        print prof.peak_D_area_wb_error
#
#    def save_diffusivities(self, folder=None, 
#                           file_ending='-diffusivities.txt'):
#        """Save diffusivities for all profiles in whole-block instance
#        to files"""
#        for prof in self.profiles:
#            prof.save_diffusivities(folder, file_ending)
#            
#    def get_diffusivities(self, folder=None, 
#                           file_ending='-diffusivities.txt'):
#        """Gets diffusivities for all profiles in whole-block instance
#        from previously saved files"""
#        if folder is None:
#            folder = self.folder
#        for prof in self.profiles:
#            prof.get_diffusivities(folder, file_ending)
#
#    def plot_diffusion(self, peak_idx=None, time_seconds=None, 
#                       diffusivities_log10D_m2s=None, 
#                       erf_or_sum='erf', figsize=(6.5, 2.5),
#                       show_plot=True, xaxis='centered',
#                       show_slice=False, 
#                       style=[None, None, None], 
#                       styles3points=[None, None, None],
#                       wb_or_3Dnpi='wb', label4legend=[None, None, None],
#                       fig_ax=None, points=50, top_spectra=1.0,
#                       top=1.2, numformat='{:.1f}', wholeblock=True, 
#                       heights_instead=False, init=1., centered=True,
#                       fin=0, approximation1D=False, labelD=True,
#                       show_errorbars=True, labelDy=None, 
#                       labelDx=[None, None, None]):
#        """For whole-block data. See above for profiles.
#        Applies 3-dimensionsal diffusion equations using equations in 
#        pynams.diffusion and plots them on 3 panels with whole-block data.
#        Requires lengths, time in seconds, and three diffusivities either
#        explicitly passed here or as attributes of the WholeBlock object.
#        Assuming whole-block diffusion (wb_or_3Dnpi='wb') but could also 
#        do 3D non-path-integrated ('npi')
#        """        
#        if self.lengths is None:
#            self.setupWB(peakfit=False, make_wb_areas=False)
#        if self.lengths is None:
#            print 'Need to setup self.lengths, which is in microns'
#            return False
#
#        if (wb_or_3Dnpi != 'npi') and (wb_or_3Dnpi != 'wb'):
#            print 'wb_or_3Dnpi only takes "wb" or "npi"'
#            return False
#
#        if self.directions is None:           
#            self.setupWB(peakfit=False, make_wb_areas=False)
#
#        if self.initial_profiles is None:
#            self.setupWB(peakfit=False, make_wb_areas=False)
#
#        if self.raypaths is None:
#            self.setupWB(peakfit=False, make_wb_areas=False)
#
#        if time_seconds is None:
#            if self.time_seconds is not None:
#                time_seconds = self.time_seconds
#            else:
#                print 'Need time information'
#                return False
#                
#        # Pick which diffusivities to use
#        D3 = None
#        if diffusivities_log10D_m2s is not None:
#            D3 = diffusivities_log10D_m2s
#        elif wb_or_3Dnpi == 'wb' and peak_idx is None:
#            D3 = self.D_area_wb
#        else:
#            D3 = []
#            for prof in self.profiles:
#                D = prof.D_picker(wholeblock, heights_instead, peak_idx)
#                D3.append(D)
#
#        if D3 is None or 0.0 in D3:
#            print 'D3:', D3
#            print '\nNeed diffusivities.'
#            print 'Input directly as diffusivities_log10D_m2s'
#            print 'or input bulk in profile.D_area_wb'
#            print 'or peak_diffusivities at specified peak_idx\n'
#            return False
#                
#        L3 = self.lengths
#        
#        # Make this either way because it gets returned for possible
#        # use in other functions
#        params = diffusion.params_setup3D(L3, D3, time_seconds, 
#                                          init, fin)
#
#        xdiff, ydiff = diffusion.diffusion3Dwb_params(params, 
#                                                      raypaths=self.raypaths, 
#                                                      erf_or_sum=erf_or_sum,
#                                                      show_plot=False)
#        if show_plot is False:
#            return params, xdiff, ydiff
#            
#        else:
#            if fig_ax is None:
#                fig, fig_ax = self.plot_areas_3panels(peak_idx=peak_idx, 
#                                                      top=top,
#                                          wholeblock=wholeblock, 
#                                          heights_instead=heights_instead,
#                                          show_line_at_1=False,
#                                          label4legend=label4legend,
#                                          styles3=styles3points,
#                                          figsize=figsize, centered=centered,
#                                          show_errorbars=show_errorbars)
#                                          
#                if fig is False:
#                    return False, False, False
#                    
#            else:
#                fig = None
#
#            if centered is True:
#                for idx_len in xrange(3):
#                    xdiff[idx_len] = xdiff[idx_len] - (self.lengths[idx_len]/2.)
#
#            styles.plot_3panels(xdiff, np.array(ydiff), show_line_at_1=True, 
#                                figaxis3=fig_ax, init=init, top=top,
#                                styles3=style, 
#                                label4legend=label4legend,
#                                centered=centered)
#
#            # label diffusivities
#            if labelD is True:
#                for k in xrange(3):
#                    if labelDx[k] is None:
#                        if centered is True:
#                            labelDx[k] = 0.
#                        else:
#                            labelDx[k] = self.lengths[k]/2.
#                    dlabel = ''.join(('logD ', str(numformat.format(D3[k])), 
#                                      ' m$^2$/s'))
#                    if labelDy is None:
#                        labelDy = top-top*0.12
#                    fig_ax[k].text(labelDx[k], labelDy, dlabel, 
#                                    horizontalalignment='center',
#                                    backgroundcolor='w')                              
#   
#            return params, fig, fig_ax
#       
#    def fitD(self, peak_idx=None, init=1., fin=0.,
#             guesses_log10D=[-13., -13., -13.], 
#             heights_instead=False, wholeblock=True,
#             vary_initials=False, vary_finals=False, 
#             vary_diffusivities=[True, True, True],
#             erf_or_sum='erf',
#             show_plot=True, wb_or_3Dnpi='wb', centered=True,
#             show_initial_guess=True, style_initial=None,
#             style_final={'color' : 'red'}, points=50, top=1.2):
#        """Forward modeling to determine diffusivities in three dimensions 
#        from whole-block data. 
#        """        
#        # x and y are the data that we will fit to, centered for fitting
#        x, y = self.xy_picker(peak_idx, wholeblock, heights_instead, 
#                              centered=True)
#                
#        # for processing results
#        bestD = [] 
#        D3 = []
#        e3 = []
#
#        # set up fitting parameters in the reqquired format
#        params = diffusion.params_setup3D(microns3=self.lengths, 
#                 log10D3=guesses_log10D, 
#                 time_seconds=self.time_seconds, 
#                 initial=init, final=fin,
#                 vinit=vary_initials, vfin=vary_finals,
#                 vD=vary_diffusivities)
#
#        # other keywords needed for forward model
#        dict_fitting = {'points' : points,
#                        'erf_or_sum' : erf_or_sum} 
#
#        if wb_or_3Dnpi == 'wb':
#            # need raypaths and don't plot twice
#            if self.raypaths is None:
#                self.setupWB()
#            dict_fitting['raypaths'] = self.raypaths
#            dict_fitting['show_plot'] = False
#            
#            # run the minimizer
#            lmfit.minimize(diffusion.diffusion3Dwb_params, 
#                           params, args=(x, y), 
#                           kws=dict_fitting)
#            resid = diffusion.diffusion3Dwb_params(params, x, y, 
#                                            raypaths=self.raypaths,
#                                            erf_or_sum=erf_or_sum,
#                                            show_plot=False)
#        elif wb_or_3Dnpi == 'npi':
#            print 'npi not working well right now, sorry'
#            lmfit.minimize(diffusion.diffusion3Dnpi_params, 
#                           params, args=(x, y), 
#                           kws=dict_fitting)
#     
#            resid = diffusion.diffusion3Dnpi(params, x, y)
#        else:
#            print 'wb_or_3Dnpi can only be wb or npi'
#            return            
#
#
#        # convert to ufloats because ufloats are fun
#        bestD.append(ufloat(params['log10Dx'].value, 
#                            params['log10Dx'].stderr))
#        bestD.append(ufloat(params['log10Dy'].value, 
#                            params['log10Dy'].stderr))
#        bestD.append(ufloat(params['log10Dz'].value, 
#                            params['log10Dz'].stderr))
#        bestinit = (ufloat(params['initial_unit_value'].value, 
#                             params['initial_unit_value'].stderr))
#
#        # Plot and print results
#        for k in range(3):
#            D3.append(bestD[k].n)
#            e3.append(bestD[k].s)
#
#        if show_plot is True:
#            if wb_or_3Dnpi == 'wb':
#                self.plot_diffusion(init=init, top=top, 
#                                    peak_idx=peak_idx,
#                                    diffusivities_log10D_m2s=D3,
#                                    heights_instead=heights_instead,
#                                    centered=centered)
#            else:
#                 print 'sorry, only plotting wholeblock right now'
#                                             
#        print '\ntime in hours:', params['time_seconds'].value / 3600.
#        print '\ninitial unit values:', bestinit
#        print '\nbestfit log10D in m2/s:'
#        for D in bestD:
#            print D
#        print 'residual sum of squares:', np.sum(np.array(resid)**2.)
#        print D3[0], e3[0], D3[1], e3[1], D3[2], e3[2]
#                             
#        # Store values in profile attributes        
#        for k in range(3):
#            self.profiles[k].D_saver(D3[k], e3[k], wholeblock, 
#                            heights_instead, peak_idx)
#        return bestD
#    
#    def invert(self, grid_xyz, symmetry_constraint=True, 
#               smoothness_constraint=True, rim_constraint=True, 
#               rim_value=None, weighting_factor_lambda=0.2, 
#               show_residuals_plot=True):
#        """Takes a list of three whole-block concentration profiles (either A/Ao 
#        or water ok but must be consistent for all three) in three orthogonal 
#        directions and list of three integers to indicate number of divisions
#        in each direction. Returns matrix of values in each grid cell. 
#        Default plot showing residuals for how well results match the whole-block
#        observations."""
#        pass
#
#def make_line_style(direction, style_marker):
#    """Take direction and marker style and return line style dictionary
#    that reflects the direction (x, y, z, or u for unoriented) with the 
#    color of the base style"""
#    if direction == 'x':
#        d = styles.style_Dx_line
#    if direction == 'y':
#        d = styles.style_Dy_line
#    if direction == 'z':
#        d = styles.style_Dz_line
#    if direction == 'u':
#        d = styles.style_unoriented_line        
#    d.update({'linewidth' : 2})
#    return d
#    
#
#
#def plot_areas_functional(wbs, filename, filefolder, savefolder, 
#                          peak_ending=None, 
#                          get_baselines=True, baseline_ending='-baseline.CSV', 
#                          h=True, a=True, hwb=True, awb=True, span=True,
#                          peak_group=None, peaks_to_plot=None, a_top=50., 
#                          wb_top=1.4, peak_top=10.):
#    """Plots and saves area profiles for all spectra in set of whole blocks wbs.
#    Use h, a, hwh, awb to set whether to include heights, areas, whole-block
#    versions of both.    
#    """
#    # get the appropriate baselines
#    if get_baselines is True:
#        for wb in wbs:
#            wb.get_baselines(folder=filefolder,
#                             baseline_ending=baseline_ending)
#
#    pp = PdfPages(savefolder + filename)
#
#    def ok_plot(wholeblock, heights, peak_idx=None):
#        if peak is None: 
#            peak_label = 'Bulk hydrogen'
#            if wholeblock is False:
#                top = a_top
#            else:
#                top = wb_top
#        else:
#            peak_label = ' '.join(('peak', str(peakpos[peak_idx]), 'cm-1'))
#            if heights is False:
#                if wholeblock is False:
#                    top = 10.
#                else:
#                    top = 1.4
#            else:
#                if wholeblock is False:
#                    top = 0.25
#                else:
#                    top = 1.4
#
#        if heights is True and peak == None:
#            pass
#        else:
#            fig, ax3 = wb.plot_areas_3panels(wholeblock=wholeblock,
#                                             peak_idx=peak_idx, 
#                                             heights_instead=heights,
#                                             top=top,
#                                             styles3=[styles.style_points]*3,
#                                             show_line_at_1=wholeblock,
#                                             pie=wholeblock)
#            ax3[1].set_title(' '.join((wb.name, '\n', peak_label)))
#            fig.set_size_inches(8, 5)
#            fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15)        
#
#            if span is True:
#                if wholeblock is True:                
#                    for ax in ax3:
#                        ax.axhspan(0.8, 1.2, facecolor='g', alpha=0.5)
#
#            fig.savefig(pp, format='pdf')
#        plt.clf()
#
##        # get the peakfits - BUT SEE CHANGES IN ABOVE SECTION
##        wb.get_peakfit(peak_ending=peak_ending, baseline_ending=baseline_ending)
##        peakpos = wb.profiles[0].spectra_list[0].peakpos       
##        if peaks_to_plot is None and peakpos is not None:
##            peaks_to_plot = peakpos
##        
##    # plot and save areas, heights, whole-block areas and whole-block heights
##    # if they are included
#    ToDo = [a, h, awb, hwb]
#    wholeblock = [False, False, True, True]
#    heights = [False, True, False, True]
#    for idx, showthis in enumerate(ToDo):
#        if showthis is True:
#            for wb in wbs:                
#                # bulk 
#                peak = None
#                ok_plot(wholeblock=wholeblock[idx],
#                        heights=heights[idx])
#                
#                # peaks
##                if peakpos is not None:
##                    for peak_idx, peak in enumerate(peakpos):
##                        if peak in peaks_to_plot:
##                            ok_plot(wholeblock=wholeblock[idx],
##                                    heights=heights[idx], peak_idx=peak_idx)
#    
#    pp.close()
#
