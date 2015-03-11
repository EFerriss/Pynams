# -*- coding: utf-8 -*-
"""
FTIR class definitions and functions for processing and plotting FTIR spectra

The main unit is a class called Spectrum, which includes attributes 
such as sample thickness and lists of wavenumbers and absorbances. 
A few defaults are set up that are geared toward H in nominally anhydrous
minerals (NAMs) such a plotting wavenumber range from 3000 to 4000 /cm and
creating baselines between 3200 and 3700 /cm. 

More detailed instructions coming soon...

Copyright (c) 2015 Elizabeth Ferriss

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This software was developed using Python 2.7.

"""
# import necessary python modules
import gc
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import os.path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import string as string

# Optional: Set default folder where to find and/or save files for use with 
# functions like save_spectrum and save_baseline.
default_folder = 'C:\\Users\\Ferriss\\Documents\\FTIR\\'

# plotting styles
style_base = {'color' : 'k',
              'linewidth' : 1,
              'linestyle' :'-'}
style_spectrum = {'color' : 'b',
                  'linewidth' : 3}
style_fitpeak = {'color' : 'g',
                 'linewidth' : 1}
style_summed = {'color' : 'orangered',
                'linewidth' : 2,
                'linestyle' : '--'}

# Define classes, attributes, and functions related to samples
class StudySample():
    mineral_name = None
    source = None
    initial_water = None


class ThinSlab(StudySample):
    sample_thick_microns = None
    thick_list = []
    twoA_length = None
    twoB_length = None


class Block(StudySample):
    sample_thick_microns = None
    twoA_list = []
    twoB_list = []
    twoC_list = []

def get_3thick(sample_name):
    """Average thickness measurements in 3 directions and return a list"""
    twoA = np.mean(sample_name.twoA_list)
    twoB = np.mean(sample_name.twoB_list)
    twoC = np.mean(sample_name.twoC_list)
    return [twoA, twoB, twoC]

def make_gaussian(pos, h, w, x=np.linspace(3000, 4000, 150)):
    """Make and return Gaussian curve over specified range"""
    y = h * np.e**(-((x-pos) / (0.6005615*w))**2)
    return y

def area2water(area_cm2, mineral='cpx'):
    """Takes area in cm-2, multiplies by absorption coefficient, 
    return water concentration in ppm H2O"""
    cpx_calib_Bell95 = 1.0 / ufloat(7.09, 0.32)
    if mineral == 'cpx':           
        w = cpx_calib_Bell95*area_cm2
    else:
        print 'only mineral=cpx supported right now'
        return
    return w
        
# Define classes and functions related to FTIR spectra        
class Spectrum():
    fname = None  # Used to set filename.CSV
    sample = None
    thick_microns = None
    # full range of measured wavenumber and absorbances
    wn_full = None
    abs_full_cm = None
    # default wavenumber range of interest
    wn_high = 4000
    wn_low = 3000
    wn = None
    abs_cm = None
    # other metadata with a few defaults
    experiment = None
    raypath = None
    polar = None
    position_microns_a = None
    position_microns_b = None
    position_microns_c = None
    instrument = None
    spot_size_x_microns = 100
    spot_size_y_microns = 100
    resolution_cm = 4
    nscans = 100
    date = None
    other_name = None
    filename = None  # generated automatically below
    # default baseline range
    base_low_wn = 3200
    base_high_wn = 3700
    base_wn = None
    base_abs = None
    # need these for quadratic baseline
    base_mid_wn = 3450
    base_mid_yshift = 0.04
    base_w_small = 0.02
    base_w_large = 0.02
    # calculated after baseline set in subtract_baseline
    abs_nobase_cm = None
    # peak fit information
    peakpos = None
    numPeaks = None
    heights = None
    widths = None
    fitpeakareas = None
    peakshape = None    
    
    def set_thick(self):
        """Set spectrum thick_microns from raypath and sample.thick_microns.
        Similar utility exists for class Profile"""
        if isinstance(self.sample.sample_thick_microns, float) is True:
            th = self.sample.sample_thick_microns
        elif isinstance(self.sample.sample_thick_microns, list) is True:
            if len(self.sample.sample_thick_microns) == 3:
                if self.raypath == 'a':
                    th = self.sample.sample_thick_microns[0]
                elif self.raypath == 'b':
                    th = self.sample.sample_thick_microns[1]
                elif self.raypath == 'c':
                    th = self.sample.sample_thick_microns[2]
                else:
                    print 'Need spectrum raypath a, b, or c'
                    return False
            else:
                print 'sample.sample_thick_microns must have 1 or 3 values'
                return False
        self.thick_microns = th
        return th

    def make_average_spectra(self, spectra_list):
        """Takes list of spectra and returns average absorbance (/cm)
        to the new specified spectrum (self)"""
        list_abs_to_average = []
        list_wn = []
        for sp in spectra_list:
            if sp.abs_full_cm is None:
                check = sp.divide_by_thickness()
                if check is False:
                    return False
            abs_to_append = sp.abs_full_cm
            list_abs_to_average.append(abs_to_append)
            list_wn.append(sp.wn_full)
        ave_abs = np.mean(list_abs_to_average, axis=0)
        waven = np.mean(list_wn, axis=0)
        self.wn_full = waven
        self.abs_full_cm = ave_abs
    
    def divide_by_thickness(self, delim=','):
        """Divide raw absorbance by thickness"""
        if self.thick_microns is None:
            check = self.set_thick()
            if check is False:
                return False
        if self.fname is None:
            print 'Need .fname to know what to call saved file'
            return False
        if self.filename is None:
            make_filenames()
        signal = np.loadtxt(default_folder + self.filename, delimiter=delim)
        # Make sure wavenumber is ascending
        if signal[0, 0] > signal[-1, 0]:
            print 'WARNING: need to reverse wavenumber progression'
            # from operator import itemgetter
            # sig = sorted(signal, key=itemgetter(1))
        # print 'Getting full range signal wn_full, abs_full_cm'
        self.wn_full = signal[:, 0]
        self.abs_raw = signal[:, 1]
        # Convert from numpy.float64 to regular python float
        # or else element-wise division doesn't work.
        if isinstance(self.thick_microns, float) is True:
            th = self.thick_microns
        else:
            th = np.asscalar(self.thick_microns)
        self.abs_full_cm = self.abs_raw * 1e4 / th
        return self.abs_full_cm

    def start_at_zero(self):
        """Divide raw absorbance by thickness and
        shift minimum to 0 within specified wavenumber range """
        if self.abs_full_cm is None:
            check = self.divide_by_thickness()
            if check is False:
                return False
        # print 'Setting to zero in range of interest wn, abs_cm'
        index_lo = (np.abs(self.wn_full-self.wn_low)).argmin()
        index_hi = (np.abs(self.wn_full-self.wn_high)).argmin()
        indices = range(index_lo, index_hi, 1)
        self.wn = self.wn_full[indices]
        self.abs_cm = self.abs_full_cm[indices] - min(self.abs_full_cm[indices])
        return self.abs_cm

    def make_baseline(self, line_order=1, shiftline=None, 
                      show_fit_values=False):
        """Make linear (1) or quadratic baseline (2) baseline 
        and return baseline absorption curve. Shiftline value determines
        how much quadratic deviates from linearity"""
        if self.abs_cm is None:
            check = self.start_at_zero()
            if check is False:
                return False
        
        if shiftline is not None:
            yshift = shiftline
        else:
            yshift = self.base_mid_yshift
            
        index_lo = (np.abs(self.wn-self.base_low_wn)).argmin()
        index_hi = (np.abs(self.wn-self.base_high_wn)).argmin()        
        self.base_wn = self.wn[index_lo:index_hi]
        
        x = np.array([self.wn[index_hi], self.wn[index_lo]])
        y = np.array([self.abs_cm[index_hi], self.abs_cm[index_lo]])
        p = np.polyfit(x, y, 1)
        
        if line_order == 1:
            base_abs = np.polyval(p, self.base_wn)
        elif line_order == 2:
            x = np.insert(x, 1, self.base_mid_wn)
            yadd = np.polyval(p, self.base_mid_wn) - yshift
            y = np.insert(y, 1, yadd)
            if show_fit_values == True:
                print 'fitting x values:', x
                print 'fitting y values:', y
            p2 = np.polyfit(x, y, 2)
            base_abs = np.polyval(p2, self.base_wn)
        else:
            print "polynomial order must be either 1 (linear)",
            "or 2 (quadratic)"
            return
        self.base_abs = base_abs
        return base_abs

    def subtract_baseline(self, polyorder=1, bline=None, shifter=None):
        """Make baseline and return baseline-subtracted absorbance.
        Preferred baseline used is (1) input bline, (2) self.base_abs, 
        (3) one made using polyorder, default is linear"""
#        base_abs = self.make_baseline(polyorder)
        if self.wn is None:
            self.start_at_zero()
        
        if bline is not None:
#            print 'using baseline directly input as bline='
            base_abs = bline
        elif shifter is not None:
            print 'making custom baseline with shiftline=shifter'
            base_abs = self.make_baseline(polyorder, shiftline=shifter)          
        elif self.base_abs is not None:
            print 'using self.base_abs as baseline'
            base_abs = self.base_abs
        else:
            print 'making baseline...'
            base_abs = self.make_baseline(polyorder)
        
        index_lo = (np.abs(self.wn-self.base_low_wn)).argmin()
        index_hi = (np.abs(self.wn-self.base_high_wn)).argmin()
            
        humps = self.abs_cm[index_lo:index_hi]
        abs_nobase_cm = humps - base_abs
        return abs_nobase_cm

    def save_spectrum(self, folder=default_folder, delim='\t',
                      file_ending='-per-cm.txt'):
        """Save entire spectrum divided by thickness to file with headings.
        1st column: wavenumber /cm; 2nd column: absorbance /cm.
        Formatted for upload to PULI"""
        if self.fname is None:
            print 'Need .fname to know what to call saved file'
            return
        if self.abs_full_cm is None:
            self.divide_by_thickness()
            
        abs_filename = folder + self.fname + file_ending
        print 'Saving:'
        print abs_filename
#        t = ['wavenumber (/cm)', 'absorbance (/cm)']
#        with open(abs_filename, 'w') as abs_file:
#            for item in t:
#                abs_file.write(item+delim)
#            abs_file.write('\n')
        a = np.transpose(np.vstack((self.wn_full, self.abs_full_cm)))
        with open(abs_filename, 'w') as abs_file:
            np.savetxt(abs_file, a, delimiter=delim)
        
    def save_baseline(self, polyorder=1, folder=default_folder, 
                 bline_file_ending='-baseline.CSV',
                 shiftline=None, basel=None, delim=',', showplots=False):
        """Save baseline with baseline-subtracted spectrum 
        (-baseline) to file. These can be retrieved using get_baseline(). 
        Use save_spectrum() to save full wavenumber and -per-cm absorbances.
        """
        if self.fname is None:
            print 'Need .fname to know what to call saved file'
            return
            
        if self.base_abs is None:
            print 'making baseline...'
            basel = self.make_baseline(polyorder, shiftline)
        else:
            print 'using existing baseline'
            basel = self.base_abs

        if showplots is True:
            self.plot_showbaseline(baseline=basel)
        self.abs_nobase_cm = self.subtract_baseline(polyorder, bline=basel)
        base_filename = folder + self.fname + bline_file_ending
        print 'Saving', self.fname + bline_file_ending
        t = ['wavenumber (/cm)', 'baseline value (/cm)', 
             'baseline-subtracted absorbance (/cm)']
        with open(base_filename, 'w') as base_file:
            for item in t:
                base_file.write(item+delim)
            base_file.write('\n')
        a = np.transpose(np.vstack((self.base_wn, self.base_abs,
                                    self.abs_nobase_cm)))
        with open(base_filename, 'a') as base_file:
            np.savetxt(base_file, a, delimiter=delim)

    def get_baseline(self, folder=default_folder, delim=',', 
                bline_file_ending='-baseline.CSV'):
        """Get baseline and -subtracted spectrum saved using save_baseline().
        Same data that goes into MATLAB FTIR_peakfit_afterpython.m"""
        filename = folder + self.fname + bline_file_ending
        if os.path.isfile(filename) is False:
            print ' '
            print self.fname            
            print 'use save_baseline() to make -baseline.CSV'
            return
        data = np.genfromtxt(filename, delimiter=',', dtype='float', 
                             skip_header=1)
        self.base_wn = data[:, 0]
        self.base_abs = data[:, 1]
        self.abs_nobase_cm = data[:, 2]
        return self.base_abs, self.abs_nobase_cm
        
    def get_3baselines(self, folder=default_folder, delim=',', 
                bline_file_ending='-3baselines.CSV'):
        """Returns block of baseline data saved by water_from_spectra()
        d[:,0] = baseline wavenumbers, d[:,1:4] = baselines, 
        d[:,4:7] = baseline-subtracted absorbances"""
        filename = folder + self.fname + bline_file_ending
        if os.path.isfile(filename) is False:
            print ' '
            print self.fname            
            print 'Run water_from_spectra() with savebaselines=True'
            return
        data = np.genfromtxt(filename, delimiter=',', dtype='float', 
                             skip_header=1)
        return data
#        print data(1,1)
#        self.base_wn = data[1:,0]
#        baseline1 = data[1:,1]
#        baseline2 = data[1:,2]
#        baseline3 = data[1:,3]
#        return baseline1, baseline2, baseline3
                            
                
    def get_peakfit(self, folder=default_folder, delim=',', 
                    file_ending='-peakfit.CSV'):
        """Get individual peaks from fname-peakfit.CSV generated
        in MATLAB using FTIR_peakfit_afterpython.m and savefit.m
        Return curves (assumes Gaussians) and summation"""
        filename = folder + self.fname + file_ending
        if os.path.isfile(filename) is True:
            previous_fit = np.loadtxt(filename, delimiter=delim)
            self.peakpos = previous_fit[:, 0]
            self.numPeaks = len(self.peakpos)
            self.heights = previous_fit[:, 1]
            self.widths = previous_fit[:, 2]
            self.fitpeakareas = previous_fit[:, 3]
            self.peakshape = previous_fit[:, 4]
        else:
            print ' '            
            print self.fname
            print 'Remember to savefit.m after FTIR_peakfit_afterpython.m'
            return False

        peakfitcurves = np.ones([len(self.peakpos), len(self.base_wn)])
        summed_spectrum = np.zeros_like(self.base_wn)
        for k in range(len(self.peakpos)):
            peakfitcurves[k] = make_gaussian(x=self.base_wn, pos=self.peakpos[k],
                                        h=self.heights[k], w=self.widths[k])
            summed_spectrum += peakfitcurves[k]            
        return peakfitcurves, summed_spectrum
    
    def plot_showpeakfit(self):
        """Show peakfitting results. For now assumes Gaussian curves.
        NOT TOTALLY DONE YET."""
        if self.peakpos is None:
            check = self.get_peakfit()
            if check is False:
                return
                
        if self.base_wn is None:
            self.make_baseline()
            

    def water(self, polyorder=1, mineral_name='cpx', 
              show_plot=True, show_water=True,
              folder=default_folder, delim=',', numformat='{:.0f}',
              shiftline=None):
        """Produce water estimate without error for a single FTIR spectrum. 
        Use water_from_spectra() for errors, lists, and 
        non-Bell calibrations."""
        
        # Get or make absorbance and wavenumber range for baseline
        if (self.base_abs is None) or (self.base_wn is None):
            print 'Making baseline...'
            abs_baseline = self.make_baseline(polyorder, shiftline)
        else:
            print ('Using previously fit baseline. ' + 
                   'Use make_baseline to change it.')
            abs_baseline = self.base_abs
        wn_baseline = self.base_wn
            
        abs_nobase_cm = self.subtract_baseline(bline=abs_baseline)

        if show_plot is True:
            fig, ax = self.plot_spectrum_outline()
            plt.plot(self.wn, self.abs_cm, **style_spectrum) 
            ax.plot(wn_baseline, abs_baseline, **style_base)
            yhigh = max(self.abs_cm) + 0.1*max(self.abs_cm)
            if min(abs_baseline) > 0:
                ylow = 0
            else:
                ylow = min(abs_baseline) + 0.1*min(abs_baseline)
            plt.ylim(ylow, yhigh)
            plt.show(fig)
            
        dx = self.base_high_wn - self.base_low_wn
        dy = np.mean(abs_nobase_cm)
        area = dx * dy
        w = area2water(area, mineral=mineral_name)
                        
        # output for each spectrum
        if show_water is True:
            print self.fname
            print 'area:', numformat.format(area), '/cm^2'
            print 'water:', numformat.format(w), 'ppm H2O'
        return area, w

    def plot_spectrum_outline(self, size_inches=(3, 3), shrinker=0.05):
        """Make standard figure outline for plotting FTIR spectra"""
        f, ax = plt.subplots(figsize=size_inches)
        ax.set_xlabel('Wavenumber (cm$^{-1})$')
        ax.set_ylabel('Absorbance (cm$^{-1})$')
        ax.set_title(self.fname)
        ax.set_xlim(self.wn_high, self.wn_low)
        if self.abs_cm is None:
            self.start_at_zero()
        ax.grid()
        pad = max(self.abs_cm)
        yhigh = pad + 0.1*pad
        ylow = min(self.abs_cm) - 0.1*pad
        ax.set_ylim(ylow, yhigh)
        plt.tight_layout
        
        box = ax.get_position()
        ax.set_position([box.x0 + box.width*shrinker, 
                         box.y0 + box.height*shrinker, 
                         box.width*(1.0-shrinker), 
                         box.height*(1.0-shrinker)])

        return f, ax

    def plot_spectrum(self):
        """Plot the raw spectrum divided by thickness"""
        fig, ax = self.plot_spectrum_outline()
        if self.wn is None:
            self.start_at_zero()
        ax.plot(self.wn, self.abs_cm, **style_spectrum)
        return fig, ax
    
    def plot_peakfit(self):
        """Plot peaks fit in MATLAB using peakfit.m"""
        # Take baseline-subtracted spectrum from saved file every time 
        # to avoid any possible funny business from playing with baselines
        gaussian, summed_spectrum = self.get_peakfit()
        if self.abs_nobase_cm is None:
            self.abs_nobase_cm = self.subtract_baseline(2)
        observed = self.abs_nobase_cm
        fig, ax = self.plot_spectrum_outline()
        ax.plot(self.base_wn, observed, **style_spectrum)        
        ax.plot(self.base_wn, summed_spectrum, **style_summed)
        ax.set_xlim(max(self.base_wn), min(self.base_wn))
        ax.set_ylim(0, max(summed_spectrum)+0.2*max(summed_spectrum))
        ax.grid()
        for k in range(len(self.peakpos)):
            ax.plot(self.base_wn, gaussian[k], **style_fitpeak)
        return fig
        
    def plot_showbaseline(self, polyorder=1, shiftline=None,
                          wn_baseline=None, abs_baseline=None):
        """Plot FTIR spectrum and show baseline. 
        Can pass in your own baseline"""
        fig, ax = self.plot_spectrum_outline()
        
        # get baseline wavenumber range
        if wn_baseline is None:
            if self.base_wn is not None:
                wn_baseline = self.base_wn
            else:
                print ('Need to pass in baseline wavenumber range too' +
                       'either here as wn_baseline or as spectrum.base_wn')
                return
                
        # Get or make absorbance for baseline
        if abs_baseline is None:          
            if self.base_abs is None:
                print 'Making baseline...'
                abs_baseline = self.make_baseline(polyorder, shiftline)
            else:
                print ('Using previously fit baseline. ' + 
                       'Use make_baseline to change it.')
                abs_baseline = self.base_abs
            
        ax.plot(self.wn, self.abs_cm, **style_spectrum)
        ax.plot(wn_baseline, abs_baseline, **style_base)        
        return fig, ax

    def plot_subtractbaseline(self, polyorder=1):
        """Make and plot baseline-subtracted spectrum"""
        abs_nobase_cm = self.subtract_baseline(polyorder)
        fig, ax = self.plot_spectrum_outline()
        pad = max(abs_nobase_cm)
        yhigh = pad + 0.1*pad
        ax.set_ylim(0, yhigh)
        ax.set_xlim(self.base_high_wn, self.base_low_wn)
        ax.plot(self.base_wn, abs_nobase_cm, **style_spectrum)
        return fig, ax


def make_filenames(classname=Spectrum, folder=default_folder, 
                   file_ending='.CSV'):
    """ Set filename attribute based on folder and fname attribute
    for all Spectra() with fname but no filename."""
    for obj in gc.get_objects():
        if isinstance(obj, classname):
            if obj.fname is not None and obj.filename is None:
                obj.filename = ''.join((obj.fname, file_ending))

def water_from_spectra(list3, mineral_name='cpx',
                       proper3=False, numformat='{:.0f}',
                       savebaselines=False, show_plots=True, 
                       bline_file_ending='-3baselines.CSV', 
                       folder=default_folder, delim=',', 
                       calibration='Bell'):
    """Produce water estimate from list of FTIR spectra; 
    Default calibration is Bell et al. 1995, ideally using 
    3 spectra that are polarized in orthogonal directions (proper3=True)
    Default (proper3=False) is to estimate area and total water from each spectrum.
    Then average and std the results."""
    if calibration == 'Bell':
        pass
    elif calibration == 'Paterson':
        print 'Paterson calibration!'
        return
    else:
        print 'Only Bell (default) or Paterson calibration so far'
        return
        
    if proper3 == True and len(list3) != 3:
        print ' '
        print 'For proper3=True, list should contain only 3 spectra'
        proper3 = False
    
    uarea_list = []
    uwater_list = []
    for spec in list3:
        if show_plots is True:
            fig, ax = spec.plot_spectrum_outline()
            plt.plot(spec.wn, spec.abs_cm, **style_spectrum) 
        
        dx = spec.base_high_wn - spec.base_low_wn

        main_yshift = spec.base_mid_yshift                        
        window_large = spec.base_w_large
        window_small = spec.base_w_small

        if spec.abs_cm is None:
            spec.start_at_zero()

        # Generate list of 3 possible areas under the curve
        area_list = np.array([])
        spec.make_baseline()
        baseline_list = np.ones([3, len(spec.base_wn)])
        pek_list = np.ones([3, len(spec.base_wn)])
        
        k = 0
        for bam in [window_small, -window_small, -window_large]:
            # Setting up 3 different baseline shapes
            main_yshift = main_yshift - bam
            base_abs = spec.make_baseline(2, shiftline=main_yshift)
            abs_nobase_cm = spec.subtract_baseline(bline=base_abs)
            if show_plots is True:
                ax.plot(spec.base_wn, base_abs, **style_base)
            dy = np.mean(abs_nobase_cm)
            area = dx * dy
            area_list = np.append(area_list, area)
            baseline_list[k,:] = base_abs
            pek_list[k,:] = abs_nobase_cm
            k+=1

        uarea = ufloat(np.mean(area_list), np.std(area_list))
        uarea_list.append(uarea)
        uwater = area2water(uarea, mineral=mineral_name)
        # Water estimate depends on if you are adding 3 areas 
        # or just multiplying by three and later averaging
        if proper3 is False:
            uwater = 3*uwater
        uwater_list.append(uwater)

        # I show() the figures explicitly in this first loop because
        # otherwise they show up inline *after* all the printed output
        if show_plots is True:
            yhigh = max(spec.abs_cm) + 0.1*max(spec.abs_cm)
            if min(base_abs) > 0:
                ylow = 0
            else:
                ylow = min(base_abs) + 0.1*min(base_abs)
            plt.ylim(ylow, yhigh)
            plt.show(fig)

        if savebaselines is True:
            base_filename = folder + spec.fname + bline_file_ending
            t = ['wavenumber (/cm)', 'upper baseline (/cm) (smaller area)', 
                 'main baseline (/cm)', 'lower baseline (/cm) (larger area)',
                 'absorbance - upper baseline (/cm)', 
                 'absorbance - main baseline (/cm)',
                 'absorbance - lower baseline (/cm)']
            with open(base_filename, 'w') as base_file:
                for item in t:
                    base_file.write(item+delim)
                base_file.write('\n')
            a = np.transpose(np.vstack((spec.base_wn, baseline_list, pek_list)))
            with open(base_filename, 'a') as base_file:
                np.savetxt(base_file, a, delimiter=delim)
            
    # output for each spectrum
    for x in range(len(list3)):
        print list3[x].fname
        print "area:", numformat.format(uarea_list[x]), "/cm^2"
        print "water:", numformat.format(uwater_list[x]), "ppm H2O"
        if savebaselines is True:
            print 'Saved baselines to', list3[x].fname+bline_file_ending
        print ' '
    # Final output
    if proper3 == True:
        a = np.sum(uarea_list)
        w = np.sum(uwater_list)
        print 'Sum of individual contributions'
    else:
        a = np.mean(uarea_list)
        w = np.mean(uwater_list)
        print 'Averaged individual (*3) estimates'
    print 'area:', numformat.format(a), '/cm^2'
    print 'water:', numformat.format(w), 'ppm H2O'
    print ' '
    return a, w  

def list_with_attribute(classname, attributename, attributevalue):
    """Gather all instances in specified class with a particular attribute"""
    my_list = []
    for obj in gc.get_objects():
        if isinstance(obj, classname):
            if getattr(obj, attributename) == attributevalue:
                my_list.append(obj.fname)
    return my_list

#
# Begin experimental section - have not done much with this yet
#

class InTheLab():
    sample_tested = None
    initial_condition = None
    temperature_celcius = None
    pressure_bars = None
    fO2_bars = None
    fO2_QFM = None
    time_minutes = None
    final_profiles = None

#
# Define classes and functions for working with profiles of spectra
#

class Profile():
    # Required
    fname_list = []
    sample = None
    raypath = None
    spectrumclass_string = None
    # Optional
    direction = None
    name_profile = None
    initial_profile = None
    spectrum_class_name = None
    # Made by make_spectra_list
    positions_microns = np.array([])
    spectra_list = []
    len_microns = None
    thick_microns = None
    # Made by make_water_list
    waters_list = []
    waters_errors = []
    areas_list = []

    def set_len(self):
        """Set profile.len_microns from profile.direction and 
        profile.sample.thick_microns""" 
        if self.sample is None:
            print 'Need sample'
            return
        if self.direction == 'a':
           self.len_microns = self.sample.sample_thick_microns[0]
        elif self.direction == 'b':
            self.len_microns = self.sample.sample_thick_microns[1]
        elif self.direction == 'c':
            self.len_microns = self.sample.sample_thick_microns[2]
        else:
            print 'Set direction of profile to a, b, or c to set len_microns\n'
        return self.len_microns

    def set_thick(self):
        """Set profile.thick_microns from profile.raypath and
        profile.sample.thick_microns"""
        if len(self.sample.sample_thick_microns) == 1:
            self.thick_microns = self.sample.sample_thick_microns
        elif len(self.sample.sample_thick_microns) == 3:
            if self.raypath == 'a':
               self.thick_microns = self.sample.sample_thick_microns[0]
            elif self.raypath == 'b':
                self.thick_microns = self.sample.sample_thick_microns[1]
            elif self.raypath == 'c':
                self.thick_microns = self.sample.sample_thick_microns[2]
            else:
                print 'Need raypath'
                return False
        else:
            print 'sample.sample_thick_microns must have 1 or 3 values'
            return False
        return self.thick_microns

    def make_spectra_list(self, class_from_module=None, 
                          my_module='my_spectra'):
        """Set profile length and generate spectra_list 
        with key attributes"""
        if len(self.spectra_list) > 0:
            print len(self.spectra_list), 'overwriting spectra list'
            self.spectra_list = []
        if len(self.fname_list) == 0:
            print 'Need fnames'
            return False
        if self.sample is None:
            print 'Need sample'
            return False
        if (self.raypath is not None) and (self.raypath == self.direction):
            print "raypath cannot be the same as profile direction"
            return False
        if self.thick_microns is None:
            check = self.set_thick()
            if check is False:
                return False
            
        if class_from_module is not None:
            classToUse = class_from_module
        elif self.spectrum_class_name is not None:
            classToUse = self.spectrum_class_name
#            print 'Using class', self.spectrum_class_name
        else:
            classToUse = Spectrum
#            print 'Using generic Spectrum() as class'       

        fspectra_list = []
        for x in self.fname_list:
            newspec = classToUse()
            newspec.fname = x
            newspec.sample = self.sample
            newspec.raypath = self.raypath
            newspec.thick_microns = self.thick_microns
            fspectra_list.append(newspec)
        self.spectra_list = fspectra_list
        
        return fspectra_list

    def make_area_list(self, polyorder=1, show_plot=True, show_water=True, 
                        set_class=None):
        """Make list of areas (no errors!) and simple water 
        concentration estimates for profile"""
        if len(self.spectra_list) < 1:
            check = self.make_spectra_list(class_from_module=set_class)
            if check is False:
                return False
            
        areas = []
        waters = []
        waterse = []
        
        for x in self.spectra_list:
            a, w = x.water(polyorder, show_plot, show_water)
            areas.append(a)
            waters.append(w.n)
            waterse.append(w.s)
            
        self.areas_list = areas
        self.waters_list = waters
        self.waters_errors = waterse
        
    def plot_area_profile(self, polyorder=1, show_bestfit=False, 
                          show_initial=False, show_FTIR=False, 
                          show_values=False, set_class=None):
        """Plot area profile"""
        if len(self.positions_microns) < 1:
            print 'Need positions_microns'
            return
        if len(self.areas_list) < 1:
            print 'making water list'
            check = self.make_area_list(polyorder, show_FTIR, 
                                        show_values, set_class)
            if check is False:
                return
            
        if np.shape(self.areas_list) != np.shape(self.positions_microns):
            print 'Area and positions lists are not the same size!'
            print 'area:', np.shape(self.areas_list)
            print 'positions:', np.shape(self.positions_microns)
            return

        f, ax = plt.subplots(1, 1)
        ax.set_xlabel('Position ($\mu$m)')
        ax.set_ylabel('Area (cm$^{-2}$)')
        if self.profile_name is not None:
            ax.set_title(self.profile_name)
        else:
            print 'Consider adding profile_name'

        ax.plot(self.positions_microns, self.areas_list, 's')
        ax.grid()
        if self.len_microns is None:
            leng = max(self.positions_microns)
        else:
            leng = self.len_microns
            
        ax.set_xlim(0, leng)
        ax.set_ylim(0, max(self.areas_list)+0.2*max(self.areas_list))
        
        if show_initial is True:
            prof_init = self.initial_profile
            if prof_init is None:
                print 'Set initial_profile'
            else:
                if len(prof_init.positions_microns) == 0:
                    print 'Need positions_microns for initial profile'
                if len(prof_init.areas_list) == 0:
                    print 'Run make_water_list for initial profile'
                else:
                    p = np.polyfit(prof_init.positions_microns, 
                                   prof_init.areas_list, 1)
                    x = np.linspace(0, leng, 100)
                    y = np.polyval(p, x)
                    ax.plot(x, y, '--r')
            
        if show_bestfit is True:
            p = np.polyfit(self.positions_microns, self.areas_list, 1)
            x = np.linspace(0, leng, 100)
            y = np.polyval(p, x)
            ax.plot(x, y, '--g')
        return f, ax

    def start_at_arbitrary(self, wn_matchup=3000, offset=0.):
        """For each spectrum in profile, divide raw absorbance by thickness 
        and set spectra abs_full_cm such that they overlap at the specified 
        wavenumber wn_matchup with specified offset up from zero"""
        for x in self.spectra_list:
            # Divide by thickness if not done already
            if x.abs_full_cm is None:
                check = x.divide_by_thickness()
                if check is False:
                    return False
            # print 'Setting to zero at wn_matchup'
            index = (np.abs(x.wn_full - wn_matchup)).argmin()
            abs_matched = (x.abs_full_cm - x.abs_full_cm[index]) + offset
            x.abs_full_cm = abs_matched
        return

#
# Other functions
#

def subtract_2spectra(list2, wn_high=4000, wn_low=3000):
    """Subtract spectrum 1 from spectrum 0 input as list between given
    wavenumbers (defaults to 4000 and 3000 cm-1). Spectra do not need 
    to be the same length or have the exact same wavenumbers."""
    # initial checks
    if len(list2) != 2:
        print 'Takes a list of exactly 2 spectra'
        print 'length:', len(list2)
        return
    for x in list2:
        # Check they are spectra
        if isinstance(x, Spectrum) is False:
            print x, 'is not a Spectrum'
            return
        # Check for and if necessary make absorbance and wavenumber full range
        if x.abs_full_cm is None:
            check = x.start_at_zero()
            if check is False:
                return False
        if x.wn_full is None:
            x.make_baseline()

    # Set up wavenumber list for just the main spectrum 0
    x = list2[0]
    index_lo = (np.abs(x.wn_full-wn_low)).argmin()
    index_hi = (np.abs(x.wn_full-wn_high)).argmin()
    
    wn_upper_spectrum = x.wn_full[index_lo:index_hi]
    abs_upper_spectrum = x.abs_full_cm[index_lo:index_hi]
    abs_difference = np.zeros_like(wn_upper_spectrum)

    idx = 0
    for wn in wn_upper_spectrum:
        # find index of nearest wavenumber in spectrum to be subtracted off
        idx_subtract = (np.abs(list2[1].wn_full-wn)).argmin()
        # subtract
#        abs_difference[idx] = (x.abs_full_cm[idx_upper_spectrum] - 
#                           list2[1].abs_full_cm[idx_subtract])
        abs_difference[idx] = (abs_upper_spectrum[idx] - 
                                list2[1].abs_full_cm[idx_subtract])
        idx += 1

    dx = wn_upper_spectrum[-1] - wn_upper_spectrum[0]
    dy = np.mean(abs_difference)
    area = dx * dy
    return area

    if len(list2[0].base_wn) != len(list2[1].base_wn):
        print 'Length problem in subtract_spectra. To be dealt with.'
        return
    
    # subtract
#    difference = np.zeros_like(list2[0].base_wn)
#    for inx in range(len(list2[0].base_wn)):
#        difference[inx] = list2[0].abs - list2[1].abs
#    return difference
    

def make_all_specta_lists(classname=Profile):
    for obj in gc.get_objects():
        if isinstance(obj, classname):
            obj.make_spectra_list()
    

#
# General plot setup
#

def plotsetup_3x3minus2(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
                        xtickgrid=250, ytickgrid=0.5):
    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
    fig = plt.figure()
    fig.set_size_inches(6.5, 6.5)
#    fig.set_size_inches(3, 3)
    gs = gridspec.GridSpec(3,3)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])
    ax5 = plt.subplot(gs[1, 1])
    ax6 = plt.subplot(gs[1, 2])
    ax8 = plt.subplot(gs[2, 1])
    ax9 = plt.subplot(gs[2, 2])
    axis_list = [ax1, ax2, ax3, ax5, ax6, ax8, ax9]
    gs.update(bottom=0.14, left=0.15, right=0.95, top=0.95)    
    xmajorLocator = MultipleLocator(xtickgrid)
    ymajorLocator = MultipleLocator(ytickgrid)
    for ax in axis_list:
        ax.grid()
        ax.set_xlim([xhi, xlo])
        ax.set_ylim([ylo, yhi])
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
    for k in range(len(axis_list)):
        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
                string.ascii_uppercase[k], fontweight='bold')
    ax1.set_ylabel('Absorbance (cm$^{-1}$)\nPolarized\nBefore heating')
    ax5.set_ylabel('Absorbance (cm$^{-1}$)\nUnpolarized\nBefore heating')
    ax8.set_ylabel('Absorbance (cm$^{-1}$)\nUnpolarized\nAfter heating')
    ax8.set_xlabel('Wavenumber (cm$^{-2}$)')
    plt.setp(ax5.get_yticklabels(), visible=True)
    plt.setp(ax9.get_xticklabels(), visible=True, rotation=45)
    for ax in [ax1, ax8]:
        plt.setp(ax.get_yticklabels(), visible=True)
        plt.setp(ax.get_xticklabels(), visible=True, rotation=45)    
    return axis_list

def plotsetup_3x3(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
                  xtickgrid=250, ytickgrid=0.5):
    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
    fig = plt.figure()
    fig.set_size_inches(6.5, 6.5)
    gs = gridspec.GridSpec(3,3)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])
    ax4 = plt.subplot(gs[1, 0])
    ax5 = plt.subplot(gs[1, 1])
    ax6 = plt.subplot(gs[1, 2])
    ax7 = plt.subplot(gs[2, 0])
    ax8 = plt.subplot(gs[2, 1])
    ax9 = plt.subplot(gs[2, 2])
    axis_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    gs.update(bottom=0.14, left=0.15, right=0.95, top=0.95)    
    xmajorLocator = MultipleLocator(xtickgrid)
    ymajorLocator = MultipleLocator(ytickgrid)
    for ax in axis_list:
        ax.grid()
        ax.set_xlim([xhi, xlo])
        ax.set_ylim([ylo, yhi])
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
    for k in range(len(axis_list)):
        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
                string.ascii_uppercase[k], fontweight='bold')
    ax8.set_xlabel('Wavenumber (cm$^{-2}$)')
    for ax in [ax1, ax4, ax7]:
        plt.setp(ax.get_yticklabels(), visible=True)
    for ax in [ax7, ax8, ax9]:
        plt.setp(ax.get_xticklabels(), visible=True, rotation=45)    

#    blue_line = mlines.Line2D([], [], label='Observed', 
#                              **sp.style_spectrum)
#    green_line = mlines.Line2D([], [], label='Fit peaks', **sp.style_fitpeak)
#    dashed_line = mlines.Line2D([], [], label='Sum of fit peaks', 
#                                **sp.style_summed)
#    leg = axis_list[7].legend(handles=[blue_line, green_line, dashed_line], ncol=3,
#            loc = 'upper center', bbox_to_anchor = (0.5, -0.5), fancybox=True)
        
    return axis_list

def plotsetup_3stacked(yhi=2, ylo=0, xlo = 3200, xhi = 3800,
                  xtickgrid=200, ytickgrid=0.5):
    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
    fig = plt.figure()
    fig.set_size_inches(3, 6.5)
    gs = gridspec.GridSpec(3, 1)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0])
    axis_list = [ax1, ax2, ax3]
    gs.update(bottom=0.14, left=0.35, right=0.95, top=0.95)    
    xmajorLocator = MultipleLocator(xtickgrid)
    ymajorLocator = MultipleLocator(ytickgrid)
    for ax in axis_list:
        ax.grid()
        ax.set_xlim([xhi, xlo])
        ax.set_ylim([ylo, yhi])
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.setp(ax.get_xticklabels(), rotation=45)    
    for k in range(len(axis_list)):
        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
                string.ascii_uppercase[k], fontweight='bold')
    ax3.set_xlabel('Wavenumber (cm$^{-2}$)')
    for ax in [ax1, ax2]:
        plt.setp(ax.get_xticklabels(), visible=False)
    return axis_list
