"""
FTIR class definitions and functions for processing and plotting FTIR spectra

The main unit is a class called Spectrum, which includes attributes 
such as sample thickness and lists of wavenumbers and absorbances. 
A few defaults are set up that are geared toward H in nominally anhydrous
minerals (NAMs) such a plotting wavenumber range from 3000 to 4000 /cm and
creating baselines between 3200 and 3700 /cm. 

More detailed instructions and examples are to be written. For an example 
of how to instantiate these classes, check out my_spectra.py in my 
CpxCode GitHub repository.

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
import matplotlib.lines as mlines
from uncertainties import ufloat
import os.path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import string as string
import lmfit
import scipy
import uncertainties
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import xlsxwriter

# Optional: Set default folder where to find and/or save files for use with 
# functions like save_spectrum and save_baseline. Eventually I will get 
# rid of this, but for now having it hard-coded here is helpful for me.
default_folder = 'C:\\Users\\Ferriss\\Documents\\FTIR\\'


## plotting style defaults
style_baseline = {'color' : 'k', 'linewidth' : 1, 'linestyle' :'-'}
style_spectrum = {'color' : 'b', 'linewidth' : 3}
style_fitpeak = {'color' : 'g', 'linewidth' : 1}
style_summed = {'color' : 'orangered', 'linewidth' : 2, 'linestyle' : '--'}
style_profile = {'markeredgecolor' : 'black', 'linestyle' : '', 'marker' : 'o', 
                 'markersize' : 10, 'markerfacecolor' : 'gold', 'alpha' : 0.5}
style_initial = {'color' : 'green', 'linewidth' : 1}
style_1 = {'linestyle' : '-', 'color' : 'k', 'marker' : None}
style_unoriented = {'fillstyle' : 'none'}
# different profile directions
style_Dx = {'fillstyle' : 'left', 'color' : 'red', 'markerfacecolor' : 'red'}
style_Dy = {'fillstyle' : 'bottom', 'color' : 'green', 
            'markerfacecolor' : 'green' }
style_Dz = {'fillstyle' : 'right', 'color' : 'blue', 
            'markerfacecolor' : 'blue'}
style_Dx_line = {'linestyle' : '--', 'color' : 'red'}
style_Dy_line = {'linestyle' : '-.', 'color' : 'green'}
style_Dz_line = {'linestyle' : ':', 'color' : 'blue'}
style_unoriented_line = {'linestyle' : '-', 'color' : 'black'}
# different ray paths
style_Rx = {'marker' : 'd'}
style_Ry = {'marker' : 'o'}
style_Rz = {'marker' : 's'}

def make_line_style(direction, style_marker):
    """Take direction and marker style and return line style dictionary
    that reflects the direction (x, y, z, or u for unoriented) with the 
    color of the base style"""
    if direction == 'x':
        d = style_Dx_line
    if direction == 'y':
        d = style_Dy_line
    if direction == 'z':
        d = style_Dz_line
    if direction == 'u':
        d = style_unoriented_line        
    d.update({'linewidth' : 2})
    return d
    

# Define classes, attributes, and functions related to samples
class Sample():
    mineral_name = None
    source = None
    initial_water = None
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

def area2water(area_cm2, mineral='cpx', calibration='Bell'):
    """Takes area in cm-2, multiplies by absorption coefficient, 
    return water concentration in ppm H2O"""
    if calibration == 'Bell':
        cpx_calib_Bell95 = 1.0 / ufloat(7.09, 0.32)
    else:
        print 'only Bell calibration so far'
        return
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
    area = None
    # peak fit information
    peakpos = None
    numPeaks = None
    peak_heights = None
    peak_widths = None
    peak_areas = None
    peakshape = None        
    # Experimental information if applicable
    temperature_celcius = None
    time_seconds = None
    diffusivity_log10m2s = None
    
    def set_thick(self):
        """Set spectrum thick_microns from raypath and sample.thick_microns.
        Similar utility exists for Profiles"""
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
        to the new spectrum (self)"""
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
                      show_fit_values=False, show_plot=False):
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
        
        if show_plot is True:
            self.plot_showbaseline()
            
        return base_abs

    def subtract_baseline(self, polyorder=1, bline=None, 
                          show_plot=False, shifter=None):
        """Make baseline and return baseline-subtracted absorbance.
        Preferred baseline used is (1) input bline, (2) self.base_abs, 
        (3) one made using polyorder, default is linear"""
        if self.wn is None:
            self.start_at_zero()
        
        if bline is not None:
            base_abs = bline
        elif self.base_abs is not None:
            print 'Subtracting baseline stored in attribute base_abs'
            base_abs = self.base_abs
        else:
            print 'Making baseline or polynomial order', polyorder
            base_abs = self.make_baseline(polyorder, shiftline=shifter)
        
        index_lo = (np.abs(self.wn-self.base_low_wn)).argmin()
        index_hi = (np.abs(self.wn-self.base_high_wn)).argmin()
            
        humps = self.abs_cm[index_lo:index_hi]
        abs_nobase_cm = humps - base_abs

        # plotting
        if show_plot is True:
            fig, ax = self.plot_spectrum_outline()
            plt.plot(self.wn, self.abs_cm, **style_spectrum) 
            ax.plot(self.base_wn, base_abs, **style_baseline)
            yhigh = max(self.abs_cm) + 0.1*max(self.abs_cm)
            if min(base_abs) > 0:
                ylow = 0
            else:
                ylow = min(base_abs) + 0.1*min(base_abs)
            plt.ylim(ylow, yhigh)
            plt.show(fig)
            
        return abs_nobase_cm

    def area_under_curve(self, polyorder=1, area_plot=True, shiftline=None,
                         printout=True, numformat='{:.1f}', 
                         require_saved_baseline=False):
        """Returns area under the curve in cm^2"""
        if require_saved_baseline is True:
            self.get_baseline()
            return

        if (self.base_high_wn is None) or (self.base_low_wn is None):
            print 'Need to specify spectrum baseline wavenumber range'
            return False        

        # Get or make absorbance and wavenumber range for baseline
        if (self.base_abs is None) or (self.base_wn is None):
            print 'Making baseline...'
            abs_baseline = self.make_baseline(polyorder, shiftline)
        else:
            print ('Using previously fit baseline. ' + 
                   'Use make_baseline to change it.')
            abs_baseline = self.base_abs

        abs_nobase_cm = self.subtract_baseline(bline=abs_baseline,
                                               show_plot=area_plot)

        dx = self.base_high_wn - self.base_low_wn
        dy = np.mean(abs_nobase_cm)
        area = dx * dy
        if printout is True:
            print self.fname
            print 'area:', numformat.format(area), '/cm^2'
            
        self.area = area
        return area

    def water(self, polyorder=1, mineral_name='cpx', numformat='{:.1f}',
              show_plot=True, show_water=True, printout_area=False,
              shiftline=None):
        """Produce water estimate without error for a single FTIR spectrum. 
        Use water_from_spectra() for errors, lists, and 
        non-Bell calibrations."""
        # determine area
        area = self.area_under_curve(polyorder, show_plot, shiftline, 
                                     printout_area)
        w = area2water(area, mineral=mineral_name)
                        
        # output for each spectrum
        if show_water is True:
            print self.fname
            print 'area:', numformat.format(area), '/cm^2'
            print 'water:', numformat.format(w), 'ppm H2O'
        return area, w

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

        print 'Baseline information taken from file -baseline.CSV'
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
                
    def get_peakfit(self, folder=default_folder, delim=',', 
                    file_ending='-peakfit.CSV'):
        """Get individual peaks for spectrum from fname-peakfit.CSV generated
        in MATLAB using FTIR_peakfit_afterpython.m and savefit.m
        Return curves (assumes Gaussians) and summation"""
        filename = folder + self.fname + file_ending
        if os.path.isfile(filename) is True:
            previous_fit = np.loadtxt(filename, delimiter=delim)
            self.peakpos = previous_fit[:, 0]
            self.numPeaks = len(self.peakpos)
            self.peak_heights = previous_fit[:, 1]
            self.peak_widths = previous_fit[:, 2]
            self.peak_areas = previous_fit[:, 3]
#            self.peakshape = previous_fit[:, 4]
        else:
            print ' '            
            print self.fname
            print 'Remember to savefit.m after FTIR_peakfit_afterpython.m'
            return False

        if self.base_wn is None:
            self.get_baseline()
            
        peakfitcurves = np.ones([len(self.peakpos), len(self.base_wn)])
        summed_spectrum = np.zeros_like(self.base_wn)
        for k in range(len(self.peakpos)):
            peakfitcurves[k] = make_gaussian(x=self.base_wn, 
                                            pos=self.peakpos[k],
                                            h=self.peak_heights[k], 
                                            w=self.peak_widths[k])
            summed_spectrum += peakfitcurves[k]            
        return peakfitcurves, summed_spectrum
    
    def plot_spectrum_outline(self, size_inches=(3., 2.5), shrinker=0.15):
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
        
        box = ax.get_position()
        ax.set_position([box.x0 + box.width*shrinker, 
                         box.y0 + box.height*shrinker, 
                         box.width*(1.0-shrinker), 
                         box.height*(1.0-shrinker)])

        return f, ax

    def plot_spectrum(self, style=style_spectrum, figaxis=None, wn=None):
        """Plot the raw spectrum divided by thickness"""
        if figaxis is None:
            fig, ax = self.plot_spectrum_outline()
        else:
            fig = None
            ax = figaxis
            
        if self.wn is None:
            self.start_at_zero()

        ax.plot(self.wn, self.abs_cm, **style)
        if wn is not None:
            ax.plot([wn, wn], [ax.get_ylim()[0], ax.get_ylim()[1]], color='r')
        return fig, ax
    
    def plot_peakfit(self, style=style_spectrum, stylesum=style_summed, 
                     stylepeaks=style_fitpeak, top=None, legloc=1):
        """Plot peaks fit in MATLAB using peakfit.m"""
        # Take baseline-subtracted spectrum from saved file every time 
        # to avoid any possible funny business from playing with baselines
        gaussian, summed_spectrum = self.get_peakfit()
        wn, observed = self.get_baseline()

        fig, ax = self.plot_subtractbaseline(style=style) # just plot setup
        ax.plot(wn, observed, label='observed', **style)
        ax.plot(wn, gaussian[0], label='fit bands', **stylepeaks)
        ax.plot(self.base_wn, summed_spectrum, label='peak sum', **stylesum)
        ax.legend(loc=legloc)
        
        if top is None:
            topnat = ax.get_ylim()[1]
            top = topnat + topnat*0.75
            
        ax.set_ylim(0., top)        
        
        for k in range(len(self.peakpos)):
            ax.plot(self.base_wn, gaussian[k], **stylepeaks)

        return fig, ax
        
    def plot_showbaseline(self, polyorder=1, shiftline=None,
                          wn_baseline=None, abs_baseline=None,
                          style=style_spectrum, style_base=style_baseline,
                          figaxis=None):
        """Plot FTIR spectrum and show baseline. 
        You can pass in your own baseline"""
        if figaxis is None:
            fig, ax = self.plot_spectrum_outline()
        else:
            fig = None
            ax = figaxis
        
        # Get or make absorbance for baseline
        if abs_baseline is None:          
            if self.base_abs is None:
                print 'Making baseline...'
                abs_baseline = self.make_baseline(polyorder, shiftline)
            else:
                abs_baseline = self.base_abs

        # get baseline wavenumber range
        if wn_baseline is None:
            if self.base_wn is not None:
                wn_baseline = self.base_wn
            else:
                print ('Need to pass in baseline wavenumber range too ' +
                       'either here as wn_baseline or as spectrum.base_wn')
                return                
            
        ax.plot(self.wn, self.abs_cm, **style)
        ax.plot(wn_baseline, abs_baseline, **style_base)
        return fig, ax

    def plot_subtractbaseline(self, polyorder=1, style=style_spectrum, 
                              figaxis=None):
        """Make and plot baseline-subtracted spectrum"""
        abs_nobase_cm = self.subtract_baseline(polyorder)
        
        if figaxis is None:
            fig, ax = self.plot_spectrum_outline()
        else:
            fig = None
            ax = figaxis
            
        pad = max(abs_nobase_cm)
        yhigh = pad + 0.1*pad
        ax.set_ylim(0, yhigh)
        ax.set_xlim(self.base_high_wn, self.base_low_wn)
        ax.plot(self.base_wn, abs_nobase_cm, **style)
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
        print 'Bell calibration'
        pass
    elif calibration == 'Paterson':
        print 'Paterson calibration not available quite yet...'
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
            print 'Showing spectrum for', spec.fname
            fig, ax = spec.plot_spectrum_outline()
            plt.plot(spec.wn, spec.abs_cm, **style_spectrum) 
        
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
#        
        k = 0
        for bam in [window_small, -window_small, -window_large]:
            # Setting up 3 different baseline shapes
            main_yshift = main_yshift - bam
            base_abs = spec.make_baseline(2, shiftline=main_yshift)
            abs_nobase_cm = spec.subtract_baseline(bline=base_abs)
            if show_plots is True:
                ax.plot(spec.base_wn, base_abs, **style_baseline)

            area = spec.area_under_curve(area_plot=False, printout=False)
            area_list = np.append(area_list, area)
            baseline_list[k,:] = base_abs
            pek_list[k,:] = abs_nobase_cm
            k+=1

        # list areas and water with uncertainties
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
    profile_name = None
    initial_profile = None
    spectrum_class_name = None
    # averaged spectra made by self.make_average()
    avespec = None
    iavespec = None
    # Made by make_spectra_list
    positions_microns = np.array([])
    spectra_list = []
    len_microns = None
    thick_microns = None
    # Made by make_water_list etc. below
    waters_list = []
    waters_errors = []
    areas_list = []
    bestfitline_areas = None
    # Whole-block 3D-WB information if applicable
    wb_areas = None # bulk hydrogen
    wb_water = None
    wb_initial_profile = None
    # for plotting
    style_base = style_profile
    style_x_marker = None
    style_y_marker = None
    style_z_marker = None
    style_x_line = None
    style_y_line = None
    style_z_line = None
    # Experimental information
    temperature_celcius = None
    time_seconds = None
    diffusivity_log10m2s = None
    diff_error = None
    # peak fitting. i = initial
    peakpos = None
    peak_heights = None
    peak_widths = None
    peak_areas = None
    peak_iheights = None
    peak_iwidths = None
    peak_iareas = None
    peak_wb_areas = None
    peak_wb_heights = None
    peak_wb_widths = None
    peak_diffusivities = None
    peak_diff_error = None
    peak_initial_wb_ratio = 1.0

    def set_len(self):
        """Set profile.len_microns from profile.direction and 
        profile.sample.thick_microns""" 

        if self.sample is None:
            print '\n', self.profile_name
            print 'Need to specify profile sample\n'
            return False
        else:
            s = self.sample

        if s.sample_thick_microns is None:
            s.sample_thick_microns = get_3thick(s)                
        
        if self.direction == 'a':
           self.len_microns = s.sample_thick_microns[0]
        elif self.direction == 'b':
            self.len_microns = s.sample_thick_microns[1]
        elif self.direction == 'c':
            self.len_microns = s.sample_thick_microns[2]
        else:
            print 'Set direction of profile to a, b, or c'

        return self.len_microns

    def set_thick(self):
        """Set profile.thick_microns from profile.raypath and
        profile.sample.thick_microns"""
        if self.sample is None:
            print 'Need to specify profile sample'
            return False
        else:
            s = self.sample
        
        if s.sample_thick_microns is None:
            s.sample_thick_microns = get_3thick(s)

        if self.raypath == 'a':
           self.thick_microns = s.sample_thick_microns[0]
        elif self.raypath == 'b':
            self.thick_microns = s.sample_thick_microns[1]
        elif self.raypath == 'c':
            self.thick_microns = s.sample_thick_microns[2]
        else:
            print 'Need raypath'
            return False
            
        return self.thick_microns

    def make_spectra_list(self, class_from_module=None, 
                          my_module='my_spectra'):
        """Set profile length and generate spectra_list 
        with key attributes"""
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
        else:
            classToUse = Spectrum    

        if len(self.spectra_list) == 0:
            if len(self.fname_list) == 0:
                print 'Need fnames'
                return False                
            fspectra_list = []
            for x in self.fname_list:
                newspec = classToUse()
                newspec.fname = x
                newspec.sample = self.sample
                newspec.raypath = self.raypath
                newspec.thick_microns = self.thick_microns
                fspectra_list.append(newspec)
            self.spectra_list = fspectra_list
            
        else:
            for spec in self.spectra_list:
                spec.sample = self.sample
                spec.raypath = self.raypath
                spec.thick_microns = self.thick_microns
                
        return

    def average_spectra(self):
        """Creates and returns averaged spectrum and stores it in
        attribute avespec"""
        spec_list = self.spectra_list
        avespec = Spectrum()
        avespec.make_average_spectra(spec_list)
        
        if self.profile_name is not None:
            avespec.fname = (self.profile_name + '\naverage profile')
        else:
            avespec.fname = 'average profile'
            
        avespec.base_high_wn = spec_list[0].base_high_wn
        avespec.base_low_wn = spec_list[0].base_low_wn
        avespec.start_at_zero()
        return avespec

    def plot_spectra(self, show_baseline=True, show_initial_ave=True, 
                     show_final_ave=True, plot_all=False, 
                     initial_and_final_together=False, style=style_spectrum, 
                     stylei=style_initial, wn=None):
        """Plot averaged spectrum across profile. Returns figure and axis."""
        if self.spectra_list is None:
            self.make_spectra_list()

        if show_initial_ave is True or initial_and_final_together is True:
            if self.initial_profile is None:
                print 'No initial_profile attribute specified'
                show_initial_ave = False
                initial_and_final_together = False
        
        f = None
        ax = None
        
        if plot_all is True:
            for spec in self.spectra_list:
                if show_baseline is True:
                    spec.plot_showbaseline(style=style)
                else:
                    spec.plot_spectrum(style=style, wn=wn)
            if show_initial_ave is True:
                for spec in self.initial_profile.spectra_list:
                    if show_baseline is True:
                        spec.plot_showbaseline(style=stylei)
                    else:
                        spec.plot_spectrum(style=stylei, wn=wn)
                        
        if show_final_ave is True or initial_and_final_together is True:
            avespec = self.average_spectra()
        
        if show_initial_ave is True or initial_and_final_together is True:
            initspec = self.initial_profile.average_spectra()
            
        # Plot initial average
        if show_initial_ave is True:
            if show_baseline is True:
                initspec.plot_showbaseline(style=stylei)
            else:
                initspec.plot_spectrum(style=stylei, wn=wn)

        # Plot final average
        if show_final_ave is True:
            if show_baseline is True:
                f, ax = avespec.plot_showbaseline(style=style)
            else:            
                f, ax = avespec.plot_spectrum(style=style, wn=wn)
        
        # Plot average spectra together
        if initial_and_final_together is True:
            f, ax = avespec.plot_spectrum_outline()
            ax.plot(avespec.wn, avespec.abs_cm, label='Final', **style)
            ax.plot(initspec.wn, initspec.abs_cm, label='Initial', **stylei)            
            ax.legend()
            tit = self.profile_name + '\nAverage profiles'
            ax.set_title(tit)
            if wn is not None:
                ax.plot([wn, wn], [ax.get_ylim()[0], ax.get_ylim()[1]], 
                        color='r')
        return f, ax

    def change_baseline(self, highwn=3800, lowwn=3000, shift=None):
        """Change baseline parameters for all spectra, final and initial"""
        for spectrum in self.spectra_list + self.initial_profile.spectra_list:
            spectrum.base_high_wn = highwn
            spectrum.base_low_wn = lowwn
            if shift is not None:
                spectrum.base_mid_yshift = shift

    def make_baselines(self, lineorder=1, shiftline=None, 
                       show_fit_values=False, show_plot=False):
        """Make baselines for all final and initial spectra""" 
        for spectrum in self.spectra_list + self.initial_profile.spectra_list:
            spectrum.make_baseline(line_order=lineorder, shiftline=shiftline,
                                    show_fit_values=show_fit_values,
                                    show_plot=False)

    def get_baselines(self, initial_too=False):
        """Get previously saved baselines for all spectra in profile"""
        for spectrum in self.spectra_list:
            spectrum.get_baseline()
            
        if initial_too is True:
            for spectrum in self.initial_profile.spectra_list:
                spectrum.get_baseline()

    def matlab_list(self):
        """Print out spectra fnames for FTIR_peakfit_loop.m"""
        string = "{"
        for spec in self.spectra_list:
            stringname = spec.fname
            string = string + "'" + stringname + "' "
        string = string + "};"
        print '\nfilenames ready for FTIR_peakfit_loop.m:'
        print string
                        
    def save_baselines(self, printnames=True):
        """Save all baselines in profile"""
        for spectrum in self.spectra_list:
            spectrum.save_baseline()
        
        if printnames is True:
            self.matlab_list()            
            
    def make_area_list(self, polyorder=1, show_plot=False, set_class=None,
                       shiftline=None, printout_area=False, peak=None):
        """Make list of areas (no error) under the curve for an FTIR profile.
        Default is bulk area. Set peak=wavenumber for peak-specific profile"""
        # You need the list of spectra for the profile
        if len(self.spectra_list) < 1:
            check = self.make_spectra_list(class_from_module=set_class)
            if check is False:
                return False
              
        areas = []
        if peak is None:
            if len(self.areas_list) < 1:
                print 'generating bulk areas under curves...'
                for x in self.spectra_list:
                    if x.abs_nobase_cm is None:
                        x.get_baseline()
                    
                    a = x.area_under_curve(polyorder, show_plot, shiftline, 
                                           printout_area, 
                                           require_saved_baseline=False)
                    areas.append(a)            
                self.areas_list = areas
            else:
                areas = self.areas_list
        else:
            peaklist = list(self.spectra_list[0].peakpos)
            print 'peak at', peak

            if peak in peaklist:
                idx = peaklist.index(peak)
                for x in self.spectra_list:
                    a = x.peak_areas[idx]
                    areas.append(a)
            else:
                print 'No peak at wavenumber', peak
                print 'peak positions:', peaklist
                return False
        return areas
        
    def get_peakfit(self):
        """Get fit peaks from MATLAB for all spectra in profile, including 
        in the initial profile. The resulting numpy arrays are in dimensions
        (number of peaks, number of spectra in profile) and stored in
        the profiles's attributes peak_heights, peak_widths, and peak_areas"""
        for spectrum in self.spectra_list:
            spectrum.get_peakfit()
        peakpos = self.spectra_list[0].peakpos

        hbig = []
        wbig = []
        abig = []
        
        if peakpos is None:
            return False
            
        for p in range(len(peakpos)):
            h = []
            w = []
            a = []        
            for spectrum in self.spectra_list:
                h.append(spectrum.peak_heights[p])
                w.append(spectrum.peak_widths[p])
                a.append(spectrum.peak_areas[p])
            hbig.append(h)
            wbig.append(w)
            abig.append(a)
                
        self.peak_heights = np.array(hbig)
        self.peak_widths = np.array(wbig)
        self.peak_areas = np.array(abig)
        self.peakpos = peakpos
        if self.peak_diffusivities is None:
            self.peak_diffusivities = np.zeros_like(peakpos)
            self.peak_diff_error = np.zeros_like(peakpos)

    def plot_peakfits(self, initial_too=False, legloc=1):
        """Show fit peaks for all spectra in profile"""
        for spectrum in self.spectra_list:
            spectrum.plot_peakfit(legloc=legloc)
        
        if initial_too is True and self.initial_profile is not None:
            for spectrum in self.initial_profile.spectra_list:
                spectrum.plot_peakfit(legloc=legloc)
                
    def make_wholeblock(self, peakfit=False, show_plot=False, bulk=True):
        """Take initial and final profiles and make self.wb_areas.
        If peakfit=True, then get_peakfit and make peak-specific ratios
        peak_wb_areas, peak_wb_heights, peak_wb_widths"""

        if self.initial_profile is None:
            print '\n profile name:', self.profile_name
            print 'Need to specify an initial profile'
            return False

        # Bulk H whole-block
        if bulk is True:
            # get both initial and final area lists
            init = self.initial_profile
            fin = self
            for prof in [init, fin]:
                if len(prof.areas_list) < 1:
                    check = prof.make_area_list(1, show_plot=False)
                    if check is False:
                        return False

                if ((len(init.areas_list) == 2) and 
                    (init.positions_microns[0] == init.positions_microns[1])):
                    p = (0, np.mean(init.areas_list))
                elif len(init.areas_list) > 1:
                    p = np.polyfit(init.positions_microns, init.areas_list, 1)
                else:
                    p = (0, init.areas_list[0])
            
            init_line = np.polyval(p, self.positions_microns)
            area_ratio = self.areas_list / init_line
            self.wb_areas = area_ratio
        if peakfit is False:
            return area_ratio
    
        # Peak-specific whole-block values
        if self.peak_areas is None:
            self.get_peakfit()
            if self.peak_areas is None:
                print '\nCould not get peak areas'
                return
            
        if self.initial_profile.peak_areas is None:
            self.initial_profile.get_peakfit()
            if self.initial_profile.peak_areas is None:
                print '\nCould not get initial peak areas'
                return
            
        ipos = self.initial_profile.positions_microns
        iareas = self.initial_profile.peak_areas
        iheights = self.initial_profile.peak_heights
        iwidths = self.initial_profile.peak_widths
        pos = self.positions_microns       
        areas = self.peak_areas
        heights = self.peak_heights
        widths = self.peak_widths        
        npeaks = len(iareas)
        
        wb_areas = []
        wb_heights = []
        wb_widths = []
        for k in range(npeaks):
            pa = np.polyfit(ipos, iareas[k], 1)
            ph = np.polyfit(ipos, iheights[k], 1)
            pw = np.polyfit(ipos, iwidths[k], 1)
            anormalizeto = np.polyval(pa, pos)
            hnormalizeto = np.polyval(ph, pos)
            wnormalizeto = np.polyval(pw, pos)
            
            if areas is not None:
                wb_areas.append(areas[k] / anormalizeto)
            else:
                print '\nProblem making whole-block area ratios\n'
                return
                
            wb_heights.append(heights[k] / hnormalizeto)
            wb_widths.append(widths[k] / wnormalizeto)
        self.peak_wb_areas = wb_areas
        self.peak_wb_heights = wb_heights
        self.peak_wb_widths = wb_widths
        return wb_areas

    def get_peak_wb_areas(self, peak_idx=0, peakwn=None, 
                          heights_instead=False):
        """Returns peak-specific whole-block areas for the profile
        AND peak wavenumber in cm-1 because usually the peak index is easier
        to pass.
        Defaults to the first peak in the peak position list"""
        # Add check that argument is actually a profile
        
        if self.peakpos is None: 
            print 'Getting peaks fit in matlab for', self.profile_name
            self.get_peakfit()

        # Getting peak from wavenumber if index not given
        # Or peak wavenumber from index
        if peak_idx is None:
            if peakwn not in self.peakpos:
                print 'There is no peak at', peakwn
                print 'Peaks are at', self.peakpos
                return
            peak_idx = np.where(self.peakpos==peakwn)[0][0]
            print 'peak at', peakwn, 'is index', peak_idx
        else:
            peakwn = self.peakpos[peak_idx]

        # Make peak-specific whole-block areas if needed
        if self.peak_wb_areas is None:
            self.make_wholeblock(peakfit=True, show_plot=False, 
                                 bulk=False)
               
        if heights_instead is False:
            returnvals = self.peak_wb_areas[peak_idx]
        else:
            returnvals = self.peak_wb_heights[peak_idx]
    
        return returnvals, peakwn

    def make_style_subtypes(self):
        """Make direction-specific marker and line style dictionaries for 
        profile based on profile's self.style_base"""
        if self.style_base is None:
            print 'Set base style (style_base) for profile first'
            return False
        self.style_x_marker =  dict(self.style_base.items() + style_Dx.items())
        self.style_y_marker =  dict(self.style_base.items() + style_Dy.items())
        self.style_z_marker =  dict(self.style_base.items() + style_Dz.items())
        self.style_x_line = make_line_style('x', self.style_base)
        self.style_y_line = make_line_style('y', self.style_base)
        self.style_z_line = make_line_style('z', self.style_base)
        
        if self.raypath == 'a':
            self.style_x_marker.update({'marker' : style_Rx['marker']})
            self.style_y_marker.update({'marker' : style_Rx['marker']})
            self.style_z_marker.update({'marker' : style_Rx['marker']})
        elif self.raypath == 'b':
            self.style_x_marker.update({'marker' : style_Ry['marker']})
            self.style_y_marker.update({'marker' : style_Ry['marker']})
            self.style_z_marker.update({'marker' : style_Ry['marker']})
        elif self.raypath == 'c':
            self.style_x_marker.update({'marker' : style_Rz['marker']})        
            self.style_y_marker.update({'marker' : style_Rz['marker']})
            self.style_z_marker.update({'marker' : style_Rz['marker']})
        return

    def choose_line_style(self):
        """Returns line style with direction information"""
        if self.style_base is None:
            print 'Using default styles. Set profile style_base to change'
            self.style_base = style_profile            
        if self.style_x_line is None:
            self.make_style_subtypes()
        if self.direction == 'a':
            style_bestfitline = self.style_x_line
        elif self.direction == 'b':
            style_bestfitline = self.style_y_line
        elif self.direction == 'c':
            style_bestfitline = self.style_z_line
        else:
            style_bestfitline = {'linestyle' : '-'}
        return style_bestfitline
    
    def choose_marker_style(self):
        """Returns marker style with direction and ray path information"""
        if self.style_base is None:
            print 'Using default styles. Set profile style_base to change'
            self.style_base = style_profile
        if self.style_x_marker is None:
            self.make_style_subtypes()
        if self.direction == 'a':
            style = self.style_x_marker
        elif self.direction == 'b':
            style = self.style_y_marker
        elif self.direction == 'c':
            style = self.style_z_marker
        return style
        
    def plot_area_profile_outline(self, centered=True, peakwn=None,
                                  figsize=(6.5, 2.5), top=None):
        """Set up area profile outline and style defaults. 
        Default is for 0 to be the middle of the profile (centered=True)."""
        if self.style_base is None:
            self.style_base = style_profile
        self.make_style_subtypes()
        
        if self.len_microns is None:
            leng = self.set_len()
        else:
            leng = self.len_microns

        f, ax = plt.subplots(1, 1)
        f.set_size_inches(figsize)
        ax.set_xlabel('Position ($\mu$m)')
        ax.set_ylabel('Area (cm$^{-2}$)')

        if top is None:
            if len(self.areas_list) > 1:
                top = max(self.areas_list)+0.2*max(self.areas_list)
            else:
                top = 1.
            
        ax.set_ylim(0, top)

        if centered is True:
            ax.set_xlim(-leng/2.0, leng/2.0)
        else:
            ax.set_xlim(0, leng)
                    
        ax.grid()
        return f, ax

    def plot_area_profile(self, polyorder=1, centered=True, 
                          heights_instead = False, figaxis=None, 
                          bestfitline=False, show_FTIR=False, 
                          show_values=False, set_class=None,
                          peakwn=None, peak_idx=None,
                          style=None, show_initial_areas=False,
                          error_percent=2.):
        """Plot area profile. Centered=True puts 0 in the middle of the x-axis.
        figaxis sets whether to create a new figure or plot on an existing axis.
        bestfitline=True draws a best-fit line through the data.
        Set peak=wavenumber for peak-specific profile"""
        # Check for or create positions and areas
        if len(self.positions_microns) < 1:
            print 'Need positions_microns for profile'
            return

        # Get list of areas            
        if peakwn is None and peak_idx is None:
            # bulk hydrogen
            self.get_baselines()
            areas = self.make_area_list(polyorder, show_FTIR, 
                                    set_class, peak=peakwn)

        else:
            # peak-specific
            if self.peak_areas is None:
                self.get_peakfit()

            if peak_idx is None:
                peak_idx = np.where(self.peakpos==peakwn)[0][0]
                print 'peak at', peakwn, 'is index', peak_idx
            else:
                peakwn = self.peakpos[peak_idx]

        if heights_instead is True and peak_idx is not None:
            areas = self.peak_heights[peak_idx]
        else:            
            areas = self.peak_areas[peak_idx]
                   
        if areas is False:
            return            
        if np.shape(areas) != np.shape(self.positions_microns):
            print 'Area and positions lists are not the same size!'
            print 'area:', np.shape(self.areas_list)
            print 'positions:', np.shape(self.positions_microns)
            return

        # Use new or old figure axes
        if figaxis is None:
            f, ax = self.plot_area_profile_outline(centered, peakwn)
        else:
            ax = figaxis

        # Set length
        if self.len_microns is None:
            leng = self.set_len()
        else:
            leng = self.len_microns

        # Set up plotting styles
        style_bestfitline = self.choose_line_style()
        if style is None:
            style = self.choose_marker_style()

        # Plot best fit line beneath data points
        if bestfitline is True:
            if ((len(areas) == 2) and 
                (self.positions_microns[0] == self.positions_microns[1])):
                p = (0, np.mean(areas))
            elif len(areas) > 1:
                p = np.polyfit(self.positions_microns, areas, 1)
            else:
                p = (0, areas[0])

            self.bestfitline_areas = p
            x = np.linspace(0, leng, 100)
            y = np.polyval(p, x)
            if centered is True:
                ax.plot(x-leng/2.0, y, **style_bestfitline)
            else:
                ax.plot(x, y, **style_bestfitline)

        # Plot data
        if centered is True:
            x = self.positions_microns - leng/2.0
        else:
            x = self.positions_microns            
        
        yerror = np.array(areas)*error_percent/100.
        ax.errorbar(x, areas, yerr=yerror, ** style)
            
        # Plot initial profile areas
        if show_initial_areas is True:
            self.initial_profile.plot_area_profile(style={'color' : 'r', 
                                                          'marker' : 'o'},
                                                   figaxis=ax)            
        ax.set_ylim(0, max(areas)+0.2*(max(areas)))

        if heights_instead is True and peak_idx is not None:
            ax.set_ylabel('Peak height (/cm)')

        if figaxis is None:
            return f, ax
        else:
            return


    def plot_wholeblock_profile(self, polyorder=1, centered=True, figaxis=None, 
                                bestfitline=False, show_FTIR=False, 
                                show_values=False, set_class=None,
                                peakwn=None, peak_idx=None,
                                style=None, error_percent=3., xerror=50.,
                                show_initial_areas=False, top=1.0, 
                                heights_instead=True):
        """Plot whole-block profile"""
        # Bulk hydrogen
        if peak_idx is None:           
            if self.wb_areas is None:
                self.make_wholeblock(peakfit=False, bulk=True)
            a = self.wb_areas
        # peak-specific selected by peak index
        else:
            if self.peak_wb_areas is None:
                self.make_wholeblock(peakfit=True, bulk=False)

            if heights_instead is False:
                a = self.peak_wb_areas[peak_idx]
            else:
                a = self.peak_wb_heights[peak_idx]
                
            peakwn = self.peakpos[peak_idx]            
        
        # Use new or old figure axes
        if figaxis is None:
            f, ax = self.plot_area_profile_outline(centered=centered,
                                                   peakwn=peakwn)
        else:
            ax = figaxis
            if heights_instead is False:
                ax.set_ylabel('Area/Area$_0$')
            else:
                ax.set_ylabel('Height/Height$_0$')
        
        ax.set_ylim(0, top)

        if self.len_microns is None:
            leng = self.set_len()
        else:
            leng = self.len_microns

        if style is None:
            style = self.choose_marker_style()

        # Plot data            
        if centered is True:
            x = self.positions_microns - (leng/2.0)
        else:
            x = self.positions_microns
        
        yerror = np.array(a)*error_percent/100.
        ax.errorbar(x, a, xerr=xerror, yerr=yerror, **style)
        
        if figaxis is None:
            return f, ax
        else:
            return
        
#        # Multiply by initial water to get scale up from area ratio        
#        if (self.sample is None) or (self.sample.initial_water is None):
#            print ('Need self.sample.initia_water for for A/A0 to water' +
#                    'Returning only area ratios')
#            return area_ratio
#        water_profile = area_ratio * self.sample.initial_water
#        return area_ratio, water_profile

#    def plot_wb_water(self, polyorder=1, centered=True, style=style_profile):
#        """Plot area ratios and scale up with initial water"""
#        if self.sample.initial_water is None:
#            print 'Set sample initial water content please.'
#            return
#        a, w = self.make_3DWB_water_list(polyorder=1)
#        
#        fig, ax, leng = self.plot_area_profile_outline(centered)
#        top = max(a) + 0.2*max(a)
#        ax.set_ylim(0, top)
#        ax.set_ylabel('Final area / Initial area')
#        
#        if centered is True:
#            ax.plot([-leng/2.0, leng/2.0], [1, 1], **style_1)
#            ax.plot(self.positions_microns - leng/2.0, a, **style)           
#        else:
#            ax.plot([0, leng], [1, 1], **style_1)
#            ax.plot(self.positions_microns, a, **style)
#        return fig, ax, leng

    def fitDiffusivity(self, initial_unit_value=1., vary_initial=False,
                       guess=-13., peak_idx=None, top=1.2, peakwn=None,
                       show_plot=True, polyorder=1, heights_instead=False):
        """Fits 1D diffusion curve to profile data.
        I still need to add almost all of the diffusion_1D keywords."""
        if self.time_seconds is None:        
            print 'Need to set profile time_seconds attribute'
            return
        if self.positions_microns is None:
            print 'Need to set profile positions'
            return
        
        # Bulk hydrogen
        if peak_idx is None:
            if self.wb_areas is None:
                self.make_wholeblock(peakfit=False, bulk=True)
            y = self.wb_areas
        # Peak-specific 
        else:
            if self.peak_wb_areas is None:
                self.make_wholeblock(peakfit=True, bulk=False)
            y = self.peak_wb_areas[peak_idx]
            peakwn = self.peakpos[peak_idx]
            initial_unit_value = self.peak_initial_wb_ratio

        x = self.positions_microns
        L = self.set_len()
        t = self.time_seconds
        D0 = guess
        init = initial_unit_value
       
        params = lmfit.Parameters()
        #               (Name,           Value,      Vary,   Min,  Max,  Expr)
        params.add_many(('length_microns',  L,  False,  0.0,    None,   None),
                        ('time_seconds',    t,  False,  0.0,    None,   None),
                        ('log10D_m2s',   D0, True,   None,    None,   None),
                ('initial_unit_value', init, vary_initial, 0.0, None, None))
        kwdict = {'show_plot' : 'False'}
                    
        if show_plot is True:
            fig, ax = self.plot_area_profile_outline(peakwn=peakwn)
            self.plot_wholeblock_profile(figaxis=ax, peak_idx=peak_idx,
                                         top=top, 
                                         heights_instead=heights_instead)
        
        lmfit.minimize(diffusion_1D, params, args=(x, y), kws=kwdict)        
        best_D = ufloat(params['log10D_m2s'].value, 
                         params['log10D_m2s'].stderr)
        best_init = ufloat(params['initial_unit_value'].value, 
                         params['initial_unit_value'].stderr)


        if show_plot is True:
            diffusion_1D(params, fig_axis=ax, style={'color' : 'green'})
        
        print '\ntime in hours:', params['time_seconds'].value / 3600.
        print 'initial unit value:', best_init
        print 'bestfit log10D in m2/s:', best_D
        
        if peak_idx is not None:
            self.peak_diffusivities[peak_idx] = best_D.n
            self.peak_diff_error[peak_idx] = best_D.s
            if vary_initial is True:
                self.peak_initial_wb_ratio = best_init.n
        else:
            self.diffusivity_log10m2s = best_D.n
            self.diff_error = best_D.s
        
        return best_D, best_init

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

#%% 3D Plot setup
def plot_3panels_outline(style=None, top=1.2, figsize=(6.5, 2.5),
                         shrinker=0.1, heights_instead=False):
    """Outline setup for 3 subplots for 3D profiles"""
    if style is None:
        style = {'color' : 'lightgreen', 'linewidth' : 4}
        fig, axis3 = plt.subplots(nrows=1, ncols=3)
    fig.set_size_inches(figsize)
    for k in range(3):
        axis3[k].set_ylim(0, top)
        axis3[k].grid()
        box = axis3[k].get_position()
        axis3[k].set_position([box.x0 + box.width*shrinker, 
                               box.y0 + box.height*shrinker, 
                               box.width*(1.0-shrinker), 
                               box.height*(1.0-shrinker)])
    if heights_instead is False:
        axis3[0].set_ylabel('Area/Area$_0$')
    else:
        axis3[0].set_ylabel('Height/Height$_0$')
    axis3[1].set_xlabel('position ($\mu$m)')
    plt.setp(axis3[1].get_yticklabels(), visible=False)
    plt.setp(axis3[2].get_yticklabels(), visible=False)
    return fig, axis3
    
def plot_3panels(positions_microns, area_profiles, lengths,
                 styles3=[None, None, None], top=1.2, figaxis3=None, show_line_at_1=True,
                 centered=True, percent_error=3., xerror=50.,
                 heights_instead=False):
    """Make 3 subplots for 3D and 3DWB profiles. The position and area profiles
    are passed in lists of three lists for a, b, and c.
    Positions are assumed to start at 0 and then are centered.
    """
    if figaxis3 is None:
        fig, axis3 = plot_3panels_outline(top=top, 
                                          heights_instead=heights_instead)
    else:
        axis3 = figaxis3

    for k in range(3):      
        x = positions_microns[k]
        y = area_profiles[k]
        
        if len(x) != len(y):
            print 'Problem in plot_3panels'
            print 'len(x):', len(x)
            print 'len(y):', len(y)

        a = lengths[k] / 2.
        axis3[k].set_xlim(-a, a)

        if show_line_at_1 is True:
            axis3[k].plot([-a, a], [1., 1.], '-k')
    
        if styles3[k] is None:
            styles3[k] = style_profile
            
        if np.isnan(y).any():
            axis3[k].text(0, axis3[k].get_ylim()[1]/2., 
                         'nan values!\n\nProbably the\ninitial area was 0',
                         horizontalalignment='center', backgroundcolor='w',
                         verticalalignment='center')
        elif np.isinf(y).any():
            infstring = ''.join(('inf values!\n\nProbably the',
                                '\ninitial area was 0\nand a peak grew'))
            axis3[k].text(0, axis3[k].get_ylim()[1]/2., infstring,
                         horizontalalignment='center', backgroundcolor='w',
                         verticalalignment='center')
        else:
           yerror = np.array(y) * percent_error/100.
           axis3[k].errorbar(x-a, y, xerr=xerror, yerr=yerror, **styles3[k])

    if figaxis3 is None:
        return fig, axis3

#%% Generate 3D whole-block area and water profiles
def make_3DWB_area_profile(final_profile, 
                           initial_profile=None, 
                           initial_area_list=None, 
                           initial_area_positions_microns=None,
                           show_plot=True, top=1.2, fig_ax=None,
                           peakwn=None):
    """Take final 3D-wholeblock FTIR profile and returns
    a profile of the ratio of the two (A/Ao). A/A0 is also saved in the 
    profile attribute wb_areas.
    Requires information about initial state, so either an initial 
    profile (best) as the profile's initial_profile attribute,
    or else an initial area list passed in with its position. 
    Defaults to making a plot with max y-value 'top'.
    Note initial_area_positions_microns is assumed to start at 0 and then 
    gets shifted for the fit.
    """
    fin = final_profile
    leng = fin.set_len()

    # If whole-block areas are already made, use them. 
    # Otherwise make them.
    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
        wb_areas = fin.wb_areas
    else:
        print fin.wb_areas
        # initial checks
        if len(fin.positions_microns) == 0:
            print 'Need position information'
            return False
        if fin.len_microns is None:
            check = fin.set_len()
            if check is False:
                print 'Need more info to set profile length'
                return False
        if fin.len_microns is None:
            fin.set_len()
        if len(fin.areas_list) == 0:
            print 'making area list for profile'
            fin.make_area_list(show_plot=False)
    
        # What to normalize to? Priority given to self.wb_initial_profile, then
        # initial_profile passed in here, then initial area list passed in here.
        if fin.initial_profile is not None:
            
            init = fin.initial_profile
            if init.len_microns != fin.len_microns:
                print 'initial and final lengths must be the same!'
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas_list) == 0:
                    print profile.profile_name
                    print 'making area list for profile'
                    profile.make_area_list(show_plot=False)
            A0 = init.areas_list
            positions0 = init.positions_microns           
    
        elif initial_profile is not None:
            init = initial_profile
            if isinstance(init, Profile) is False:
                print 'initial_profile argument must be a pynams Profile.'
                print 'Consider using initial_area_list and positions instead'
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas_list) == 0:
                    print 'making area list for profile'
                    profile.make_area_list(show_plot=False)
            A0 = init.areas_list
            positions0 = init.positions_microns
    
        elif initial_area_list is not None:
            if initial_area_positions_microns is None:
                print 'Need initial_area_positions_microns for initial_area_list'
                return False
            A0 = initial_area_list
            positions0 = initial_area_positions_microns
            if len(fin.areas_list) == 0:
                print 'making area list for final profile'
                fin.make_area_list(show_plot=False)
        else:
            print 'Need some information about initial state'
            return False
    
        # More initial checks
        if len(fin.areas_list) != len(fin.positions_microns):
            print 'area and position lists do not match'
            print 'length areas_list:', len(fin.areas_list)
            print 'length positions list:', len(fin.positions_microns)
            return False    
        if len(A0) < 1:        
            print 'Nothing in initial area list'
            return False
        if len(positions0) < 1:
            print 'Nothing in initial positions list'
            return False
        if len(A0) == 1:
            print 'Using single point to generate initial line'
            A0.extend([A0[0], A0[0]])
            positions0.extend([0, fin.len_microns])
            
        # Use best-fit line through initial values to normalize final data
        p = np.polyfit(positions0-(leng/2.), A0, 1)
        
        normalizeto = np.polyval(p, fin.areas_list)
        wb_areas = fin.areas_list / normalizeto
         
        # Save whole-block areas as part of profile
        fin.wb_areas = wb_areas    
    
    if show_plot is True:
        if fig_ax is None:
            f, ax = final_profile.plot_area_profile_outline()
        else:
            ax = fig_ax
        ax.set_ylim(0, top)
        ylabelstring = 'Final area / Initial area'
        if peakwn is not None:
            print 'NOT READY FOR PEAKFITTING YET'
#            extrabit = '\n for peak at ' + str(peakwn) + ' cm$^{-1}$'
#            ylabelstring = ylabelstring + extrabit
        ax.set_ylabel(ylabelstring)
        style = fin.choose_marker_style()
        ax.plot([-leng/2.0, leng/2.0], [1, 1], **style_1)
        ax.plot(fin.positions_microns - (leng/2.0), wb_areas, **style)
        if fig_ax is None:
            return wb_areas, f, ax
        else:
            return wb_areas
    else:
        return wb_areas


def make_3DWB_water_profile(final_profile, water_ppmH2O_initial=None,
                            initial_profile=None, 
                            initial_area_list=None, 
                            initial_area_positions_microns=None,
                            show_plot=True, top=1.2, fig_ax=None):
    """Take a profile and initial water content.
    Returns the whole-block water concentration profile based on
    the profile's attribute wb_areas. If wb_areas have not been made, 
    some initial profile information and various options are passed
    to make_3DWB_area_profile().
    Default makes a plot showing A/Ao and water on parasite y-axis
    """
    fin = final_profile
    init = initial_profile

    # Set initial water
    if water_ppmH2O_initial is not None:
        w0 = water_ppmH2O_initial
    else:
        if fin.sample is not None:
            if fin.sample.initial_water is not None:
                w0  = fin.sample.initial_water
        elif init is not None:
            if init.sample is not None:
                if init.sample.initial_water is not None:
                    w0 = init.sample.initial_water
        else:
            print 'Need initial water content.'
            return False
    
    # Set whole-block areas
    if (fin.wb_areas is not None) and (len(fin.wb_areas) > 0):
        wb_areas = fin.wb_areas
    else:  
        wb_areas = make_3DWB_area_profile(fin, initial_profile, 
                                          initial_area_list, 
                                          initial_area_positions_microns)
    water = wb_areas * w0
    if show_plot is True:
        # Use a parasite y-axis to show water content
        fig = plt.figure()
        ax_areas = SubplotHost(fig, 1,1,1)
        fig.add_subplot(ax_areas)
        area_tick_marks = np.arange(0, 100, 0.2)
        ax_areas.set_yticks(area_tick_marks)
        ax_water = ax_areas.twin()
        ax_water.set_yticks(area_tick_marks)
        if isinstance(w0, uncertainties.Variable):
            ax_water.set_yticklabels(area_tick_marks*w0.n)
        else:
            ax_water.set_yticklabels(area_tick_marks*w0)
        ax_areas.axis["bottom"].set_label('Position ($\mu$m)')
        ax_areas.axis["left"].set_label('Final area / Initial area')
        ax_water.axis["right"].set_label('ppm H$_2$O')
        ax_water.axis["top"].major_ticklabels.set_visible(False)
        ax_water.axis["right"].major_ticklabels.set_visible(True)
        ax_areas.grid()
        ax_areas.set_ylim(0, 1.2)
        if fin.len_microns is not None:
            leng = fin.len_microns
        else:
            leng = fin.set_len()
        ax_areas.set_xlim(-leng/2.0, leng/2.0)
            
        style = fin.choose_marker_style()
        ax_areas.plot([-leng/2.0, leng/2.0], [1, 1], **style_1)
        ax_areas.plot(fin.positions_microns-leng/2.0, wb_areas, **style)
        return water, fig, ax_areas
    else:
        return water

#%% 
#
#
# Writing diffusion equations as suitable ojective functions that can be 
# called by lmfit.minimize.
#
#
def diffusion_1D(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', show_plot=True, fig_axis=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=50):
    """Function set up to follow lmfit fitting requirements.
    Requires input as lmfit parameters value dictionary 
    passing in key information as 'length_microns',
    'time_seconds', 'log10D_m2s', and 'initial_unit_value'. 
    
    Here is an example of the setup:
    params = lmfit.Parameters()
    #           (Name,                 Value,   Vary,   Min,  Max,  Expr)
    params.add_many(('length_microns',  1000,   False,  0.0, None,  None),
                ('time_seconds',      3600*10., False,  0.0,  None,  None),
                ('log10D_m2s',         -12.,    True,   None, None,  None),
                ('initial_unit_value',  1.,     False,  0.0,   1.0, None))

    Fitting is a profile function fitDiffusivity().
    For fitting, include data_x_microns and data_y_unit_areas. 
    If need_to_center_x_data is True (default), x data starts at 0 and
    will get moved. Also, lmfit screws up if you try setting a minimum
    on log10D_m2s.

    Return values
    If data are None (default), returns 1D unit diffusion profile 
    x_microns and y as vectors of length points (default 50). 

    Optional keywords:    
     - erf_or_sum: whether to use python's error functions (default) 
       or infinite sums
     - show_plot: whether to plot results (default True, so plot)
     - fig_ax: which figure axis to plot onto (default None makes new fig)
     - style: curve plotting style dictionary
#    - whether diffusion is in or out of sample (default out), 
#    - equilibrium concentration (default 0 for diffusion out; 1 for in), 
     - points sets how many points to calculate in profile. Default is 50.
     - what 'infinity' is if using infinite sum approximation
    """
    # extract important values from parameter dictionary passed in
    p = params.valuesdict()
    L_meters = p['length_microns'] / 1e6
    t = p['time_seconds']
    D = 10.**p['log10D_m2s']
    initial_value = p['initial_unit_value']
    a_meters = L_meters / 2.

    if t < 0:
        print 'no negative time'
        return           

    # Fitting to data or not? Default is not
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        if len(data_x_microns) == len(data_y_unit_areas):
            fitting = True
        else:
            print 'x and y data must be the same length'
            print 'x', len(data_x_microns)
            print 'y', len(data_y_unit_areas)
        
    # x is in meters and assumed centered around 0
    if fitting is True:
        # Change x to meters and center it
        x = np.array(data_x_microns) / 1e6
        if need_to_center_x_data is True:
            x = x - a_meters
    else:
        x = np.linspace(-a_meters, a_meters, points)
    
    if erf_or_sum == 'infsum':
        xsum = np.zeros_like(x)
        for n in range(infinity):
            xadd =  ((((-1.)**n) / ((2.*n)+1.))    + 
                    (np.exp((-D * (((2.*n)+1.)**2.) * (np.pi**2.) * t) / 
                        (L_meters**2.)))    + 
                    (np.cos(((2.*n)+1.) * np.pi * x / L_meters)))
            xsum = xsum + xadd        
        model = xsum * 4. / np.pi        

    elif erf_or_sum == 'erf':
        sqrtDt = (D*t)**0.5
        model = ((scipy.special.erf((a_meters+x)/(2*sqrtDt))) + 
                   (scipy.special.erf((a_meters-x)/(2*sqrtDt))) - 1) 

    else:
        print ('erf_or_sum must be set to either "erf" for python built-in ' +
               'error function approximation (defaul) or "sum" for infinite ' +
               'sum approximation with infinity=whatever, defaulting to ' + 
               str(infinity))
        return False

    # Revisit for in_or_out
    model = model * initial_value

    x_microns = x * 1e6

    if show_plot is True:
        a_microns = a_meters * 1e6
        if fig_axis is None:
            fig, fig_axis = plt.subplots()
            fig_axis.grid()
            fig_axis.set_ylim(0, 1.2)
            fig_axis.set_xlim(-a_microns, a_microns)
            fig_axis.set_ylabel('Area/Area$_0$')
            fig_axis.set_xlabel('position ($\mu$m)')

        if style is None:
            if fitting is True:
                style = {'linestyle' : 'none', 'marker' : 'o'}
            else:
                style = {'color' : 'lightgreen', 'linewidth' : 4}
            
        fig_axis.plot([-a_microns, a_microns], [initial_value, initial_value], 
                      '-k')
        fig_axis.plot(x_microns, model, **style)

    # If not including data, just return the model values
    # With data, return the residual for use in fitting.
    if fitting is False:
        return x_microns, model
    return model-data_y_unit_areas


def diffusion_3D(params, data_x_microns=None, data_y_unit_areas=None, 
                 erf_or_sum='erf', show_plot=True, fig_ax=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=50, show_1Dplots=False):
    """ Diffusion in 3 dimensions without path integration.
    lmfit parameter list input must include:
    length_microns_a, length_microns_b, length_microns_c, 
    time_seconds, initial unit value, and 
    diffusivities in log10 m2/s logDx, logDy, logDz
    
    Example parameter setup for input here:
    params = lmfit.Parameters()
    #  (Name,           Value,      Vary,   Min,  Max,  Expr)
    params.add('microns_twoa', L3[0], False, 0.0, None, None)
    params.add('microns_twob', L3[1], False, 0.0, None, None)
    params.add('microns_twoc', L3[2], False, 0.0, None, None)
    params.add('initial_unit_value_a', 1., False, 0.0, None, None)
    params.add('initial_unit_value_b', 1., True, 0.0, None, None)
    params.add('initial_unit_value_c', 1., False, 0.0, None, None)
    params.add('log10Dx', D3[0], True, 0.0, None, None)
    params.add('log10Dy', D3[1], True, 0.0, None, None)
    params.add('log10Dz', D3[2], False, 0.0, None, None)
    params.add('time_seconds', t, False, 0.0, None, None)    
    """
    # Fitting to data or not? Default is not
    # Add appropriate x and y data to fit
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_data = np.array(data_x_microns)
        y_data = np.array(data_y_unit_areas)
        if np.shape(x_data) == np.shape(y_data):
            fitting = True
            print 'fitting to data'
        else:
            print 'x and y data must be the same shape'
            print 'x', np.shape(x_data)
            print 'y', np.shape(y_data)

    p = params.valuesdict()
    L3_microns = np.array([p['microns_twoa'], p['microns_twob'], 
                          p['microns_twoc']])
    t = p['time_seconds']
    init = [p['initial_unit_value_a'], 
            p['initial_unit_value_b'],
            p['initial_unit_value_c']]
    vary_init = [params['initial_unit_value_a'].vary, 
                 params['initial_unit_value_b'].vary,
                 params['initial_unit_value_c'].vary]
    log10D3 = [p['log10Dx'], p['log10Dy'], p['log10Dz']]
    vary_D = [params['log10Dx'].vary, 
              params['log10Dy'].vary, 
              params['log10Dz'].vary]

    # First create 3 1D profiles, 1 in each direction
    xprofiles = []    
    yprofiles = []
    kwdict = {'show_plot' : show_1Dplots, 'points' : points}
    
    for k in range(3):
        p1D = lmfit.Parameters()
        p1D.add('length_microns', L3_microns[k], False)
        p1D.add('time_seconds', t, params['time_seconds'].vary)
        p1D.add('log10D_m2s', log10D3[k], vary_D[k])
        p1D.add('initial_unit_value', init[k], vary_init[k])
        
        x, y = diffusion_1D(p1D, **kwdict)
        xprofiles.append(x)
        yprofiles.append(y)
                                      
    # Then multiply them together to get a 3D matrix
    v = np.ones((points, points, points))
    for d in range(0, points):
        for e in range(0, points):
            for f in range(0, points):
                v[d][e][f] = yprofiles[0][d]*yprofiles[1][e]*yprofiles[2][f]

    mid = points/2        
    aslice = v[:, mid][:, mid]
    bslice = v[mid][:, mid]
    cslice = v[mid][mid]
    sliceprofiles = [aslice, bslice, cslice]
                
    if show_plot is True:
        if fig_ax is None:
            f, fig_ax = plot_3panels_outline()
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}
        positions = []
        for k in range(3):
            a = L3_microns[k] / 2.
            x = np.linspace(0, a*2., points)
            positions.append(x)
        plot_3panels(positions, sliceprofiles, L3_microns, style=style,
                     figaxis3=fig_ax)
            
    # Returning full matrix and 
    # slice profiles in one long list for use in fitting
    sliceprofiles = [aslice, bslice, cslice]
    
    if fitting is False:
        return v, sliceprofiles
    else:
        ### Still need to set up residuals! ###
        residuals = np.zeros_like(sliceprofiles)
        return residuals


def diffusion_3DWB(params, data_x_microns=None, data_y_unit_areas=None, 
                 raypaths=None, erf_or_sum='erf', show_plot=True, 
                 fig_ax=None,
                 style=None, need_to_center_x_data=True,
                 infinity=100, points=200, show_1Dplots=False):
    """ Diffusion in 3 dimensions with path integration.
    lmfit parameter list input must include:
    length_microns_a, length_microns_b, length_microns_c, 
    time_seconds, initial unit value, 
    diffusivities in log10 m2/s logDx, logDy, logDz, 
    
    Also must pass in a keyword list of raypaths in an order consistent
    with the length directions as, e.g., ['c', 'c', 'b']
    
    Example parameter setup for input here:
    params = lmfit.Parameters()
    #  (Name,           Value,      Vary,   Min,  Max,  Expr)
    params.add('microns_twoa', L3[0], False, 0.0, None, None)
    params.add('microns_twob', L3[1], False, 0.0, None, None)
    params.add('microns_twoc', L3[2], False, 0.0, None, None)
    params.add('initial_unit_value_a', 1., False, 0.0, None, None)
    params.add('initial_unit_value_b', 1., True, 0.0, None, None)
    params.add('initial_unit_value_c', 1., False, 0.0, None, None)
    params.add('log10Dx', D3[0], True, 0.0, None, None)
    params.add('log10Dy', D3[1], True, 0.0, None, None)
    params.add('log10Dz', D3[2], False, 0.0, None, None)
    params.add('time_seconds', t, False, 0.0, None, None)
    """
    if raypaths is None:
        print 'raypaths must be in the form of a list of three abc directions'
        return

    # v is the model 3D array of internal concentrations
    ### Need to add in all the keywords ###
    v, sliceprofiles = diffusion_3D(params, show_plot=False, points=points)

    # Fitting to data or not? Default is not
    # Add appropriate x and y data to fit
    fitting = False
    if (data_x_microns is not None) and (data_y_unit_areas is not None):
        x_array = np.array(data_x_microns)
        y_array = np.array(data_y_unit_areas)
        if np.shape(x_array) == np.shape(y_array):
            print 'fitting to data'
            fitting = True
        else:
            print 'x and y data must be the same shape'
            print 'x', np.shape(x_array)
            print 'y', np.shape(y_array)
            
    # Whole-block measurements can be obtained through any of the three 
    # planes of the whole-block, so profiles can come from one of two ray path
    # directions. These are the planes.
    raypathA = v.mean(axis=0)
    raypathB = v.mean(axis=1)
    raypathC = v.mean(axis=2)

    # Specify whole-block profiles in model
    mid = points/2
    if raypaths[0] == 'b':
        wbA = raypathB[:, mid]
    elif raypaths[0] == 'c':
        wbA = raypathC[:, mid]       
    else:
        print 'raypaths[0] for profile || a must be "b" or "c"'
        return
        
    if raypaths[1] == 'a':
        wbB = raypathA[:, mid]
    elif raypaths[1] == 'c':
        wbB = raypathC[mid]       
    else:
        print 'raypaths[1] for profile || b must be "a" or "c"'
        return

    if raypaths[2] == 'a':
        wbC = raypathA[mid]
    elif raypaths[2] == 'b':
        wbC = raypathB[mid]       
    else:
        print 'raypaths[2] for profile || c must be "a" or "b"'
        return

    p = params.valuesdict()
    L3 = [p['microns_twoa'], p['microns_twob'], p['microns_twoc']]
#    if fitting is True:
#        wb_profiles = data_y_unit_areas
#        wb_positions = np.array(data_x_microns)
#    else:
    wb_profiles = [wbA, wbB, wbC]
    wb_positions = []
    for k in range(3):
        a = L3[k] / 2.
        x_microns = np.linspace(0., 2.*a, points)
        wb_positions.append(x_microns)
        
    if show_plot is True:
        if style is None:
            style = {'color' : 'lightgreen', 'linewidth' : 4}

        if fig_ax is None:
            f, fig_ax = plot_3panels(wb_positions, wb_profiles, L3, style)
        else:
            plot_3panels(wb_positions, wb_profiles, L3, style, 
                         figaxis3=fig_ax)                         
    if fitting is False:        
        return wb_positions, wb_profiles
    
    if fitting is True:
        # Return residuals 
        y_model = []
        y_data = []
        residuals = []
        for k in range(3):
            for pos in range(len(x_array[k])):
                # wb_positions are centered, data are not
                microns = x_array[k][pos]
                # Find the index of the full model whole-block value 
                # closest to the data positions
                idx = (np.abs(wb_positions[k]-microns).argmin())
                
                model = wb_profiles[k][idx]
                data = y_array[k][pos]
                res = model - data
                
                y_model.append(model)
                y_data.append(data)
                residuals.append(res)                
        return residuals


                
#%% Group profiles together as whole-block unit
class WholeBlock():
    profiles = []
    name = None
    # generated by setupWB below
    directions = None
    raypaths = None
    initial_profiles = None
    lengths = None
    # optional for diffusion work and plotting
    style_base = None
    time_seconds = None
    diffusivities_log10_m2s = None
    diffusivity_errors = None
    peakfit_filename = None
    worksheetname = None
                
    def setupWB(self, peakfit=False, make_wb_areas=False):
        """Sets up and checks WholeBlock instance
        - Check that profiles list contains a list of three (3) profiles
        - Generate list of initial profiles
        - Generate list of profile directions
        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
        - Generate list of ray paths
        - Verify three ray path directions are compatible with directions list
        """
        if len(self.profiles) != 3:
            print 'For now, only a list of 3 profiles is allowed'
            return False
        d = []
        r = []
        ip = []
        L = []
        
        for prof in self.profiles:
            if isinstance(prof, Profile) is False:
                print 'Only profiles objects allowed in profile list!'
                return False
            d.append(prof.direction)
            r.append(prof.raypath)
            L.append(prof.set_len())

            if prof.positions_microns is None:
                print 'profile', prof.profile_name
                print 'Needs positions_microns attribute'

            if prof.initial_profile is None:
                print 'Need to set initial_profile attribute for each profile.'
                print 'Assuming these are initial profiles...'
                prof.initial_profile = prof
            ip.append(prof.initial_profile)

            
            if make_wb_areas is True:
                check = prof.make_wholeblock(peakfit=peakfit)
                if check is False:
                    return False

        self.directions = d
        self.raypaths = r
        self.initial_profiles = ip 
        self.lengths = L

        if peakfit is True:
            for prof in self.profiles:
                prof.get_peakfit()

        self.get_baselines()
            
        return True

    def get_peakfit(self):
        """Get peakfit information for all profiles"""
        for prof in self.profiles:
            prof.get_peakfit()

    def print_baseline_limits(self, initial_too=True):
        """Print out baseline wavenumber range for each profile"""
        for prof in self.profiles:
            print '\n', prof.profile_name
            print prof.spectra_list[0].base_low_wn, 
            print prof.spectra_list[0].base_high_wn
            if initial_too is True:
                print prof.initial_profile.profile_name
                print prof.initial_profile.spectra_list[0].base_low_wn, 
                print prof.initial_profile.spectra_list[0].base_high_wn


    def get_baselines(self, initial_too=False):
        """Get baselines for all spectra"""
        if self.initial_profiles is None:
            self.setupWB()
        for prof in self.profiles:
            for spectrum in prof.spectra_list:
                spectrum.get_baseline()
        if initial_too is True:
            for prof in self.initial_profiles:
                for spectrum in prof.spectra_list:
                    spectrum.get_baseline()
            

    def plot_spectra(self, profile_idx=None, show_baseline=True, 
                     show_initial_ave=True, 
                     show_final_ave=True, plot_all=False, 
                     initial_and_final_together=False, style=style_spectrum, 
                     stylei=style_initial, wn=None):
        """Plot all spectra in all or specified profile in whole-block"""
        if profile_idx is None:
            proflist = ([self.profiles[0]] + 
                        [self.profiles[1]] + 
                        [self.profiles[2]])
        else:
            proflist = [self.profiles[profile_idx]]
            
        for prof in proflist:
            print prof.profile_name
            prof.plot_spectra(show_baseline=show_baseline, 
                  show_initial_ave=show_initial_ave, plot_all=plot_all,
                  show_final_ave=show_final_ave,
                  initial_and_final_together=initial_and_final_together,
                  style=style, stylei=stylei, wn=None)

    def plot_3panels_ave_spectra(self, peak_idx=None, peakwn=None, top=1.,
                                style=style_spectrum, stylei=style_initial,
                                figsize=(6., 2.5), show_initial=True,
                                legloc=5):
        """Three suplots showing average initial and final spectra in each
        direction"""
        if self.initial_profiles is None:
            self.setupWB()

        f, ax3 = plt.subplots(1,3)
        f.set_size_inches(figsize)
        for k in range(3):
            prof = self.profiles[k]
            avespec = prof.average_spectra()
            ax3[k].plot(avespec.wn, avespec.abs_cm, label='Final', **style)

            if peak_idx is not None:
                if self.profiles[k].peakpos is None:
                    self.profiles[k].get_peakfit()
                peakpos = self.profiles[k].peakpos
                peakwn = peakpos[peak_idx]
                ax3[k].plot([peakwn, peakwn], [0, top], color='r')

            if show_initial is True:
                iprof = self.initial_profiles[k]
                if iprof is not None:
                    initspec = iprof.average_spectra()
                    ax3[k].plot(initspec.wn, initspec.abs_cm, label='Initial', 
                            **stylei)

            ax3[k].set_xlim(4000., 3000.)
            ax3[k].set_ylim(0., top)
            
            raypathstring = 'profile || ' + self.profiles[k].direction + \
                            '\nray path || ' + self.profiles[k].raypath
            ax3[k].text(3500, top-top*0.22, raypathstring, 
                        backgroundcolor='w', horizontalalignment='center')

        plt.setp(ax3[1].get_yticklabels(), visible=False)
        plt.setp(ax3[2].get_yticklabels(), visible=False)
        ax3[1].set_xlabel('wavenumber (cm$^{-1}$)')
        ax3[0].set_ylabel('absorbance m$^{-1}$)')
        
        tit = ' '.join(('Averaged profiles for', self.name))
        if peak_idx is not None:
            tit = str.join(tit, ', peak at ', str(peakpos[peak_idx]), '/cm')
        ax3[1].text(3500, top+0.07*top, tit, horizontalalignment='center')

        if show_initial is True:
            ax3[1].legend(loc=legloc)
        plt.tight_layout()
        f.autofmt_xdate()
        plt.subplots_adjust(top=0.85, bottom=0.25)
        return f, ax3

    def plot_wholeblock(self, peak_idx=None, peakwn=None, fig_ax3=None, 
                        top=1.2, f=None, wn=None, figsize=(6.5, 3.5), 
                        sshow_spectra=True, percent_error=3., 
                        styles3=[None, None, None],
                        use_area_profile_styles=True,
                        heights_instead=False):
        """Plot whole-block ratio of Area/Initial Area on three panels"""
        if ((self.directions is None) or (self.raypaths is None) or
            (self.initial_profiles is None) or (self.lengths is None)):
                check = self.setupWB(make_wb_areas=False, peakfit=False)
                if check is False:
                    print 'Problem setting up whole block in setupWB'                    
                    return False
        self.get_baselines()
        # concatenate positions and areas across three profiles to 
        # send in to the plotting function
        positions = []
        areas = []        
        for prof in self.profiles:
            positions.append(prof.positions_microns)

            # Bulk hydrogen            
            if peak_idx is None and peakwn is None:
                
                if prof.wb_areas is None:
                    print '\nMaking whole block area ratios'
                    check = prof.make_wholeblock(peakfit=False, bulk=True,
                                                 show_plot=False)
                    if check is False:
                        return False
    
                areas.append(prof.wb_areas)
                
                if max(prof.wb_areas > top):
                    top = max(prof.wb_areas) + 0.1*max(prof.wb_areas)

            # Peak-specific                
            else:
                peak_wb_areas, peakwn = prof.get_peak_wb_areas(peak_idx, 
                                           peakwn, 
                                           heights_instead=heights_instead)
                    
                peakpos = self.profiles[0].peakpos
                
                areas.append(peak_wb_areas)
                    
                if np.isinf(peak_wb_areas).any():
                    pass
                elif max(peak_wb_areas > top):
                    top = max(peak_wb_areas) + 0.1*max(peak_wb_areas)

        if use_area_profile_styles is True:
            for k in range(3):
                styles3[k] = self.profiles[k].choose_marker_style()
                styles3[k]['markersize'] = 10

        # Sent positions and areas to plotting command above
        if fig_ax3 is not None:
            plot_3panels(positions, areas, self.lengths, figaxis3=fig_ax3,
                         styles3=styles3, top=top, 
                         heights_instead=heights_instead)
        else:
            f, fig_ax3 = plot_3panels(positions, areas, self.lengths,
                         styles3=styles3, top=top, 
                         heights_instead=heights_instead)

        # Change title if peak-specific rather than bulk
        if peak_idx is not None:
            tit = ' '.join(('Peak at', str(peakpos[peak_idx]), '/cm'))
        else:
            tit = 'Bulk hydrogen'
        fig_ax3[1].text(0, top+0.07*top, tit, horizontalalignment='center')
        
        f.set_size_inches(figsize)
        f.autofmt_xdate()
#        plt.subplots_adjust(top=0.9, bottom=0.35)
        return f, fig_ax3

    def print_spectra_names(self, show_initials=True):
        """Print out fnames of all spectra associated with the whole-block 
        instance."""
        if show_initials is True:
            if self.initial_profiles is None:
                self.setupWB()
                
            if self.initial_profiles is not None:
                print '--Initial profiles--'
                for prof in self.initial_profiles:
                    print prof.profile_name
                    spec_list = []
                    for spectrum in prof.spectra_list:
                        spec_list.append(spectrum.fname)
                    print spec_list
                    print ' '
            else:
                print 'No initial profiles given'

        if self.profiles is not None:

            print '--Final profiles--'
            for prof in self.profiles:
                print prof.profile_name
                spec_list = []
                for spectrum in prof.spectra_list:
                    spec_list.append(spectrum.fname)
                print spec_list
                print ' '

    def make_profile_list(self, initial_too=False):
        """Return False or a list of profiles"""
        if initial_too is True:
            if self.initial_profiles is None:
                self.setupWB()
            if self.initial_profiles is None:
                print 'No initial profiles'
                profile_list = self.profiles
            else:
                profile_list = self.initial_profiles + self.profiles
        else:
            profile_list = self.profiles
        return profile_list

    def make_baselines(self, initial_too=False, line_order=1, shiftline=None, 
                       show_fit_values=False, show_plot=False): 
        """Make spectra baselines for all spectra in all profiles in 
        whole-block"""
        profile_list = self.make_profile_list(initial_too)        
        for prof in profile_list:
            for spectrum in prof.spectra_list:
                spectrum.make_baseline(line_order, shiftline, show_fit_values,
                                       show_plot)                

    def save_baselines(self, initial_too=True):
        """Make and save spectra baselines for all spectra."""
        profile_list = self.make_profile_list(initial_too)
        for prof in profile_list:
            for spectrum in prof.spectra_list:
                spectrum.save_baseline()

    def print_spectra4matlab(self, initial_too=False):
        """Print out a list of spectra names in a matlab-friendly way"""        
        print '\nFor use in FTIR_peakfit_loop.m\n'        
        
        if initial_too is True:
            if self.initial_profiles is None:
                self.setupWB()
            string = "{"
            for prof in self.initial_profiles:
                for spec in prof.spectra_list:
                    stringname = spec.fname
                    string = string + "'" + stringname + "' "
            string = string + "};"
            print string, '\n'
        
        string = "{"
        for prof in self.profiles:
            for spec in prof.spectra_list:
                stringname = spec.fname
                string = string + "'" + stringname + "' "
        string = string + "}"
        print string

    def plot_areas(self, profile_index=None, peak_idx=None, peakwn=None, 
                   show_initials=False, show_finals=True, show_legend=True,
                   legloc=1, frame=False, top=None, bestfitlines=False,
                   heights_instead=False):
        """Plots profiles on one figure.
        Set initial_instead_of_final to True to see initial areas.
        Need to add legend and checks is wavenumber not in list."""

        # Which profiles to plot
        if show_initials is True:
            if self.initial_profiles is None:
                self.setupWB()
            if self.initial_profiles is None:
                print 'Need initial profile'
                return
            if self.initial_profiles is None:
                self.setupWB(False, False)

        # get wavenumber if only peak_idx is givin
        if peak_idx is not None:
            if profile_index is not None:
                idx = profile_index
            else:
                idx = 0            
            prof = self.profiles[idx]
            if prof.peakpos is None:
                prof.get_peakfit()
            if peak_idx is None:
                peak_idx = np.where(prof.peakpos==peakwn)[0][0]

        f, ax = self.profiles[0].plot_area_profile_outline(peakwn=peakwn)

        if profile_index is None:
            if show_finals is True:
                ai = self.profiles[0]
                bi = self.profiles[1]
                ci = self.profiles[2]
                for prof in [ai, bi, ci]:
                    prof.plot_area_profile(figaxis=ax, peakwn=peakwn, 
                                           peak_idx=peak_idx,
                                           bestfitline=bestfitlines,
                                           heights_instead=heights_instead)

            if show_initials is True:
                ai = self.initial_profiles[0]
                bi = self.initial_profiles[1]
                ci = self.initial_profiles[2]            
                ai.plot_area_profile(figaxis=ax, peakwn=peakwn, 
                                     peak_idx=peak_idx,
                                     heights_instead=heights_instead)
                bi.plot_area_profile(figaxis=ax, peakwn=peakwn, 
                                     peak_idx=peak_idx,
                                     heights_instead=heights_instead)
                ci.plot_area_profile(figaxis=ax, peakwn=peakwn, 
                                     peak_idx=peak_idx,
                                     heights_instead=heights_instead)
        else:
            if show_finals is True:
                self.profiles[profile_index].plot_area_profile(
                        figaxis=ax, peakwn=peakwn, peak_idx=peak_idx,
                        heights_instead=heights_instead)
            if show_initials is True:
                self.initial_profiles[profile_index].plot_area_profile(
                        figaxis=ax, peakwn=peakwn, peak_idx=peak_idx,
                        heights_instead=heights_instead)

        if show_legend is True:
            leg_handle_list = []
            descript = ['profile || a*', 'raypath || a*', 'profile || b', 
                        'raypath || b',  'profile || c', 'raypath || c']
            
            bstylelinebase = {'marker' : 's', 'color' : 'black', 'alpha' : 0.5,
                         'markersize' : 10, 'linestyle': 'none'}
            bstyleline = [None, None, None, None, None, None]
            bstyleline[0] = dict(bstylelinebase.items() + style_Dx.items())
            bstyleline[1] = dict(bstylelinebase.items() + style_Rx.items())
            bstyleline[2] = dict(bstylelinebase.items() + style_Dy.items())
            bstyleline[3] = dict(bstylelinebase.items() + style_Ry.items())
            bstyleline[4] = dict(bstylelinebase.items() + style_Dz.items())
            bstyleline[5] = dict(bstylelinebase.items() + style_Rz.items()) 
            
            for k in range(6):
                add_marker = mlines.Line2D([], [], label=descript[k], 
                                           **bstyleline[k])
                leg_handle_list.append(add_marker)
            
            # Shrink current axis's height by 10% on the bottom
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1,
                             box.width, box.height * 0.9])
            
            # Put a legend below current axis
            ax.legend(handles=leg_handle_list,
                      loc='upper center', bbox_to_anchor=(0.5, -0.3),
                      fancybox=True, ncol=3)

        all_areas = np.array([])
        if top is not None:
            ax.set_ylim(0, top)
        else:
            top = ax.get_ylim()[1]
            if peak_idx is None:
                for prof in [ai, bi, ci]:
                    all_areas = np.append(all_areas, prof.make_area_list())
            else:
                for prof in [ai, bi, ci]:
                    if heights_instead is False:
                        all_areas = np.append(all_areas, 
                                              prof.peak_areas[peak_idx])
                    else:
                        all_areas = np.append(all_areas, 
                                              prof.peak_heights[peak_idx])
                                              
            maxtop = np.max(all_areas)
            if maxtop > top:
                top = maxtop + 0.15*maxtop
        ax.set_ylim(0, top)
                    
        if peak_idx is None:
            tit = 'Bulk hydrogen'
        else:
            tit = ' '.join(('Peak at',
                            str(self.profiles[idx].peakpos[peak_idx]),'/cm'))
        ax.text(0, top + 0.05*top, tit, horizontalalignment='center')
                   
        f.set_size_inches((6.5, 3.5))
        plt.subplots_adjust(top=0.9, bottom=0.35)
                   
        return f, ax

    def plot_peakfits(self, initial_too=False, profile_idx=None, legloc=1):
        """Plot peakfits for all spectra in all profiles"""
        if profile_idx is None:
            for prof in self.profiles:
                prof.plot_peakfits(initial_too, legloc=legloc)
        else: 
            self.profiles[profile_idx].plot_peakfits(initial_too, 
                                                    legloc=legloc)


    def excelify(self, filename=None, exceldpi=150, top=1.0, 
                 bestfitlines_on_areas=True):
        """Save all peak fit info AND peak figures to an excel spreadsheet """
        # file setup
        if filename is None and self.worksheetname is None:
            print 'Need filename information here'
            print 'or in worksheetname attribute of whole block instance'
            return
        if filename is None:
            filename = ''.join((self.worksheetname, '.xlsx'))            
        if self.profiles[0].peak_wb_areas is None:
            self.setupWB(peakfit=True, make_wb_areas=True)
        workbook = xlsxwriter.Workbook(filename)
        worksheetname = self.worksheetname        
        worksheet = workbook.add_worksheet(worksheetname)

        # Column locations
        col1 = 0
        col2 = 4
        col2pt5 = 8
        col3 = 12
        col4 = 20
        col_heights = 28
        col_wb_heights = 36
        col5 = 44

        # formats        
        worksheet.set_column(col1, col1, width=13)
        worksheet.set_column(col5, col5+50, width=15)
        worksheet.set_column(col1+1, col2, width=11)
        
        wraptext = workbook.add_format()
        wraptext.set_text_wrap()
        wraptext.set_align('center')
        wraptext.set_bold()
        
        boldtext = workbook.add_format()
        boldtext.set_bold()

        italictext = workbook.add_format()
        italictext.set_italic()

        one_digit_after_decimal = workbook.add_format()
        one_digit_after_decimal.set_num_format('0.0')
        one_digit_after_decimal.set_italic()

        two_digits_after_decimal = workbook.add_format()
        two_digits_after_decimal.set_num_format('0.00')

        highlight = workbook.add_format()
        highlight.set_bg_color('yellow')


        peakpos = self.profiles[0].spectra_list[0].peakpos
        
        worksheet.write(0, 0, self.name)
        worksheet.write(1, 0, 'peak position (/cm)', wraptext)
        worksheet.write(1, 1, 'peak height (/cm)', wraptext)
        worksheet.write(1, 2, 'peak width (/cm)', wraptext)
        worksheet.write(1, 3, 'peak area (/cm2)', wraptext)
        worksheet.write(1, 5, 'Baselines', wraptext)        
        worksheet.write(1, 9, 'Baseline-subtracted', wraptext)
        worksheet.set_column(9, 9, width=11)
        
        # Bulk area and peak fits in 1st column for individual fits
        row = 2        
        pic_counter = 0
        for prof in self.profiles:
            row = row + 1
            worksheet.write(row, col1, prof.profile_name)
            row = row + 1

            pos_idx = 0
            for spec in prof.spectra_list:                
                row = row + 1
                worksheet.write(row, 0, spec.fname)                               
                worksheet.write(row, 1, prof.positions_microns[pos_idx],
                                one_digit_after_decimal)
                worksheet.write(row, 2, 'microns from edge', italictext)                   
                
                # Show baselines in 2nd column
                spec.get_baseline()
                f, ax = spec.plot_showbaseline()
                newfig = 'baseline{:d}.png'.format(pic_counter)
                pic_counter = pic_counter + 1
                f.savefig(newfig, dpi=exceldpi)
                worksheet.insert_image(row, col2, newfig, 
                                       {'x_scale' : 0.5, 'y_scale' : 0.5})

                # Peak fit figures in 2nd column
                f, ax = spec.plot_peakfit()
                newfig = 'peakfit{:d}.png'.format(pic_counter)
                pic_counter = pic_counter + 1
                f.savefig(newfig, dpi=exceldpi)
                worksheet.insert_image(row, col2pt5, newfig, 
                                       {'x_scale' : 0.5, 'y_scale' : 0.5})

                pos_idx = pos_idx + 1
                row = row + 1

                if spec.peakpos is None:
                    check = spec.get_peakfit()
                    if check is False:
                        print 'trouble with getting peakfit'
                sumarea = 0                    
                for k in range(len(spec.peakpos)):
                    worksheet.write(row, 0, spec.peakpos[k])
                    worksheet.write(row, 1, spec.peak_heights[k])
                    worksheet.write(row, 2, spec.peak_widths[k])
                    worksheet.write(row, 3, spec.peak_areas[k])                   
                    
                    sumarea = sumarea + spec.peak_areas[k]
                    row = row + 1

                worksheet.write(row, 0, 'Sum of peaks areas')
                worksheet.write(row, 3, sumarea, two_digits_after_decimal)
                row = row + 1

                worksheet.write(row, 0, 'Observed total area')
                worksheet.write(row, 3, spec.area, two_digits_after_decimal)
                row = row + 1

        # 3 panel averaged spectra in 3rd column
        f, ax = self.plot_3panels_ave_spectra(top=top)
        f.savefig('panel3.png', dpi=exceldpi)
        worksheet.insert_image(0, col3, 'panel3.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
            
        # area profiles in 3rd column and whole-block profiles in 4th column
        # Then same thing for heights
        worksheet.write(11, col3, 'Peak areas profiles', boldtext)
        worksheet.write(11, col4, 'Whole-block area profiles', boldtext)
        worksheet.write(11, col_heights, 'Peak height profiles', boldtext)
        worksheet.write(11, col_wb_heights, 'Whole-block height profiles', 
                        boldtext)
        
        errorstring1 = ', '.join(('+/- 2% errors in area', 
                                  '+/- 50 microns errors in position'))
        errorstring2 = ', '.join(('+/- 3% errors in whole-block area ratio', 
                                  '+/- 50 microns errors in position'))
        errorstring3 = ', '.join(('assuming +/- 2% errors in heights', 
                                  '+/- 50 microns errors in position'))
        errorstring4 = ', '.join(('+/- 3% errors in whole-block height ratio', 
                                  '+/- 50 microns errors in position'))
        worksheet.write(12, col3, errorstring1)
        worksheet.write(12, col4, errorstring2)
        worksheet.write(12, col_heights, errorstring3)
        worksheet.write(12, col_wb_heights, errorstring4)

        worksheet.write(13, col3, 'Errors typically plot in or near symbols')
        worksheet.write(13, col4, 'Errors typically plot in or near symbols')
        worksheet.write(13, col_heights, 
                                'Errors typically plot in or near symbols')
        worksheet.write(13, col_wb_heights, 
                                'Errors typically plot in or near symbols')
        
        f, ax = self.plot_areas(bestfitlines=bestfitlines_on_areas)
        f.savefig('bulk_areas.png', dpi=exceldpi)
        worksheet.insert_image(15, col3, 'bulk_areas.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
        
        f, ax = self.plot_wholeblock()
        f.savefig('wbbulk.png', dpi=exceldpi)
        worksheet.insert_image(15, col4, 'wbbulk.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
                       
        for peak_idx in range(len(self.profiles[0].spectra_list[0].peakpos)):
            f, ax = self.plot_areas(bestfitlines=bestfitlines_on_areas, 
                                    peak_idx=peak_idx)
            newfig = 'area{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi)
            worksheet.insert_image(29+(14*peak_idx), col3, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            

            f, ax = self.plot_wholeblock(peak_idx=peak_idx)
            newfig = 'wbareas{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi)
            worksheet.insert_image(29+(14*peak_idx), col4, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})
                   
            f, ax = self.plot_areas(bestfitlines=bestfitlines_on_areas, 
                                    peak_idx=peak_idx, heights_instead=True)
            newfig = 'height{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi)
            worksheet.insert_image(29+(14*peak_idx), col_heights, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            

            f, ax = self.plot_wholeblock(peak_idx=peak_idx, 
                                         heights_instead=True)
            newfig = 'wbareas{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi)
            worksheet.insert_image(29+(14*peak_idx), col4, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})

        # Initial profiles and baseline information in 4rd column
        row = 1
        worksheet.write(row, col4, 'Initial profiles', boldtext) 
        worksheet.write(row+1, col4, self.initial_profiles[0].profile_name)
        worksheet.write(row+2, col4, self.initial_profiles[1].profile_name)
        worksheet.write(row+3, col4, self.initial_profiles[2].profile_name)   
        
        worksheet.write(row+5, col4, 'Baseline wavenumber ranges (/cm)', 
                        boldtext)
        for k in range(3):
            prof = self.profiles[k]
            iprof = self.initial_profiles[k]
            spec = prof.spectra_list[0]
            ispec = iprof.spectra_list[0]
            worksheet.write(row+6+k, col4, 
                            ''.join(('final || ', prof.direction)))
            worksheet.write(row+6+k, col4+3, 
                            ''.join(('initial || ',iprof.direction)))
            
            if spec.base_high_wn != ispec.base_high_wn:
                worksheet.write(row+6+k, col4+1, spec.base_high_wn, highlight)
                worksheet.write(row+6+k, col4+4, ispec.base_high_wn, highlight)
            else:
                worksheet.write(row+6+k, col4+1, spec.base_high_wn)
                worksheet.write(row+6+k, col4+4, ispec.base_high_wn)
            
            if spec.base_low_wn != ispec.base_low_wn:
                worksheet.write(row+6+k, col4+2, spec.base_low_wn, highlight)
                worksheet.write(row+6+k, col4+5, ispec.base_low_wn, highlight)
            else:
                worksheet.write(row+6+k, col4+2, spec.base_low_wn)
                worksheet.write(row+6+k, col4+5, ispec.base_low_wn)


        ### Summary list at the end - positions, areas, whole-block areas
        # labels
        worksheet.write(1, col5, 'centered peak position (/cm)', wraptext)
        worksheet.write(1, col5+1, 'bulk area (/cm2)', wraptext)
        col = col5 + 2
        for peak in peakpos:
            label1 = ''.join((str(peak), ' area (/cm2)'))
            worksheet.write(1, col, label1, wraptext)
            col = col + 1
        worksheet.write(1, col, 'bulk whole-block ratio (/cm2)', wraptext)
        col = col + 1
        for peak in peakpos:
            label2 = ''.join((str(peak), ' whole-block area ratio'))
            worksheet.write(1, col, label2, wraptext)
            col = col + 1
        worksheet.write(1, col, 'file label', wraptext)
        col = col + 1

        # values filling in the rows
        row = 2
        for prof in self.profiles:
            col = col5
  
            # label profile name with spaces on either side
            worksheet.write(row, col, prof.profile_name)            
            halflen = prof.set_len() / 2.            
            pos_idx = 0
            row = row + 1

            for spec in prof.spectra_list:
                col = col5 # Restart at first column each time
                worksheet.write(row, col, 
                                prof.positions_microns[pos_idx]-halflen,
                                two_digits_after_decimal)
                col = col + 1

                # bulk and peak fit areas
                worksheet.write(row, col, prof.areas_list[pos_idx],
                                two_digits_after_decimal)

                col = col + 1
                for k in range(len(peakpos)):
                    area = prof.peak_areas[k][pos_idx]
                    worksheet.write(row, col, area, 
                                    two_digits_after_decimal)
                    col = col + 1

                # whole-block bulk and peakfit ratios
                worksheet.write(row, col, prof.wb_areas[pos_idx],
                                two_digits_after_decimal)

                col = col + 1
                for k in range(len(peakpos)):
                    wb = prof.peak_wb_areas[k][pos_idx]
                    if np.isnan(wb):                        
                        worksheet.write(row, col, 'nan')
                    elif np.isinf(wb): 
                        worksheet.write(row, col, 'inf')
                    else:
                        worksheet.write(row, col, wb, 
                                    two_digits_after_decimal)
                    col = col + 1

                # filenames at the end of summary
                worksheet.write(row, col, spec.fname)

                row = row + 1
                pos_idx = pos_idx + 1

        workbook.close()        

    def print_peakfits(self, initial_too=False, excelfriendly=True):
        """Print out all peakfit information"""
        if initial_too is True and self.initial_profiles is not None:
            proflist = self.profiles + self.initial_profiles
        else:
            proflist = self.profiles

        if excelfriendly is True:
            print 'position height width area'
            
        for prof in proflist:
            print '\n\n', prof.profile_name
            poscounter = 0
            for spectrum in prof.spectra_list:
                print '\n', spectrum.fname, \
                        prof.positions_microns[poscounter], 'microns'
                poscounter += 1
                if spectrum.peakpos is None:
                    check = spectrum.get_peakfit()
                    if check is False:
                        print 'trouble with getting peakfit'
                if excelfriendly is False:
                    print 'position, height, width, area'
                if spectrum.peakpos is not None:
                    for k in range(len(spectrum.peakpos)):
                        print spectrum.peakpos[k], spectrum.peak_heights[k], \
                              spectrum.peak_widths[k], spectrum.peak_areas[k]        

    def print_diffusivities(self, peak_idx=None, profile_idx=None,
                            show_plot=True, top=1.0):
        """Print diffusivities for each profile"""
        if show_plot is True:
            self.plot_3panels_ave_spectra(peak_idx=peak_idx, top=top)
                     
        if peak_idx is None and profile_idx is None:
            for prof in self.profiles:
                if prof.peakpos is None:
                    prof.get_peakfit()
                print '\n', prof.profile_name, 'log10(D in m2/s)'
                print 'bulk H :  ', prof.diffusivity_log10m2s, '+/-', \
                        prof.diff_error
                for k in range(len(prof.peakpos)):
                    print prof.peakpos[k], ':  ', prof.peak_diffusivities[k],\
                            '+/-', prof.peak_diff_error[k]
            return
        
        if profile_idx is None:
            for prof in self.profiles:
                k = peak_idx
                print '\n', prof.profile_name
                print 'peak position and log10(diffusivity in m2/s)'
                print 'bulk H :', prof.diffusivity_log10m2s, '+/-', \
                    prof.diff_error
                print prof.peakpos[k], ':', prof.peak_diffusivities[k],\
                    '+/-', prof.peak_diff_error[k]
            return
            
        prof = self.profiles[profile_idx]
        print '\n', prof.profile_name
        print 'peak positions, log10(D in m2/s), D errors'
        print 'bulk H :  ', prof.diffusivity_log10m2s, '+/-', \
                prof.diff_error
        print prof.peakpos
        print prof.peak_diffusivities
        print prof.peak_diff_error

    def plot_diffusion(self, peak_idx=None, peakwn=None, 
                       time_seconds=None, 
                       list_of_log10D_m2s=None, 
                       initial_value=None, in_or_out='out', erf_or_sum='erf', 
                       equilibrium_value=None, show_plot=True, 
                       show_slice=False, style=None, wb_or_3Dnpi='wb', 
                       fig_ax=None, points=50, top_spectra=1.0,
                        top=1.2, numformat='{:.1f}', show_spectra=True,
                        heights_instead=False):
        """Applies 3-dimensionsal diffusion equations to instance shape.
        Requires lengths, time in seconds, and three diffusivities either
        explicitly passed here or as attributes of the WholeBlock object.
        Assuming whole-block diffusion (wb_or_3Dnpi='wb') but could also 
        do 3D non-path-integrated ('npi')
        """
        if time_seconds is None:
            if self.time_seconds is not None:
                time_seconds = self.time_seconds
            else:
                print 'Need time information'
                return        
        
        ## Which diffusivities to use
        if list_of_log10D_m2s is not None:
            D3 = list_of_log10D_m2s
        elif peak_idx is not None:
            if self.profiles[0].peak_diffusivities is not None:
                D3 = []
                for k in range(3):
                    D3.append(self.profiles[k].peak_diffusivities[peak_idx])
            else:
                D3 = [-13., -13., -13.]
        elif self.profiles[0].diffusivity_log10m2s is not None:
            D3 = []
            for k in range(3):
                D3.append(self.profiles[k].diffusivity_log10m2s)
        else:
            D3 = [-13., -13., -13.]
        
        if (wb_or_3Dnpi != 'npi') and (wb_or_3Dnpi != 'wb'):
            print 'wb_or_3Dnpi only takes "wb" or "npi"'
            return

        if self.directions is None:           
            self.setupWB(peakfit=False, make_wb_areas=False)

        if self.initial_profiles is None:
            self.setupWB(peakfit=False, make_wb_areas=False)

        if self.lengths is None:
            self.setupWB(peakfit=False, make_wb_areas=False)            

        if self.raypaths is None:
            self.setupWB(peakfit=False, make_wb_areas=False)
        # end initial checks and setup

        if peak_idx is not None:
            init3 = []
            for k in range(3):
                initwb = self.profiles[k].peak_initial_wb_ratio
                init3.append(initwb)
                print 'initial ratio', initwb
        else:
            init3 = [1., 1., 1.]

        # Set up parameters to pass into equations
        L3 = self.lengths
        params = lmfit.Parameters()
        # (Name, Value, Vary, Min, Max, Expr)
        params.add('microns_twoa', L3[0], False, None, None, None)
        params.add('microns_twob', L3[1], False, None, None, None)
        params.add('microns_twoc', L3[2], False, None, None, None)
        params.add('initial_unit_value_a', init3[0], False, None, None, None)
        params.add('initial_unit_value_b', init3[1], False, None, None, None)
        params.add('initial_unit_value_c', init3[2], False, None, None, None)
        params.add('log10Dx', D3[0], True, None, None, None)
        params.add('log10Dy', D3[1], True, None, None, None)
        params.add('log10Dz', D3[2], True, None, None, None)
        params.add('time_seconds', time_seconds, False, None, None, None)

        if show_spectra is True:
            if peak_idx is not None:
                self.plot_3panels_ave_spectra(peak_idx=peak_idx)

        if show_plot is True:
            if fig_ax is None:
                fig, fig_ax = plot_3panels_outline(top=top)
            
            kws = {'show_plot' : show_plot, 'fig_ax' : fig_ax,
                   'style' : style, 'points' : points}
    
            # send to diffusion equation
            # diffusion line plots get made here
            if init3 == [1., 1., 1.]:
                if wb_or_3Dnpi == 'npi':
                    diffusion_3D(params, **kws)
                else:
                    R3 = self.raypaths
                    kws['raypaths'] = R3
                    diffusion_3DWB(params, **kws)
            else:
                for prof in self.profiles:
                    prof.plot_wholeblock_profile(
                    heights_instead=heights_instead)
    
            # Add data on top of diffusion line in plot
            self.plot_wholeblock(fig_ax3=fig_ax, peak_idx=peak_idx, 
                              peakwn=peakwn)
            for k in range(3):
                dlabel = str(numformat.format(D3[k])) + ' m$^2$s$^{-1}$'
                fig_ax[k].text(0, top-top*0.12, dlabel, 
                                horizontalalignment='center',
                                backgroundcolor='w')                              
        return params
       
    def fitD(self, peak_idx=None, peakwn=None,
             initial_unit_values=[1., 1., 1.], 
             vary_initials=['False', 'False', 'False'], 
             guesses=[-13., -13., -13.], show_plot=True,
             vary_D=['True', 'True', 'True'],
             approx_with_1D=False, wb_or_3Dnpi='wb', polyorder=1,
             show_initial_guess=True, style_initial=None,
             style_final={'color' : 'red'}, points=100):
        """Forward modeling to determine diffusivities in three dimensions 
        from whole-block data. 
        """
        if wb_or_3Dnpi != 'wb' and wb_or_3Dnpi != 'npi':
            print 'wb_or_3Dnpi can only be wb or npi'
            return
            
        if wb_or_3Dnpi == 'npi':
            print 'NPI not supported yet'
            return
            
        # Plot setup and initial guess lines
        if show_plot is True:
            fig, ax3 = plot_3panels_outline()
        else: 
            ax3 = None

        # Set up parameters and keywords to go into diffusion equations
        # This will plot the initial guess diffusivities if show_initial_guess
        # is True
        dict_plotting = {'show_plot' : show_initial_guess, 'fig_ax' : ax3,
                 'wb_or_3Dnpi' : wb_or_3Dnpi, 'style' : style_initial,
                 'points' : points, 'peakwn' : peakwn, 'peak_idx' : peak_idx}
        dict_fitting = {'show_plot' : False, 'points' : points}

        params = self.plot_diffusion(list_of_log10D_m2s=guesses,
                                     **dict_plotting)


        # x and y are the data that will be fit to
        x = []
        y = []
        for prof in self.profiles:
            # Get positions
            if prof.positions_microns is None:
                print ''
                print prof.profile_name
                print 'Need to set profile positions'
                return
            x.append(prof.positions_microns)

            # Get whole-block areas
            if peakwn is None and peak_idx is None:
                # bulk hydrogen is the default
                if prof.wb_areas is None:
                    prof.make_3DWB_water_list(polyorder)
                y.append(prof.wb_areas)
            else:
                # peak-specific
                wb_areas, peakwn = prof.get_peak_wb_areas(peak_idx, peakwn)
                y.append(wb_areas)

        # Fitting takes place here        
        if wb_or_3Dnpi == 'wb':
            dict_fitting['raypaths'] = self.raypaths            
            lmfit.minimize(diffusion_3DWB, params, args=(x, y), 
                           kws=dict_fitting)
        elif wb_or_3Dnpi == 'npi':
            diffusion_3D(params, x, y, **dict_fitting)        
            lmfit.minimize(diffusion_3D, params, args=(x, y), kws=dict_fitting)

        # results  
        bestD = []
        bestinit = []
        D3 = []
        bestD.append(ufloat(params['log10Dx'].value, params['log10Dx'].stderr))
        bestD.append(ufloat(params['log10Dy'].value, params['log10Dy'].stderr))
        bestD.append(ufloat(params['log10Dz'].value, params['log10Dz'].stderr))
        bestinit.append(ufloat(params['initial_unit_value_a'].value, 
                             params['initial_unit_value_a'].stderr))
        bestinit.append(ufloat(params['initial_unit_value_b'].value, 
                             params['initial_unit_value_b'].stderr))
        bestinit.append(ufloat(params['initial_unit_value_c'].value, 
                             params['initial_unit_value_c'].stderr))     
        
        for k in range(3):
            D3.append(bestD[k].n)
            if peak_idx is not None:
                if self.profiles[k].peak_diffusivities is None:
                    self.profiles[k].peak_diffusivities = (
                    np.zeros_like(self.profiles[k].peakpos)) 
                    self.profiles[k].peak_diff_error = (
                    np.zeros_like(self.profiles[k].peakpos)) 
                self.profiles[k].peak_diffusivities[peak_idx] = bestD[k].n
                self.profiles[k].peak_diff_error[peak_idx] = bestD[k].s
            else:
                self.profiles[k].diffusivity_log10m2s = bestD[k].n
                self.profiles[k].diff_error = bestD[k].s

        # and then plot after
        if show_plot is True:
            dict_plotting['show_plot'] = True
            dict_plotting['style'] = style_final
            self.plot_diffusion(list_of_log10D_m2s=D3, show_spectra=False,
                                **dict_plotting)

        print '\ntime in hours:', params['time_seconds'].value / 3600.
        print '\ninitial unit values:'
        for initval in bestinit:
            print initval
        print '\nbestfit log10D in m2/s:'
        for D in bestD:
            print D

        return bestD
    
    def invert(self, grid_xyz, symmetry_constraint=True, 
               smoothness_constraint=True, rim_constraint=True, 
               rim_value=None, weighting_factor_lambda=0.2, 
               show_residuals_plot=True):
        """Takes a list of three whole-block concentration profiles (either A/Ao 
        or water ok but must be consistent for all three) in three orthogonal 
        directions and list of three integers to indicate number of divisions
        in each direction. Returns matrix of values in each grid cell. 
        Default plot showing residuals for how well results match the whole-block
        observations."""
        pass
