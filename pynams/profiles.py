"""
Code for processing and plotting *groups* of FTIR spectra.

The central concept is an object class called Profile, which creates groups 
of Spectrum objects, which are defined in a separate module.
"""

from __future__ import print_function, division, absolute_import
from . import styles
from . import diffusion
from . import pynams
from .spectra import Spectrum
import gc
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
from .uncertainties import ufloat
import os.path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import string as string
import pynams.diffusion.lmfit as lmfit
from . import uncertainties
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
#from matplotlib.backends.backend_pdf import PdfPages
#import xlsxwriter
import json
#from scipy import signal as scipysignal
#import scipy.interpolate as interp


class Profile():
    def __init__(self, profile_name=None, time_seconds=None, folder='',
                 fname_list=[], positions_microns = np.array([]),
                 sample=None, direction=None, raypath=None, short_name=None,
                 spectra=[], set_thickness=False,
                 initial_profile=None, base_low_wn=None, base_high_wn=None,
                 diffusivity_log10m2s=None, diff_error=None, length_microns=None,
                 peak_diffusivities=[], peak_diff_error=[], thick_microns=None):
        """fname_list = list of spectra filenames without the .CSV extension.
        Raypath and direction expressed as 'a', 'b', 'c' with thickness/length
        info contained in sample's length_a_microns, length_b_microns, and length_c_microns.
        base_low_wn and base_high_wn can be used to set the wavenumber
        range of the baseline for the spectra.
        
        """
        self.profile_name = profile_name
        self.spectra = spectra
        self.folder = folder
        self.fname_list = fname_list
        self.positions_microns = positions_microns
        self.sample = sample
        self.length_microns = length_microns
        self.direction = direction
        self.raypath = raypath
        self.initial_profile = initial_profile
        self.short_name = short_name
        self.time_seconds = time_seconds
        self.diffusivity_log10m2s = diffusivity_log10m2s
        self.diff_error = diff_error
        self.peak_diffusivities = peak_diffusivities
        self.peak_diff_error = peak_diff_error
        self.thick_microns = thick_microns
        
#        if (self.fname_list is not None) and (self.sample is not None):
        if base_low_wn is not None:
            for spectrum in self.spectra:
                spectrum.base_low_wn = base_low_wn

        if base_high_wn is not None:
            for spectrum in self.spectra:
                spectrum.base_low_wn = base_high_wn

        self.make_spectra(set_thickness=set_thickness)

    short_name = None # short string for saving diffusivities, etc.
    thick_microns_list = None

    
    # for constructing whole-block profiles

    # The actual list of spectra is made automatically by make_spectra.
    # Set spectrum_class_name to use non-default spectra classes
    # e.g., to set different baseline limits consistently
    spectra = []
    spectrum_class_name = None 
    avespec = None # averaged spectra made by self.average_spectra()
    iavespec = None # initial averaged spectra
    length_microns = None # length, but I usually use set_len() directly each time
    waters_list = []
    waters_errors = []
    areas_list = None
    bestfitline_areas = None
    # Bulk whole-block 3D-WB information and diffusivities if applicable
    wb_areas = None 
    wb_water = None
    # for plotting 
    style_base = styles.style_profile
    style_x_marker = None
    style_y_marker = None
    style_z_marker = None
    style_x_line = None
    style_y_line = None
    style_z_line = None
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
    peak_initial_wb_ratio = 1.0    
    # Maximum values to scale up to with diffusion curves
    maximum_area = None
    maximum_wb_area = 1.0
    peak_maximum_areas = None
    peak_maximum_heights = None
    peak_maximum_areas_wb = None
    peak_maximum_heights_wb = None
    
    # Diffusivities (D) and associated errors all in log10 m2/s
    # D's can be derived from any or all of the following
    # - absolute area bulk H: D_area
    # - absolute area peak-specific: D_peakarea
    # - absolute height peak-specific: D_height
    # - whole-block area bulk: D_area_wb 
    # - whole-block area peak-specific: D_peakarea_wb
    # - whole-block height peak-specific: D_height_wb

    # absolute value-derived diffusivities
    D_area = 0.
    D_peakarea = None
    D_height = None

    D_area_error = 0.
    D_peakarea_error = None
    D_height_error = None

    # whole-block-derived diffusivities
    D_area_wb = 0.
    D_peakarea_wb = None
    D_height_wb = None

    D_area_wb_error = 0.
    D_peakarea_wb_error = None
    D_height_wb_error = None
    
    def set_all_thicknesses_from_SiO(self):
        """Individually set thicknesses for all spectra based on the area
        under their Si-O overtone peaks"""
        self.make_spectra()

    def set_len(self):
        """Set profile.length_microns from profile.direction and 
        profile.sample.thick_microns""" 

        if self.sample is None:
            print('\n', self.profile_name)
            print('Need to specify profile sample\n')
            return False
        else:
            s = self.sample

        if self.direction == 'a':
           self.length_microns = s.thickness_microns[0]
        elif self.direction == 'b':
            self.length_microns = s.thickness_microns[1]
        elif self.direction == 'c':
            self.length_microns = s.thickness_microns[2]
        else:
            print('Set direction of profile to a, b, or c')

        return self.length_microns

    def set_thick(self):
        """Set profile.thick_microns from profile.raypath and
        profile.sample.thick_microns"""
        if self.sample is None:
            print('Need to specify profile sample or thickness')
            return False
        else:
            s = self.sample

        if self.raypath == 'a':
           self.thick_microns = s.thickness_microns[0]
        elif self.raypath == 'b':
            self.thick_microns = s.thickness_microns[1]
        elif self.raypath == 'c':
            self.thick_microns = s.thickness_microns[2]
        else:
            print('Need raypath')
            return False
            
        return self.thick_microns

    def plot_thicknesses(self, figaxis=None):
        """Plot thickness across profile"""
        if self.length_microns is None:
            self.length_microns = max(self.positions_microns)+1.
        if figaxis is not None:
            ax = figaxis
        else:
            fig, ax, ax_right = styles.plot_area_profile_outline(self,
                                                            centered=False)
        ax.plot(self.positions_microns, self.thick_microns_list, 'o')
        ax.set_ylabel('thickness ($\mu$m)')
        ax.set_title(self.profile_name)
        ax.set_ylim(min(self.thick_microns_list)-0.05*min(self.thick_microns_list), 
                    max(self.thick_microns_list)+0.05*max(self.thick_microns_list))
        return fig, ax            

    def make_spectra(self, set_thickness=True):
        """Set profile length and generate spectra 
        with key attributes"""
        try:
            if (self.raypath is not None) and (self.raypath == self.direction):
                print("raypath cannot be the same as profile direction")
                return False
        except AttributeError:
            self.raypath = None

        # construct each spectrum from fnames
        if len(self.spectra) == 0:
            if len(self.fname_list) == 0:
                print('Need fnames')
                return False                
            fspectra_list = []
            for x in self.fname_list:
                newspec = Spectrum(fname=x, folder=self.folder)
                newspec.fname = x
                newspec.thick_microns = self.thick_microns
                fspectra_list.append(newspec)
            self.spectra = fspectra_list

        # set sample, raypath for all
        for spec in self.spectra:
            spec.sample = self.sample
            spec.raypath = self.raypath

        if set_thickness is True:
            self.set_thicknesses()
        return

    def set_thicknesses(self):
        """Sets thickness for each spectrum and makes list of thickness for profile"""
#        if self.thick_microns_list is None:
#            self.thick_microns_list = []
        thickness_list = []
        for spec in self.spectra:
            thickness_list.append(spec.thickness_microns)
        self.thick_microns_list = thickness_list
#           if len(self.thick_microns_list) < len(self.spectra):
#            if self.thick_microns is not None:      
#                spec.thick_microns = self.thick_microns
#            elif spec.thick_microns is None:
#                spec.get_thickness_from_SiO()
#            self.thick_microns_list.append(spec.thick_microns)   

    def set_length(self):
        if self.length_microns is None:
            self.length_microns = max(self.thick_microns_list) + 50.

    def average_spectra(self):
        """Creates and returns averaged spectrum and stores it in
        attribute avespec"""
        spec_list = self.spectra
        avespec = Spectrum(folder=None, fname=None)
        avespec.make_average_spectra(spec_list, folder=self.folder)
        
        if self.profile_name is not None:
            avespec.fname = (self.profile_name + '\naverage profile')
        else:
            avespec.fname = 'average profile'
            
        avespec.base_high_wn = spec_list[0].base_high_wn
        avespec.base_low_wn = spec_list[0].base_low_wn
        avespec.start_at_zero()
        self.avespec = avespec
        return avespec

    def plot_spectra(self, show_baseline=True, show_initial_ave=True, 
                     show_final_ave=True, plot_all=False, 
                     initial_and_final_together=False, style=styles.style_spectrum, 
                     stylei=styles.style_initial, wn=None,
                     figsize=(3.2, 3.2)):
        """Plot averaged spectrum across profile. Returns figure and axis."""
        for spec in self.spectra:
            spec.plot_spectrum()
#        if self.spectra is None or len(self.spectra) < 1:
#            self.make_spectra()
#
#        if show_initial_ave is True or initial_and_final_together is True:
#            if self.initial_profile is None:
#                print 'No initial_profile attribute specified'
#                show_initial_ave = False
#                initial_and_final_together = False
#        
#        f = None
#        ax = None
#        
#        if plot_all is True:
#            for spec in self.spectra:
#                if show_baseline is True:
#                    spec.plot_showbaseline(style=style)
#                else:
#                    spec.plot_spectrum(style=style, wn=wn)
#            if show_initial_ave is True:
#                for spec in self.initial_profile.spectra:
#                    if show_baseline is True:
#                        spec.plot_showbaseline(style=stylei)
#                    else:
#                        spec.plot_spectrum(style=stylei, wn=wn)
#                        
#        if show_final_ave is True or initial_and_final_together is True:
#            avespec = self.average_spectra()
#        
#        if show_initial_ave is True or initial_and_final_together is True:
#            initspec = self.initial_profile.average_spectra()
#            
#        # Plot initial average
#        if show_initial_ave is True:
#            if show_baseline is True:
#                initspec.plot_showbaseline(style=stylei)
#            else:
#                initspec.plot_spectrum(style=stylei, wn=wn)
#
#        # Plot final average
#        if show_final_ave is True:
#            if show_baseline is True:
#                f, ax = avespec.plot_showbaseline(style=style)
#            else:            
#                f, ax = avespec.plot_spectrum(style=style, wn=wn)
#        
#        # Plot average spectra together
#        if initial_and_final_together is True:
#            f, ax = avespec.plot_spectrum_outline(self)
#            ax.plot(avespec.wn_full, avespec.abs_full_cm, label='Final', **style)
#            ax.plot(initspec.wn_full, initspec.abs_full_cm, label='Initial', **stylei)            
#            ax.legend()
#            tit = self.profile_name + '\nAverage profiles'
#            ax.set_title(tit)
#            if wn is not None:
#                ax.plot([wn, wn], [ax.get_ylim()[0], ax.get_ylim()[1]], 
#                        color='r')
#        return f, ax

    def change_baseline(self, highwn=3800, lowwn=3000, shift=None):
        """Change baseline parameters for all spectra, final and initial"""
        for spectrum in self.spectra + self.initial_profile.spectra:
            spectrum.base_high_wn = highwn
            spectrum.base_low_wn = lowwn
            if shift is not None:
                spectrum.base_mid_yshift = shift

    def make_composite_peak(self,peak_idx_list):
        """Make composite peaks for all spectra in profile"""
        for spec in self.spectra:
            spec.make_composite_peak(peak_idx_list)
        self.get_peak_info()

    def make_baselines(self, linetype='line', shiftline=None, 
                       wn_high=3700., wn_low=3200., wn_mid=None,
                       show_fit_values=False, show_plot=False,
                       size_inches=(3., 2.5), abs_high=None):
        """Make baselines for all final and initial spectra in profile"""
        if len(self.spectra) < 1:
            self.make_spectra()
        for spectrum in self.spectra:
            spectrum.base_high_wn = wn_high
            spectrum.base_low_wn = wn_low
            spectrum.make_baseline(linetype=linetype, shiftline=shiftline,
                                    show_fit_values=show_fit_values,
                                    show_plot=show_plot, wn_mid=wn_mid,
                                    size_inches=size_inches)

    def get_baselines(self, initial_too=False, folder=None, delim=',', 
                      baseline_ending='-baseline.CSV'):
        """Get previously saved baselines for all spectra in profile"""
        for spectrum in self.spectra:
            spectrum.get_baseline(baseline_ending=baseline_ending,
                                  folder=folder, delim=delim)
            
        if initial_too is True:
            for spectrum in self.initial_profile.spectra:
                spectrum.get_baseline()

    def matlab(self):
        """Print out spectra fnames for FTIR_peakfit_loop.m"""
        string = "{"
        for spec in self.spectra:
            stringname = spec.fname
            string = string + "'" + stringname + "' "
        string = string + "};"
        print('\nfilenames ready for FTIR_peakfit_loop.m:')
        print(string)
                        
    def save_baselines(self, printnames=True):
        """Save all baselines in profile"""
        for spectrum in self.spectra:
            spectrum.save_baseline()
        
        if printnames is True:
            self.print_names4matlab()
            
    def make_area_list(self, polyorder=1, show_plot=False, set_class=None,
                       shiftline=None, printout_area=False, peak=None):
        """Make list of areas (no error) under the curve for an FTIR profile.
        Default is bulk area. Set peak=wavenumber for peak-specific profile"""
        # You need the list of spectra for the profile
        if len(self.spectra) < 1:
            check = self.make_spectra(class_from_module=set_class)
            if check is False:
                return False
              
        areas = []
        if peak is None:
            if self.areas_list is not None:
                print('Overwriting previous areas list')
                            
            if self.spectra[0].area is None:                
                print('generating bulk areas under curves...')
                for spec in self.spectra:
                    a = spec.area_under_curve(polyorder, show_plot, shiftline, 
                                           printout_area, 
                                           require_saved_baseline=False)
                    areas.append(a)            
                self.areas_list = np.array(areas)
            else:
                for spec in self.spectra:
                    areas.append(spec.area)
                
            self.areas_list = areas
            
        else:
            peaklist = list(self.spectra[0].peakpos)
            print('peak at', peak)

            if peak in peaklist:
                idx = peaklist.index(peak)
                for x in self.spectra:
                    a = x.peak_areas[idx]
                    areas.append(a)
            else:
                print('No peak at wavenumber', peak)
                print('peak positions:', peaklist)
                return False
        return areas
        
    def get_peakfit(self, peak_ending='-peakfit.CSV',
                    baseline_ending='-baseline.CSV'):
        """Get fit peaks from MATLAB for all spectra in profile, including 
        in the initial profile. The resulting numpy arrays are in dimensions
        (number of peaks, number of spectra in profile) and stored in
        the profiles's attributes peak_heights, peak_widths, and peak_areas"""
        for spectrum in self.spectra:
            spectrum.get_peakfit(peak_ending=peak_ending,
                                 baseline_ending=baseline_ending)
        self.get_peak_info()
            
    def get_peak_info(self):
        """Pull peak info from individual spectra into a single profile
        attribute"""
        if len(self.spectra) < 1:
            self.make_spectra()
        peakpos = self.spectra[0].peakpos

        if peakpos is None:
            return False
        hbig = []
        wbig = []
        abig = []
        
        for p in range(len(peakpos)):
            h = []
            w = []
            a = []        
            for spectrum in self.spectra:
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
        
        self.peak_wb_areas = np.zeros_like(peakpos)
        self.D_peakarea = np.zeros_like(peakpos)
        self.D_peakarea_error = np.zeros_like(peakpos)
        self.D_peakarea_wb = np.zeros_like(peakpos)
        self.D_peakarea_wb_error = np.zeros_like(peakpos)
        self.peak_maximum_areas = np.zeros_like(peakpos)
        self.peak_maximum_areas_wb = np.ones_like(peakpos)        
        
        self.D_height = np.zeros_like(peakpos)
        self.D_height_error = np.zeros_like(peakpos)
        self.D_height_wb = np.zeros_like(peakpos)
        self.D_height_wb_error = np.zeros_like(peakpos)
        self.peak_maximum_heights = np.zeros_like(peakpos)
        self.peak_maximum_heights_wb = np.ones_like(peakpos)

    def get_area_total(self):
        """Sum up total area of all peaks"""
        total_peak_area = []
        total_height = []
        
        for spec in self.spectra:
            spec.get_peakareas()
            total_peak_area.append(sum(spec.peak_areas))
            total_height.append(sum(spec.peak_heights))
        total_area = np.mean(total_peak_area)
        total_h = np.mean(total_height)
        self.total_peak_area = total_area
        self.total_peak_height = total_h
        return total_area, total_h
            
    def print_peakfits(self):
        """Print out peakfit information for all spectra in profile"""
        print('\n', self.profile_name)

        poscounter = 0
        for spectrum in self.spectra:
            print('\n', spectrum.fname, \
                    self.positions_microns[poscounter], 'microns')
            poscounter += 1
            
            if spectrum.peakpos is None:
                check = spectrum.get_peakfit()
                if check is False:
                    print('trouble with getting peakfit')
            
            if spectrum.peakpos is not None:
                for k in range(len(spectrum.peakpos)):
                    print(spectrum.peakpos[k], spectrum.peak_heights[k], \
                          spectrum.peak_widths[k], spectrum.peak_areas[k])        


    def print_peakfits_ave(self, printout=True):
        """Computes, prints, and returns average values in profile for
        peak areas, peak heights, and the summed average peak areas"""
        # get peak information for all spectra
        for spectrum in self.spectra:
            if spectrum.peak_areas is None:
                self.get_peakfit()

        average_peakareas = np.average(self.peak_areas, axis=1)
        average_peakheights = np.average(self.peak_heights, axis=1)
        total_area = np.sum(average_peakareas)
        
        max_peakareas = np.max(self.peak_areas, axis=1)
        max_peakheights = np.max(self.peak_heights, axis=1)

        if printout is True:
            print('\n', self.profile_name)
            print('peak positions (cm-1)\n', self.spectra[0].peakpos)
            print('average peak areas (cm-2)\n', average_peakareas)
            print('average peak heights (cm-1)\n', average_peakheights)
            print('summed areas (cm-2)', total_area)
            print('max peak areas (cm-2)\n', max_peakareas)
            print('max peak heights (cm-1)\n', max_peakheights)
            
        return average_peakareas, average_peakheights, total_area, \
                max_peakareas, max_peakheights

    def plot_peakfits(self, initial_too=False, legloc=1, top=1.2):
        """Show fit peaks for all spectra in profile"""
        if len(self.spectra) < 1:
            self.make_spectra()
            
        for spectrum in self.spectra:
            spectrum.plot_peakfit(legloc=legloc, top=top)
        
        if initial_too is True and self.initial_profile is not None:
            for spectrum in self.initial_profile.spectra:
                spectrum.plot_peakfit(legloc=legloc)
                
    def make_wholeblock(self, peakfit=True, show_plot=False, bulk=True):
        """Take initial and final profiles and make self.wb_areas.
        If peakfit=True, then get_peakfit and make peak-specific ratios
        peak_wb_areas, peak_wb_heights, peak_wb_widths"""

        if self.initial_profile is None:
            print('\n profile name:', self.profile_name)
            print('Need to specify an initial profile')
            return False

        # Bulk H whole-block
        if bulk is True:
            # get both initial and final area lists
            init = self.initial_profile
            fin = self
            for prof in [init, fin]:
                if prof.areas_list is None:
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
                print('\nCould not get peak areas')
                return False
            
        if self.initial_profile.peak_areas is None:
            self.initial_profile.get_peakfit()
            if self.initial_profile.peak_areas is None:
                print('\nCould not get initial peak areas')
                return False
            
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
            
            wb_areas.append(areas[k] / anormalizeto)                
            wb_heights.append(heights[k] / hnormalizeto)
            wb_widths.append(widths[k] / wnormalizeto)

#            print 'peak #', k
#            print 'normalizing to areas', anormalizeto
#            print 'areas', areas[k]
#            print 'whole-block areas', areas[k] / anormalizeto

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
            print('Getting peaks fit in matlab for', self.profile_name)
            self.get_peakfit()

        # Getting peak from wavenumber if index not given
        # Or peak wavenumber from index
        if peak_idx is None:
            if peakwn not in self.peakpos:
                print('There is no peak at', peakwn)
                print('Peaks are at', self.peakpos)
                return
            peak_idx = np.where(self.peakpos==peakwn)[0][0]
            print('peak at', peakwn, 'is index', peak_idx)
        else:
            peakwn = self.peakpos[peak_idx]

        # Make peak-specific whole-block areas
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
            print('Set base style (style_base) for profile first')
            return False
        self.style_x_marker =  dict(list(self.style_base.items()) + list(styles.style_Dx.items()))
        self.style_y_marker =  dict(list(self.style_base.items()) + list(styles.style_Dy.items()))
        self.style_z_marker =  dict(list(self.style_base.items()) + list(styles.style_Dz.items()))
        self.style_x_line = styles.make_line_style('x', self.style_base)
        self.style_y_line = styles.make_line_style('y', self.style_base)
        self.style_z_line = styles.make_line_style('z', self.style_base)
        
        if self.raypath == 'a':
            self.style_x_marker.update({'marker' : styles.style_Rx['marker']})
            self.style_y_marker.update({'marker' : styles.style_Rx['marker']})
            self.style_z_marker.update({'marker' : styles.style_Rx['marker']})
        elif self.raypath == 'b':
            self.style_x_marker.update({'marker' : styles.style_Ry['marker']})
            self.style_y_marker.update({'marker' : styles.style_Ry['marker']})
            self.style_z_marker.update({'marker' : styles.style_Ry['marker']})
        elif self.raypath == 'c':
            self.style_x_marker.update({'marker' : styles.style_Rz['marker']})        
            self.style_y_marker.update({'marker' : styles.style_Rz['marker']})
            self.style_z_marker.update({'marker' : styles.style_Rz['marker']})
        return

    def choose_line_style(self):
        """Returns line style with direction information"""
        if self.style_base is None:
            print('Using default styles. Set profile style_base to change')
            self.style_base = styles.style_profile            
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
            print('Using default styles. Set profile style_base to change')
            self.style_base = styles.style_profile

        if self.style_x_marker is None:
            self.make_style_subtypes()

        if self.direction == 'a':
            style = self.style_x_marker
        elif self.direction == 'b':
            style = self.style_y_marker
        elif self.direction == 'c':
            style = self.style_z_marker
        else:
            style = self.style_base

        return style
        
    def plot_area_profile(self, polyorder=1, centered=True, 
                          heights_instead=False, figaxis=None, 
                          bestfitline=False, style_bestfitline=None,
                          show_FTIR=False, show_water_ppm=True,
                          show_values=False, set_class=None,
                          peakwn=None, peak_idx=None,
                          style=styles.style_points, show_initial_areas=False,
                          error_percent=0, wholeblock=False,
                          label=None, initial_style=None,
                          initial_label=None, phase='olivine',
                          calibration='Bell', orientation_factor=3.):
        """Plot area profile. Centered=True puts 0 in the middle of the x-axis.
        figaxis sets whether to create a new figure or plot on an existing axis.
        bestfitline=True draws a best-fit line through the data.
        Set peak=wavenumber for peak-specific profile"""
        if wholeblock is True and self.initial_profile is None:
            print('Need to specify an initial profile')
            return False, False
        
        if len(self.spectra) < 1:
            self.make_spectra()
        
        # Check for or create positions and areas
        if len(self.positions_microns) < 1:
            print('Need positions_microns for profile')
            return

        # Get list of areas
        if wholeblock is False:
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
                    print('peak at', peakwn, 'is index', peak_idx)
                else:
                    peakwn = self.peakpos[peak_idx]

            if heights_instead is True and peak_idx is not None:
                areas = self.peak_heights[peak_idx]
            elif peak_idx is not None:
                areas = self.peak_areas[peak_idx]

        # whole-block
        else:
            # bulk hydrogen
            if peak_idx is None:
                if self.wb_areas is None:
                    self.make_wholeblock(peakfit=False, bulk=True)
                areas = self.wb_areas
            # peak-specific
            else:
                if self.peak_wb_areas is None:
                    self.make_wholeblock(peakfit=True, bulk=False)
    
                if heights_instead is False:
                    areas = self.peak_wb_areas[peak_idx]
                else:
                    areas = self.peak_wb_heights[peak_idx]
                    
        if areas is False:
            return            

        if np.shape(areas) != np.shape(self.positions_microns):
            print('Area and positions lists are not the same size!')
            print('area:', np.shape(self.areas_list))
            print('positions:', np.shape(self.positions_microns))
            return

        # Use new or old figure axes
        if figaxis is None:
            f, ax, ax_ppm = styles.plot_area_profile_outline(self, centered, 
                            peakwn, heights_instead=heights_instead, 
                            wholeblock=wholeblock,
                            show_water_ppm=show_water_ppm)
        else:
            ax = figaxis
            show_water_ppm = False

        # Set length
        if self.length_microns is None:
            leng = self.set_len()
            if self.length_microns is None:
                print('Need some information about profile length')
                return
        else:
            leng = self.length_microns

        # Set up plotting styles
        if style_bestfitline is None:
            style_bestfitline = self.choose_line_style()
        if style is None:
            style = self.choose_marker_style()
        if label is None:
            style['label'] = self.profile_name
        else:
            style['label'] = label

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
            x = np.array(self.positions_microns) - leng/2.0
        else:
            x = self.positions_microns            
        
        yerror = np.array(areas)*error_percent/100.
        
        if error_percent == 0:
            ax.plot(x, areas, **style)
        else:
            ax.errorbar(x, areas, yerr=yerror, **style)
            
        # Plot initial profile areas
        if show_initial_areas is True:
            if initial_style is None:
                initial_style = styles.style_initial_point
            if initial_label is None:
                initial_style['label'] = 'initial'
            else:
                initial_style['label'] = label
            self.initial_profile.plot_area_profile(style=initial_style,
                                                   figaxis=ax)            
        ax.set_ylim(0, max(areas)+0.2*(max(areas)))

        # Title
        if peak_idx is None:
            tit = 'Bulk hydrogen'
        else:
            peakwn = self.peakpos[peak_idx]
            tit = ' '.join(('Peak at', str(peakwn) ,'/cm'))
        ax.set_title(tit)

        if show_water_ppm is True:
            ax_ppm.set_ylabel(''.join(('ppm H2O in ', phase, ', ', calibration, 
                                       ' calibration *', str(orientation_factor))))
            parasite_tick_locations = np.linspace(ax.get_ylim()[0],
                                                  ax.get_ylim()[1], 5)
            abs_coeff = pynams.absorption_coefficients(phase=phase, 
                                                       calibration=calibration
                                                       )
#            ppm_labels = parasite_tick_locations * abs_coeff * orientation_factor
            ppm_labels = parasite_tick_locations * abs_coeff.n * orientation_factor

            ax_ppm.set_yticks(parasite_tick_locations)
#            ax_ppm.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
            np.set_printoptions(precision=1)
            ax_ppm.set_yticklabels(['{:.1f}'.format(i) for i in ppm_labels])
#    parasite_tick_locations = 1e4/(celsius_labels + 273.15)
#    ax_celsius.set_xticks(parasite_tick_locations)
#    ax_celsius.set_xticklabels(celsius_labels)
#    fig.add_subplot(ax)
#    ax.axis["bottom"].set_label("10$^4$/Temperature (K$^{-1}$)")
#    ax.axis["left"].set_label("log$_{10}$diffusivity (m$^{2}$/s)")
#    ax_celsius.axis["top"].set_label("Temperature ($\degree$C)")
#    ax_celsius.axis["top"].label.set_visible(True)
#    ax_celsius.axis["right"].major_ticklabels.set_visible(False)

#            max_water = area2water(ax.get_ylim()[1]*3, phase=phase, 
#                                   calibration=calibration)
#            ax_ppm.axis["right"].major_ticklabels.set_visible(False) 
#            ax_ppm.axis["right"].set_label(''.join(("max ppm H$_2$O ~ ", 
#                                                    str(max_water))))
         
        if figaxis is None:
            return f, ax
        else:
            return


#    def plot_wb_water(self, polyorder=1, centered=True, style=styles.style_profile):
#        """Plot area ratios and scale up with initial water"""
#        if self.sample.initial_water is None:
#            print 'Set sample initial water content please.'
#            return
#        a, w = self.make_3DWB_water_list(polyorder=1)
#        
#        fig, ax, leng = styles.plot_area_profile_outline(self, centered)
#        top = max(a) + 0.2*max(a)
#        ax.set_ylim(0, top)
#        ax.set_ylabel('Final area / Initial area')
#        
#        if centered is True:
#            ax.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
#            ax.plot(self.positions_microns - leng/2.0, a, **style)           
#        else:
#            ax.plot([0, leng], [1, 1], **styles.style_1)
#            ax.plot(self.positions_microns, a, **style)
#        return fig, ax, leng

    def y_data_picker(self, wholeblock, heights_instead, peak_idx=None):
        """Pick out and return peak area or height data of interest."""
        if wholeblock is True:
            if peak_idx is None:
                if self.wb_areas is None:
                    self.make_wholeblock(peakfit=False, bulk=True)
                y = self.wb_areas
            # Peak-specific 
            else:
                if self.peak_wb_areas is None:
                    self.make_wholeblock(peakfit=True, bulk=False)
                    
                if heights_instead is False:
                    y = self.peak_wb_areas[peak_idx]
                else:
                    y = self.peak_wb_heights[peak_idx]
        else:
            if self.areas_list is None:
                self.make_area_list()
            y = self.areas_list
        return y

    def D_picker(self, wholeblock=True, heights_instead=False, peak_idx=None):
        """Returns attribute name and diffusivity of interest. 
        Consider using get_diffusivity() first"""
        log10D_m2s = None
        
        if peak_idx is not None:
            if heights_instead is True:
                if wholeblock is False and self.D_height[peak_idx] != 0.0:
                    log10D_m2s = self.D_height[peak_idx]
                elif wholeblock is True and self.D_height_wb[peak_idx] != 0.0:
                    log10D_m2s = self.D_height_wb[peak_idx]
            else:
                if wholeblock is False and self.D_peakarea[peak_idx] != 0.0:
                    log10D_m2s = self.D_peakarea[peak_idx]
                elif wholeblock is True and self.D_peakarea_wb[peak_idx] != 0.0:
                    log10D_m2s = self.D_peakarea_wb[peak_idx]
        elif wholeblock is False and self.D_area != 0.0:
            log10D_m2s = self.D_area
        elif self.D_area_wb != 0.0:
            log10D_m2s = self.D_area_wb
            
        if log10D_m2s is None or log10D_m2s == 0.:
            print('\nNeed to set profile bulk or peak diffusivity')
            return None
            
        return log10D_m2s

    def D_saver(self, D, error, wholeblock=True, heights_instead=False, 
                peak_idx=None):
        """Saves a diffusivity in the profile attribute of interest.
        Consider using get_diffusivity() first"""
        if peak_idx is not None:
            if heights_instead is True:
                if wholeblock is False:
                    self.D_height[peak_idx] = D
                    self.D_height_error[peak_idx] = error
                elif wholeblock is True:
                    self.D_height_wb[peak_idx] = D
                    self.D_height_wb_error[peak_idx] = error
            else:
                if wholeblock is False:
                    self.D_peakarea[peak_idx] = D
                    self.D_peakarea_error[peak_idx] = error
                elif wholeblock is True:
                    self.D_peakarea_wb[peak_idx] = D
                    self.D_peakarea_wb_error[peak_idx] = error
        elif wholeblock is False:
            self.D_area = D
            self.D_area_error = error
        else:
            self.D_area_wb = D
            self.D_area_wb_error = error

    def scaling_factor_picker(self, maximum_value=None, wholeblock=True,
                              heights_instead=False, peak_idx=None):
        """Returns value to scale diffusion to"""
        if maximum_value is not None:
            scaling_factor = maximum_value
            
        elif wholeblock is True:
            if peak_idx is None:
                scaling_factor = self.maximum_wb_area
            else:
                if heights_instead is False:
                    scaling_factor = self.peak_maximum_areas_wb[peak_idx]
                else:
                    scaling_factor = self.peak_maximum_heights_wb[peak_idx]
        
        else: 
            # bulk water                
            if peak_idx is None:
                if len(self.areas_list) < 1:
                    self.make_area_list()
                if self.maximum_area is not None:
                    scaling_factor = self.maximum_area
                else:
                    scaling_factor = np.max(self.areas_list)
                    self.maximum_area = scaling_factor
        
            # peak-specific heights                
            elif heights_instead is True: 
                if self.peak_maximum_heights[peak_idx] != 0.0:
                    scaling_factor = self.peak_maximum_heights[peak_idx]
                else:                
                    scaling_factor = np.max(self.peak_heights[peak_idx])
                    self.peak_maximum_heights[peak_idx] = scaling_factor
            # peak-specific areas
            else:
                if self.peak_maximum_areas[peak_idx] != 0.0:
                    scaling_factor = self.peak_maximum_areas[peak_idx]
                else:
                    scaling_factor = np.max(self.peak_areas[peak_idx])
                    self.peak_maximum_areas[peak_idx] = scaling_factor
    
        return scaling_factor

    def plot_diffusion(self, log10D_m2s=None, time_seconds=None, 
                       peak_idx=None, top=1.2, wholeblock=False,
                       maximum_value=None, style=None, centered=False,
                       heights_instead=False, points=200., symmetric=True,
                       labelD=True, labelDx=None, labelDy=None,
                       erf_or_sum='erf', label4legend=None,
                       initial_unit_value=1., final_unit_value=0.):
        """Plot diffusion curve with profile data."""
        if wholeblock is True and self.initial_profile is None:
            print('Need to specify an initial profile')
            return False, False
        
        if symmetric is False:
            centered = False
            
        fig, ax = self.plot_area_profile(wholeblock=wholeblock, 
                                         peak_idx=peak_idx, centered=centered,
                                         heights_instead=heights_instead)
                                                       
        if self.length_microns is not None:
            microns = self.length_microns
        elif self.sample.thickness_microns is not None:
            microns = self.set_len()
        else:
            print('\nNeed to set profile length directly or with a sample')
            return False, False

        # Get time
        if time_seconds is not None:
           pass
        elif self.time_seconds is not None:
            time_seconds = self.time_seconds
        else:
            print('\nNeed to set profile time_seconds attribute')
            return False, False, False

        # Get diffusivity - input directly or stored in an attribute
        if log10D_m2s is None:
           log10D_m2s = self.D_picker(wholeblock, heights_instead, 
                                      peak_idx)   
        if log10D_m2s is None:
            print('\nNeed to set profile bulk or peak diffusivity')
            return False, False, False

        if peak_idx is not None and self.peakpos is None:
            self.get_peakfit()

        # Get scaling factor for diffusion curve
        scaling_factor = self.scaling_factor_picker(maximum_value, 
                                        wholeblock, heights_instead, peak_idx)
                        
        ax.plot(ax.get_xlim(), [scaling_factor, scaling_factor], '--k')
        print('\nScaling to {:.2f} maximum value'.format(scaling_factor))
#        print 'You can set maximum_value to scale diffusion curve'
       
        # Setup and plot diffusion curves
        if symmetric is True:
            params = diffusion.params_setup1D(microns, log10D_m2s, time_seconds,
                                              init=initial_unit_value, 
                                              fin=final_unit_value)
            x_diffusion, y_diffusion = diffusion.diffusion1D_params(params, 
                                                                    points=points)        
            if centered is False:
                x_diffusion = x_diffusion + (self.length_microns/2.)
        else:
            params = diffusion.params_setup1D(microns*2, log10D_m2s, 
                                              time_seconds,
                                              init=initial_unit_value, 
                                              fin=final_unit_value)
            x_diffusion, y_diffusion = diffusion.diffusion1D_params(params, 
                                                                points=points)        
            x_diffusion = x_diffusion[int(points/2):]
            y_diffusion = y_diffusion[int(points/2):]

        ### FINALLY PLOTTING
        ax.plot(x_diffusion, y_diffusion*scaling_factor, label=label4legend)
            
            
        # label diffusivity on plot
        if labelD is True:
            if labelDx is None:
                if centered is True:
                    labelDx = -microns/2.
                else:
                    labelDx = self.length_microns * 0.05
            if labelDy is None:
                top = ax.get_ylim()[1]
                labelDy = top-0.15*top
            strD = "{:.2f}".format(log10D_m2s) 
            
            ax.text(labelDx, labelDy,
                    ''.join(('  D=10$^{', strD, '}$ m$^2$/s')), 
                    ha='left', fontsize=16, backgroundcolor='w') 
                  
        return fig, ax

    def diffusion_residuals(self, time_seconds, log10D_m2s, wholeblock=True,
                            heights_instead=False, peak_idx=None,
                            initial_unit_value=1., final_unit_value=0.,
                            show_plot=True, top=1.2, 
                            maximum_value=None):
        """Compare 1D diffusion curve with profile data.
        Returns vector containing the residuals and 
        and the variance = sqrt(sum of squares of the residuals)"""
        if self.positions_microns is None:
            print('Need to set profile positions')
            return

        if peak_idx is not None:
            if peak_idx is not None and self.peakpos is None:
                self.get_peakfit()

        y = self.y_data_picker(wholeblock, heights_instead, peak_idx)
        scaling_factor = self.scaling_factor_picker(maximum_value, 
                                        wholeblock, heights_instead, peak_idx)
        
        x = self.positions_microns
        L = self.set_len()
        t = time_seconds

        # Need to work on this        
        init = scaling_factor
            
        fin = final_unit_value
        
        params = diffusion.params_setup1D(L, log10D_m2s, t, init, fin, 
                                          False, False, False)

        xdif, model = diffusion.diffusion1D_params(params)
        resid = diffusion.diffusion1D_params(params, x, y)
        RSS = np.sum(resid**2)
#        plt.plot(x-L/2., y, '+k')
#        plt.plot(xdif, model, '-r')
        
        if show_plot is True:
            f, ax = self.plot_diffusion(log10D_m2s, t, peak_idx, top,
                                        wholeblock=wholeblock, 
                                        heights_instead=heights_instead,
                                        maximum_value=maximum_value)
        return resid, RSS
            
    def fitD(self, time_seconds=None, points=200, 
             initial_unit_value=1., vary_initial=True,
             varyD=True, guess=-13., peak_idx=None, top=1.2, 
             peakwn=None, wholeblock=True, centered=False,
             show_plot=True, polyorder=1, heights_instead=False,
             final_unit_value=0., vary_final=False, 
             maximum_value=None, min_water=0., symmetric=False):
        """Fits 1D diffusion curve to profile data."""
        if self.time_seconds is None and time_seconds is None:
            print('Need time_seconds')
            return
        elif time_seconds is None:
            time_seconds = self.time_seconds            

        if self.length_microns is None:
            print('Need to set profile attribute length_microns')
            return
            
        if self.positions_microns is None:
            print('Need to set profile positions')
            return

        if self.areas_list is None:
            self.make_area_list()
        
        ### Choose y data to fit to ###
        y = self.y_data_picker(wholeblock, heights_instead, peak_idx)
        scaling_factor = max(y)
        y = y / scaling_factor # scale down to unit, so y range between 0 and 1
        
        # Set up x data and other parameters
        x = self.positions_microns
        L = self.length_microns
        t = time_seconds
        D0 = guess
        init = initial_unit_value
        fin = final_unit_value
        params = diffusion.params_setup1D(L, D0, t, init, fin, vD=varyD, 
                                          vinit=vary_initial, vfin=vary_final)

        dict_fitting = {'points' : points, 
#                        'symmetric' : symmetric,
#                        'centered' : centered
                        }

#        # minimization
        lmfit.minimize(diffusion.diffusion1D_params, params, args=(x, y), 
                       kws=dict_fitting)
        best_D = ufloat(params['log10D_m2s'].value, 
                        params['log10D_m2s'].stderr)
        best_init = ufloat(params['initial_unit_value'].value, 
                         params['initial_unit_value'].stderr)   
        print('best-fit log10D m2/s', best_D)
        print('best-fit initial    ', best_init*scaling_factor)

#        # save results as attributes
#        # Use save_diffusivities to save to a file
#        if wholeblock is True:
#            if peak_idx is not None:
#                if heights_instead is False:
#                    self.D_peakarea_wb[peak_idx] = best_D.n
#                    self.D_peakarea_wb_error[peak_idx] = best_D.s
#                    self.peak_maximum_areas_wb[peak_idx] = best_init.n
#                else:
#                    self.D_height_wb[peak_idx] = best_D.n
#                    self.D_height_wb_error[peak_idx] = best_D.s
#                    self.peak_maximum_heights_wb[peak_idx] = best_init.n
#
#            else:
#                self.D_area_wb = best_D.n
#                self.D_area_wb_error = best_D.s
#                self.maximum_wb_area = best_init.n
#        
#        resid, RSS = self.diffusion_residuals(best_D.n, wholeblock, 
#                          heights_instead, peak_idx, best_init.n, 
#                          final_unit_value, show_plot=False, 
#                          maximum_value=best_init.n)

        resid, RSS = self.diffusion_residuals(time_seconds=time_seconds, 
                                              log10D_m2s=best_D.n,
                                              wholeblock=wholeblock,
                                              heights_instead=heights_instead,
                                              peak_idx=peak_idx,
                                              initial_unit_value=best_init.n,
                                              final_unit_value=final_unit_value,
                                              show_plot=False,
                                              maximum_value=best_init.n)
        # report results
        print('\ntime in hours:', params['time_seconds'].value / 3600.)
        print('initial unit value:', '{:.2f}'.format(best_init*scaling_factor))
        print('bestfit log10D in m2/s:', '{:.2f}'.format(best_D))
        print('residual sum of squares', '{:.2f}'.format(RSS))
        if show_plot is True:
            fig, ax = self.plot_diffusion(time_seconds=time_seconds,
                                          log10D_m2s=best_D.n,
                                          peak_idx=peak_idx, top=top, 
                                          wholeblock=wholeblock,
                                          centered=centered,
                                          symmetric=symmetric,
                                          heights_instead=heights_instead, 
                                          maximum_value=best_init.n*scaling_factor,
                                          final_unit_value=final_unit_value
                                          )
            return fig, ax
#        else:
#            return 1, 2
##        return best_D, best_init, RSS

    def save_diffusivities(self, folder=None, 
                           file_ending='-diffusivities.txt'):
        """Save diffusivities for profile to a file"""
        if folder is None:
            folder = self.folder
            
        if self.short_name is None:
            print('Need profile short_name attribute to label file')
            return
            
        if self.peakpos is None:
            self.get_peakfit()

        # Same format as output by print_diffusivities:
        # Diffuvities from WB areas - area - WB heights - heights 
        # Diffusivity - error - maximum value to scale up to for each

        # Saving diffusivities in first column, errors in 2nd, 
        # Bulk H in first row, then each peak-specific value
        a = []
        a.append([self.D_area_wb, self.D_area_wb_error, self.maximum_wb_area])
        a.append([self.D_area,    self.D_area_error,    self.maximum_area])
        
        for k in range(len(self.peakpos)):
            a.append([self.D_peakarea_wb[k], self.D_peakarea_wb_error[k],
                     self.peak_maximum_areas_wb[k]])
            a.append([self.D_peakarea[k], 
                      self.D_peakarea_error[k], 
                      self.peak_maximum_areas[k]])
            a.append([self.D_height_wb[k], self.D_height_wb_error[k],
                     self.peak_maximum_heights_wb[k]])
            a.append([self.D_height[k],
                     self.D_height_error[k], self.peak_maximum_heights[k]])
        
        workfile = ''.join((folder, self.short_name, file_ending))
        with open(workfile, 'w') as diff_file:
            diff_file.write(json.dumps(a))


    def get_diffusivities(self, folder=None, 
                           file_ending='-diffusivities.txt'):
        """Get saved diffusivities for profile from a file"""
        if folder is None:
            folder = self.folder
        
        # Look into json format        
        if self.short_name is None:
            print('Need profile short_name attribute to label file')
            return
            
        if self.peakpos is None:
            self.get_peakfit()

        workfile = ''.join((folder, self.short_name, file_ending))
        if os.path.isfile(workfile) is False:
            print(' ')
            print(self.fname)      
            print('use save_diffusivities() to make -diffusivities.txt')
            return

        with open(workfile, 'r') as diff_file:
            diffusivities_string = diff_file.read()

        diffusivities = json.loads(diffusivities_string)

        self.D_area_wb = diffusivities[0][0]
        self.D_area_wb_error = diffusivities[0][1]
        self.maximum_wb_area = diffusivities[0][2]
        
        self.D_area = diffusivities[1][0]
        self.D_area_error = diffusivities[1][1]
        self.maximum_area = diffusivities[1][2]
        
        if self.peakpos is None:
            self.get_peakfit()
        if self.peakpos is None:
            print('Having trouble getting peakfit info to grab diffusivities')
            return            
                    
        npeaks = len(self.peakpos)
        for k in range(npeaks):
            self.D_peakarea_wb[k] = diffusivities[2+4*k][0]
            self.D_peakarea_wb_error[k] = diffusivities[2+4*k][1]
            self.peak_maximum_areas_wb[k] = diffusivities[2+4*k][2]
            self.D_peakarea[k] = diffusivities[3+4*k][0]
            self.D_peakarea_error[k] = diffusivities[3+4*k][1]
            self.peak_maximum_areas[k] = diffusivities[3+4*k][2]
            self.D_height_wb[k] = diffusivities[4+4*k][0]
            self.D_height_wb_error[k] = diffusivities[4+4*k][1]
            self.peak_maximum_heights_wb[k] = diffusivities[4+4*k][2]
            self.D_height[k] = diffusivities[5+4*k][0]
            self.D_height_error[k] = diffusivities[5+4*k][1]
            self.peak_maximum_heights[k] = diffusivities[5+4*k][2]
        
        return diffusivities
        
    def print_diffusivities(self, show_on_screen=True):
        """Print out all diffusivities, including peak-specific, in profile"""
        if self.peakpos is None:
            self.get_peakfit()

        D_area_wb = ufloat(self.D_area_wb, self.D_area_wb_error)   
        wbmax = '{:.2f}'.format(self.maximum_wb_area)

        print('\n', self.profile_name)
        print(' Diffusivities as log10(D in m2/s), errors, then max value A/A0')
        print('bulkH', D_area_wb, wbmax)

        # peak-specific
        for k in range(len(self.peakpos)):
            D_area_wb = ufloat(self.D_peakarea_wb[k], 
                               self.D_peakarea_wb_error[k])   
            awb = '{:.2f}'.format(self.peak_maximum_areas_wb[k])            
            Da_wb = '{:.2f}'.format(D_area_wb)            
            str1 = ' '.join((Da_wb, str(awb)))
            string = ' '.join((str(self.peakpos[k]), str1))
            print(string)

#        ### Showing all ways to generate ###
#        # bulk H diffusivities and initial values
#        D_area_wb = ufloat(self.D_area_wb, self.D_area_wb_error)   
#        D_area = ufloat(self.D_area, self.D_area_error)        
#        wb = '{:.2f}'.format(D_area_wb)
#        a = '{:.2f}'.format(D_area)
#        wbmax = '{:.2f}'.format(self.maximum_wb_area)
#        if self.maximum_area is not None:
#            amax = '{:.1f}'.format(self.maximum_area)  
#        else:
#            amax = 'n/a'
#        st1 = ''.join((wb, '(', wbmax, ')'))
#        st2 = ''.join((a, '(', amax, ')'))
#        bulkstring = ''.join(('bulk H : ', st1, ';  ', st2, 
#                              ';  n/a;         n/a'))
#
#        if show_on_screen is True:
#            print '\n', self.profile_name
#            print ' Diffusivities as log10(D in m2/s) (max value)'
#            print '         wb areas;         areas;       wb heights;     heights'
#            print bulkstring
#
#        # peak-specific
#        peakstrings = []
#        for k in range(len(self.peakpos)):
#            D_area_wb = ufloat(self.D_peakarea_wb[k], 
#                               self.D_peakarea_wb_error[k])   
#            D_area = ufloat(self.D_peakarea[k], self.D_peakarea_error[k])   
#            D_h_wb = ufloat(self.D_height_wb[k], self.D_height_wb_error[k])   
#            D_h = ufloat(self.D_height[k], self.D_height_error[k])
#            a = '{:.2f}'.format(self.peak_maximum_areas[k])
#            awb = '{:.2f}'.format(self.peak_maximum_areas_wb[k])
#            h = '{:.2f}'.format(self.peak_maximum_heights[k])
#            hwb = '{:.2f}'.format(self.peak_maximum_heights_wb[k])
#            
#            Da = '{:.2f}'.format(D_area)
#            Da_wb = '{:.2f}'.format(D_area_wb)
#            Dh = '{:.2f}'.format(D_h)
#            Dh_wb = '{:.2f}'.format(D_h_wb)   
#            
#            st1 = ''.join((Da_wb, '(', str(awb), ')'))
#            st2 = ''.join((Da, '(', str(a), ')'))
#            st3 = ''.join((Dh_wb, '(', str(hwb), ')'))
#            st4 = ''.join((Dh, '(', str(h), ')'))
#            
#            string0 = str(self.peakpos[k])
#            stringD = ';  '.join((st1, st2, st3, st4))
#            string = ' '.join((string0, ':', stringD))
#            if show_on_screen is True:
#                print string
#            peakstrings.append(string)
#
    def start_at_arbitrary(self, wn_matchup=3000, offset=0.):
        """For each spectrum in profile, divide raw absorbance by thickness 
        and set spectra abs_full_cm such that they overlap at the specified 
        wavenumber wn_matchup with specified offset up from zero"""
        for x in self.spectra:
            # Divide by thickness if not done already
            if x.abs_full_cm is None:
                check = x.divide_by_thickness(folder=self.folder)
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

class TimeSeries(Profile):
    def __init__(self, sample=None, fname_list=[], time_hours=[], folder='',
                 thick_microns=None, style_base=styles.style_points):
        self.sample = sample
        self.thick_microns = thick_microns
        self.fname_list = fname_list
        self.times_hours = time_hours
        self.style_base = style_base        
        self.folder = folder
        self.make_spectra()
    
    def plot_timeseries(self, y=None, peak_idx=None, tit=None, D_list=[], 
                        thickness_microns=None, max_hours=None,
                        style=None, idx_C0=0, figaxis=None): 
        """Plot and return figure and axis of time-series data and
        all diffusivities in D_list in m2/s
        """
        if max_hours is None:
            max_hours = max(self.times_hours)
            
        if thickness_microns is None:
            thickness_microns = self.thick_microns
            
        if style is None:
            style = self.style_base
    
        if figaxis is None:
            fig = plt.figure()
            fig.set_size_inches(3.2, 3.2)
            ax = fig.add_subplot(111)
            fig.autofmt_xdate()
        else:
            ax = figaxis
            
        ax.set_xlabel('Time (hours)', fontsize=12)
        ax.set_ylabel('Concentration/\nMaximum Concentration', fontsize=12)
        ax.set_xlim(0, max_hours)

        # curves for diffusivities in D_list    
        for D in D_list:
             t, cc = diffusion.diffusionThinSlab(log10D_m2s=D, 
                                        thickness_microns=thickness_microns, 
                                        max_time_hours=max_hours)
             ax.plot(t, cc, '-k', linewidth=1)
    
        # plot area data
        x = self.times_hours
        if y is None:
            if peak_idx is None:
                if len(self.areas_list) < 1:
                    self.make_area_list()
                C = self.areas_list
            else:
                print('not ready for peak-specific yet')
            C0 = C[idx_C0]
            y = np.array(C) / C0
        
        ax.plot(x, y, **style)
        
        if figaxis is None:
            return fig, ax



def subtract_2spectra(list2, wn_high=4000, wn_low=3000):
    """Subtract spectrum 1 from spectrum 0 input as list between given
    wavenumbers (defaults to 4000 and 3000 cm-1). Spectra do not need 
    to be the same length or have the exact same wavenumbers."""
    # initial checks
    if len(list2) != 2:
        print('Takes a list of exactly 2 spectra')
        print('length:', len(list2))
        return
    for x in list2:
        # Check they are spectra
        if isinstance(x, Spectrum) is False:
            print(x, 'is not a Spectrum')
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
        print('Length problem in subtract_spectra. To be dealt with.')
        return
    
    # subtract
#    difference = np.zeros_like(list2[0].base_wn)
#    for inx in range(len(list2[0].base_wn)):
#        difference[inx] = list2[0].abs - list2[1].abs
#    return difference
    

def make_all_specta_lists(classname=Profile):
    for obj in gc.get_objects():
        if isinstance(obj, classname):
            obj.make_spectra()
    

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
    ax8.set_xlabel('Wavenumber (cm$^{-1}$)')
    plt.setp(ax5.get_yticklabels(), visible=True)
    plt.setp(ax9.get_xticklabels(), visible=True, rotation=45)
    for ax in [ax1, ax8]:
        plt.setp(ax.get_yticklabels(), visible=True)
        plt.setp(ax.get_xticklabels(), visible=True, rotation=45)    
    return axis_list

def plotsetup_3x3(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
                  xtickgrid=250, ytickgrid=0.5,
                  fig_size_inches=(6.5, 6.5)):
    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
    fig = plt.figure()
    fig.set_size_inches(fig_size_inches)
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
#        ax.grid()
        ax.set_xlim([xhi, xlo])
        ax.set_ylim([ylo, yhi])
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.setp(ax.get_xticklabels(), rotation=45)    

     # Place label A., B. C.
#    for k in range(len(axis_list)):
#        axis_list[k].text(xhi-(0.02*xhi), yhi-(0.2*yhi), 
#                string.ascii_uppercase[k], fontweight='bold')

    ax3.set_xlabel('Wavenumber (cm$^{-1}$)')
    for ax in [ax1, ax2]:
        plt.setp(ax.get_xticklabels(), visible=False)
    return axis_list

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
        print(fin.wb_areas)
        # initial checks
        if len(fin.positions_microns) == 0:
            print('Need position information')
            return False
        if fin.length_microns is None:
            check = fin.set_len()
            if check is False:
                print('Need more info to set profile length')
                return False
        if fin.length_microns is None:
            fin.set_len()
        if len(fin.areas_list) == 0:
            print('making area list for profile')
            fin.make_area_list(show_plot=False)
    
        # What to normalize to? Priority given to self.wb_initial_profile, then
        # initial_profile passed in here, then initial area list passed in here.
        if fin.initial_profile is not None:
            
            init = fin.initial_profile
            if init.length_microns != fin.length_microns:
                print('initial and final lengths must be the same!')
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas_list) == 0:
                    print(profile.profile_name)
                    print('making area list for profile')
                    profile.make_area_list(show_plot=False)
            A0 = init.areas_list
            positions0 = init.positions_microns           
    
        elif initial_profile is not None:
            init = initial_profile
            if isinstance(init, Profile) is False:
                print('initial_profile argument must be a pynams Profile.')
                print('Consider using initial_area_list and positions instead')
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas_list) == 0:
                    print('making area list for profile')
                    profile.make_area_list(show_plot=False)
            A0 = init.areas_list
            positions0 = init.positions_microns
    
        elif initial_area_list is not None:
            if initial_area_positions_microns is None:
                print('Need initial_area_positions_microns for initial_area_list')
                return False
            A0 = initial_area_list
            positions0 = initial_area_positions_microns
            if len(fin.areas_list) == 0:
                print('making area list for final profile')
                fin.make_area_list(show_plot=False)
        else:
            print('Need some information about initial state')
            return False
    
        # More initial checks
        if len(fin.areas_list) != len(fin.positions_microns):
            print('area and position lists do not match')
            print('length areas_list:', len(fin.areas_list))
            print('length positions list:', len(fin.positions_microns))
            return False    
        if len(A0) < 1:        
            print('Nothing in initial area list')
            return False
        if len(positions0) < 1:
            print('Nothing in initial positions list')
            return False
        if len(A0) == 1:
            print('Using single point to generate initial line')
            A0.extend([A0[0], A0[0]])
            positions0.extend([0, fin.length_microns])
            
        # Use best-fit line through initial values to normalize final data
        p = np.polyfit(positions0-(leng/2.), A0, 1)
        
        normalizeto = np.polyval(p, fin.areas_list)
        wb_areas = fin.areas_list / normalizeto
         
        # Save whole-block areas as part of profile
        fin.wb_areas = wb_areas    
    
    if show_plot is True:
        if fig_ax is None:
            f, ax = styles.plot_area_profile_outline(fin)
        else:
            ax = fig_ax
        ax.set_ylim(0, top)
        ylabelstring = 'Final area / Initial area'
        if peakwn is not None:
            print('NOT READY FOR PEAKFITTING YET')
#            extrabit = '\n for peak at ' + str(peakwn) + ' cm$^{-1}$'
#            ylabelstring = ylabelstring + extrabit
        ax.set_ylabel(ylabelstring)
        style = fin.choose_marker_style()
        ax.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
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
            print('Need initial water content.')
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
        if fin.length_microns is not None:
            leng = fin.length_microns
        else:
            leng = fin.set_len()
        ax_areas.set_xlim(-leng/2.0, leng/2.0)
            
        style = fin.choose_marker_style()
        ax_areas.plot([-leng/2.0, leng/2.0], [1, 1], **styles.style_1)
        ax_areas.plot(fin.positions_microns-leng/2.0, wb_areas, **style)
        return water, fig, ax_areas
    else:
        return water

                
