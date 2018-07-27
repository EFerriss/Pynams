"""
Code for processing and plotting *groups* of FTIR spectra.

The central concept is an object class called Profile, which creates groups 
of Spectrum objects, which are defined in a separate module.

@author Elizabeth Ferriss
"""

from __future__ import print_function, division, absolute_import
import numpy as np
from . import styles
from . diffusion import models
from . import pynams
from .spectra import Spectrum
import uncertainties
from uncertainties import ufloat
import gc
import matplotlib.pyplot as plt
import os.path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import string as string
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import json
import lmfit
import copy

class Profile():
    def __init__(self, name=None,
                 fnames=[], 
                 folder='',
                 positions_microns = np.array([]),
                 length_microns=None,
                 thicknesses_microns=None,
                 time_seconds=None, 
                 sample=None, 
                 direction=None, 
                 raypath=None, 
                 initial_profile=None, 
                 base_low_wn=None, 
                 base_high_wn=None,
                 diffusivity_log10m2s=None, 
                 diff_error=None, 
                 peak_diffusivities=[], 
                 peak_diff_error=[], 
                 ):
        """
        Creates a group of FTIR Spectrum objects that can be handled
        and interpreted together.
        
        name is a description automatically used for labeling plots.
        
        fnames must be a *list* of spectra filenames without the .CSV or 
        .txt extension, just like the fname when making a Spectrum.
        
        folder specifies the location of the FTIR files.
        
        positions_microns is a *list* of locations in the same order as the
        fnames
        
        The thickness of the spectra can be specified now with 
        thicknesses_microns, either using a float or integer, or a list 
        of thicknesses, 
        
        Raypath and direction expressed as 'a', 'b', 'c' with thickness/length
        info contained in sample's length_a_microns, length_b_microns, 
        and length_c_microns.
        
        """
        self.name = name
        self.folder = folder
        self.fnames = fnames
        self.positions_microns = positions_microns
        self.sample = sample
        self.length_microns = length_microns
        self.direction = direction
        self.raypath = raypath
        self.initial_profile = initial_profile
        self.time_seconds = time_seconds
        self.diffusivity_log10m2s = diffusivity_log10m2s
        self.diff_error = diff_error
        self.peak_diffusivities = peak_diffusivities
        self.peak_diff_error = peak_diff_error
        
        try:
            if (self.raypath is not None) and (self.raypath == self.direction):
                print("raypath cannot be the same as profile direction")
                return False
        except AttributeError:
            self.raypath = None

        # Create list of spectra in profile.spectra
        # construct each spectrum from fnames
        fspectra_list = []
        for x in self.fnames:
            newspec = Spectrum(fname=x, folder=self.folder)
            newspec.fname = x
            fspectra_list.append(newspec)
        self.spectra = fspectra_list

        # set sample, raypath, low and high baseline wavenumber for all
        for spec in self.spectra:
            spec.sample = self.sample
            spec.raypath = self.raypath

        # set thickness from direct input or sample + raypath
        if thicknesses_microns is not None:
            self.thicknesses_microns = thicknesses_microns
        elif sample is not None and raypath is not None:
            idx = styles.get_iorient(raypath)
            self.thicknesses_microns = np.mean(sample.thicknesses_microns[idx])
        else:
            self.thicknesses_microns = None

        # set thicknesses if thickness information is available
        self.update_spectra_thicknesses_from_profile()
            
        # set or guess profile length
        if self.length_microns is None:
            if (self.sample is not None) and (self.direction is not None):
                self.set_length_from_sampledirection()
            else:
                if self.positions_microns is not None:
                    maxlen = max(list(self.positions_microns))
                    self.length_microns = maxlen + 0.1*maxlen
                
                
    def get_thicknesses_from_SiO(self, show_plots=False, printout=False,
                               accept_thickness=True):
        """
        Individually set thicknesses for all spectra based on the area
        under their Si-O overtone peaks
        """
        for spec in self.spectra:
            spec.get_thickness_from_SiO(show_plot=show_plots,
                                        printout=printout,
                                        accept_thickness=accept_thickness)
        if accept_thickness is True:
            self.update_profile_thicknesses_from_spectra()


    def update_profile_thicknesses_from_spectra(self):
        """
        Gathers thicknesses from spectra into list profile.thicknesses_microns
        """
        thicknesses = []
        for spec in self.spectra:
            try:
                thick = spec.thickness_microns
            except AttributeError:
                print('Spectrum has no thickness. Cannot update profile.')
                return
            thicknesses.append(thick)            
            self.thicknesses_microns = thicknesses        


    def update_spectra_thicknesses_from_profile(self):
        """
        Updates spectra thickness information based on what is in 
        profile.thicknesses_microns
        """
        th = self.thicknesses_microns
        if (isinstance(th, float)) or (isinstance(th, int)):
            for spec in self.spectra:
                spec.thickness_microns = th
        elif isinstance(th, list):
            if len(th) == len(self.fnames):
                for idx, spec in enumerate(self.spectra):
                    spec.thickness_microns = th[idx]
            else:
                print('length of thicknesses_microns not equal length fnames')

            
    def set_length_from_sampledirection(self):
        """
        Set profile.length_microns from profile.direction and 
        profile.sample.thickness_microns
        """ 
        if self.sample is None:
            print('\n', self.name)
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


    def set_thickness_from_sample(self):
        """
        Set profile.thickness_microns from profile.raypath and
        profile.sample.thickness_microns
        """
        if self.sample is None:
            print('Need to specify profile sample or thickness')
            return False
        else:
            s = self.sample

        if self.raypath == 'a':
           self.thickness_microns = s.thickness_microns[0]
        elif self.raypath == 'b':
            self.thickness_microns = s.thickness_microns[1]
        elif self.raypath == 'c':
            self.thickness_microns = s.thickness_microns[2]
        else:
            print('Need raypath')
            return False
            
        return self.thickness_microns


    def plot_thicknesses(self, axes=None):
        """Plot thickness across profile"""
        if self.length_microns is None:
            self.length_microns = max(self.positions_microns)+1.
        if axes is not None:
            ax = axes
        else:
            fig, ax, ax_right = styles.plot_area_profile_outline(
                                                            centered=False)
        ax.plot(self.positions_microns, self.thickness_microns_list, 'o')
        ax.set_ylabel('thickness ($\mu$m)')
        ax.set_title(self.name)
        ax.set_ylim(min(self.thickness_microns_list)-0.05*min(self.thickness_microns_list), 
                    max(self.thickness_microns_list)+0.05*max(self.thickness_microns_list))
        return fig, ax            


    def average_spectra(self, ignore_idx=[]):
        """
        Creates and returns a spectrum that is the average of all spectra
        in the profile.
        
        Pass a list of the indexes of any spectra that you don't included 
        to ignore_idx 
        """
        spec_list = copy.deepcopy(self.spectra)
        
        # ignore certain data points
        if ignore_idx is None:
            pass
        elif type(ignore_idx) != list:
            print('ignore_idx takes a list')
        elif len(ignore_idx) == 0:
            pass
        else:
             spec_list = [i for j, i in enumerate(spec_list) if j not in ignore_idx]
             
        avespec = Spectrum(folder=None, fname=None)
        avespec.make_average_spectra(spec_list, folder=self.folder)
        
        if self.name is not None:
            avespec.fname = (self.name + '\naveraged across profile')
        else:
            avespec.fname = 'averaged across profile'
            
        avespec.thickness_microns = np.mean(self.thicknesses_microns)
        return avespec


    def plot_spectra(self):
        """
        Plot all spectra across profile with default spectrum.plot_spectrum()
        method. For more detailed inputs or to return figure and axes handles,
        loop through profile.spectra
        """
        for spec in self.spectra:
            spec.plot_spectrum()
            
            
    def plot_showbaselines(self, axes=None, style=styles.style_spectrum, 
                          style_base=styles.style_baseline,
                          label=None, label_baseline=False,
                          offset=0.0):
        """
        Plot all spectra with baselines in profile
        """
        for spec in self.spectra:
                spec.plot_showbaseline(axes=axes, style=style, 
                                       style_base=style_base, label=label,
                                       label_baseline=label_baseline,
                                       offset=offset)
            
            
    def plot_subtractbaselines(self):
        """
        Plot all spectra.plot_subtractbaseline() with default settings
        """
        for spec in self.spectra:
            spec.plot_subtractbaseline()


    def change_baseline(self, highwn=3800, lowwn=3000, shift=None):
        """Change baseline parameters for all spectra, final and initial"""
        for spectrum in self.spectra + self.initial_profile.spectra:
            spectrum.base_high_wn = highwn
            spectrum.base_low_wn = lowwn
            if shift is not None:
                spectrum.base_mid_yshift = shift


    def make_composite_peak(self, peak_idx_list):
        """
        Takes a list of peak indices and makes a new composite
        peak for all spectra in profile        
        """
        for spec in self.spectra:
            spec.make_composite_peak(peak_idx_list)
        self.get_peak_info()


    def make_baselines(self,
                      raw_data=False, 
                      wn_low=3200, 
                      wn_high=3700, 
                      linetype='line', 
                      spline_kind='cubic', 
                      spline_wn_low=3000,
                      spline_wn_high=4000,
                      curvature=None, 
                      force_through_wn=None,
                      polynomial_order=None,
                      show_fit_values=False, 
                      show_plot=False,
                      abs_high=None, 
                      abs_low=None,
                      abs_smear_high=0, 
                      abs_smear_low=0,
                      store_baseline=True
                      ):
        """Make baselines for all final and initial spectra in profile"""
        for spectrum in self.spectra:
            spectrum.base_high_wn = wn_high
            spectrum.base_low_wn = wn_low
            spectrum.make_baseline(raw_data=raw_data, wn_low=wn_low, 
                                   wn_high=wn_high, linetype=linetype, 
                                   spline_kind=spline_kind, 
                                   spline_wn_high=spline_wn_high,
                                   spline_wn_low=spline_wn_low,
                                   curvature=curvature,
                                   force_through_wn=force_through_wn,
                                   polynomial_order=polynomial_order,
                                   show_fit_values=show_fit_values,
                                   show_plot=show_plot, abs_high=abs_high,
                                   abs_low=abs_low, 
                                   abs_smear_high=abs_smear_high,
                                   abs_smear_low=abs_smear_low,
                                   store_baseline=store_baseline)


    def get_baselines(self, folder=None, delim=',', 
                     baseline_ending='-baseline.CSV',
                     print_confirmation=True):
        """
        Get previously saved baselines for all spectra in profile.
        Automatically updates the profile.areas for these baselines.
        """
        for spectrum in self.spectra:
            spectrum.get_baseline(folder=folder, delim=delim, 
                                  baseline_ending=baseline_ending,
                                  print_confirmation=print_confirmation)
        self.make_areas()

                        
    def save_baselines(self, folder=None, delim=',',
                      baseline_ending='-baseline.CSV'):
        """Save all baselines in profile"""
        for spectrum in self.spectra:
            spectrum.save_baseline(folder=folder, delim=delim,
                                   baseline_ending=baseline_ending)
        

    def make_areas(self, peak=None, show_plot=False, 
                       printout_area=False, wn_low=None, wn_high=None):
        """
        Make list of areas under the curve for an FTIR profile.
        Default is bulk area. Set peak=wavenumber for peak-specific profile.
        Change wn_low and wn_high for wavenumber range other than 
        full baseline range, as in spectrum.make_area.
        
        If show_plot and printout_area are set to True, then the plot and
        areas will show up when it calculates the areas under the curve
        for each spectrum.
        """
        areas = []
        if peak is None:
            for spec in self.spectra:
                a = spec.make_area(show_plot=show_plot, 
                                              printout=printout_area,
                                              wn_low=wn_low, wn_high=wn_high)
                areas.append(a)            
            self.areas = np.array(areas)
            
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

        
    def make_peakheights(self, peaks=[3600, 3525, 3356, 3236]):
        """
        Requires a list of peak wavenumber locations in cm-1
            (default peaks=[3600, 3525, 3356, 3236])
        Creates or overwrites profiles peak positions and peak_heights using 
        that baseline, not heights from gaussian curves
        """
        self.peakpos = peaks
        self.peak_heights = []
        for peak in (peaks):
            heights = []
            for spec in self.spectra:
                idx = np.abs(peak - spec.base_wn).argmin()
                height_base = spec.base_abs[idx]
                idx = np.abs(peak - spec.wn_full).argmin()
                height_abs = spec.abs_full_cm[idx]                        
                height = height_abs - height_base
                heights.append(height)
            self.peak_heights.append(heights)


    def make_peakfit_like(self, spectrum):
        """
        Exactly like spectrum.make_peakfit_like() but applies the 
        input spectrum's peakfitting information to all spectra in 
        the profile. Any deviations along the profile will have
        to be handfit spectrum by spectrum for now. To do this for
        the first individual spectrum within the profile:
            profile.spectra[0].make_peakfit()
            profile.spectra[0].save_peakfit()
        and then change the index from 0 to whatever spectrum
        you want to look at next
        """
        for spec in self.spectra:
            spec.make_peakfit_like(spectrum)

            
    def save_peakfits(self, folder=None, peak_ending='-peakfit.CSV'):
        """
        Saves all peakfits in the file. Same keywords as 
        spectrum.save_peakfit
        """
        for spec in self.spectra:
            spec.save_peakfit(folder=folder, peak_ending=peak_ending)


    def make_peakfits(self, sensitivity=40, peak_positions=None,
                      peak_heights=None, peak_widths=None,
                      show_plots=False):
        """ 
        Applies Spectrum method make_peakfit() to all spectra in profile
        """
        for spec in self.spectra:
            spec.make_peakfit(sensitivity=sensitivity, 
                              peak_positions=peak_positions,
                              peak_heights=peak_heights, 
                              peak_widths=peak_widths,
                              show_plot=show_plots)

        self.get_peak_info()
        

    def get_peakfits(self, peak_ending='-peakfit.CSV'):
        """
        Get saved fit peak information for all spectra in profile and
        any initial profile that have been saved using the 
        spectrum.save_peakfit or profile.save_peakfits methods.
        
        Change the peak_ending keyword (default='-peakfit.CSV') to 
        specify a particular saved peakfit.
        """
        for spectrum in self.spectra:
            spectrum.get_peakfit(peak_ending=peak_ending)
        try:
            for spectrum in self.initial_profile.spectra:
                spectrum.get_peakfit(peak_ending=peak_ending)
        except AttributeError:
            pass

        if len(self.spectra) > 0:        
            self.get_peak_info()

           
    def get_peak_info(self):
        """
        Pull peak info from individual spectra into a single profile
        attribute
        """
        try:
            peakpos = self.spectra[0].peakpos
        except AttributeError:
            return

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
                try:
                    a.append(spectrum.peak_areas[p])
                except AttributeError:
                    print('Problem with peak areas')
                    a = 0.
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
        print('\n', self.name)

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
            print('\n', self.name)
            print('peak positions (cm-1)\n', self.spectra[0].peakpos)
            print('average peak areas (cm-2)\n', average_peakareas)
            print('average peak heights (cm-1)\n', average_peakheights)
            print('summed areas (cm-2)', total_area)
            print('max peak areas (cm-2)\n', max_peakareas)
            print('max peak heights (cm-1)\n', max_peakheights)
            
        return average_peakareas, average_peakheights, total_area, \
                max_peakareas, max_peakheights


    def plot_peakfits(self, initial_too=False, legloc=1, ytop=1.2):
        """Show fit peaks for all spectra in profile"""
        if len(self.spectra) < 1:
            self.make_spectra()
            
        for spectrum in self.spectra:
            spectrum.plot_peakfit(legloc=legloc, ytop=ytop)
        
        if initial_too is True and self.initial_profile is not None:
            for spectrum in self.initial_profile.spectra:
                spectrum.plot_peakfit(legloc=legloc)
                

    def make_wholeblock(self, peakfit=False, show_plot=False,
                        heights_or_areas='both'):
        """
        Take initial and final profiles and make self.wb_areas.
        
        If peakfit=True, then make whole-block ratios for heights and
        peak areas (heights_or_areas='both' by default) or just 'heights' or
        'areas'. The results are stored in peak_wb_areas, peak_wb_heights, 
        peak_wb_widths
        """
        if self.initial_profile is None:
            self.initial_profile = self
        init=self.initial_profile

        if peakfit is False:        
            iareas = self.initial_profile.make_areas()
            areas = self.make_areas()
    
            # best-fit line through initial areas        
            if ((len(iareas) == 2) and 
                (init.positions_microns[0] == init.positions_microns[1])):
                p = (0, np.mean(init.areas))
            elif len(iareas) > 1:
                p = np.polyfit(init.positions_microns, iareas, 1)
            else:
                try:
                    p = (0, iareas[0])
                except IndexError:
                    p = (0, 1)
            init_line = np.polyval(p, self.positions_microns)
            
            # whole-block areas relative to that initial line
            area_ratio = self.areas / init_line            
            self.wb_areas = area_ratio
            return

        # Peak-specific whole-block values
        ipos = self.initial_profile.positions_microns
        pos = self.positions_microns       
        
        if heights_or_areas == 'heights' or heights_or_areas == 'both':
            wb_heights = []
            if self.peak_heights is None:
                self.get_peakfit()
                
            try:
                heights = self.peak_heights
            except AttributeError:
                self.get_peakfit()
            try:
                iheights = self.initial_profile.peak_heights
            except AttributeError:
                self.initial_profile.get_peakfit()
                iheights = self.initial_profile.peak_heights
                
            for height, iheight in zip(heights, iheights):
                ph = np.polyfit(ipos, iheight, 1)
                hnormalizeto = np.polyval(ph, pos)
                wb_heights.append(height / hnormalizeto)
            self.peak_wb_heights = wb_heights
                
        if heights_or_areas == 'areas' or heights_or_areas == 'both':
            wb_areas = []
            try:
                areas = self.peak_areas
            except AttributeError:
                print('need profile.peak_areas')
                return
            iareas = self.initial_profile.peak_areas
            for area, iarea in zip(areas, iareas):
                pa = np.polyfit(ipos, iarea, 1)
                anormalizeto = np.polyval(pa, pos)
                wb_areas.append(area / anormalizeto)                
            self.peak_wb_areas = wb_areas            


    def get_peak_wb_areas(self, peak_idx=0, peakwn=None, 
                          heights_instead=False):
        """Returns peak-specific whole-block areas for the profile
        AND peak wavenumber in cm-1 because usually the peak index is easier
        to pass.
        Defaults to the first peak in the peak position list"""
        # Add check that argument is actually a profile
        
        if self.peakpos is None: 
            print('Getting peaks fit in matlab for', self.name)
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
        if heights_instead is True:
            h_or_a = 'heights'
        elif heights_instead is False:
            h_or_a = 'areas'
        else:
            h_or_a = 'both'
        self.make_wholeblock(peakfit=True, show_plot=False,
                             heights_or_areas=h_or_a)
               
        if heights_instead is False:
            returnvals = self.peak_wb_areas[peak_idx]
        else:
            returnvals = self.peak_wb_heights[peak_idx]
    
        return returnvals, peakwn


    def plot_area_profile(self, 
                          axes=None, 
                          show_water_ppm=True,
                          wholeblock=False,
                          peak_idx=None, 
                          centered=False, 
                          normalize_areas=False,
                          normalize_positions=False,
                          heights_instead=False, 
                          bestfitline=False, 
                          style_bestfitline=None,
                          show_FTIR=False, 
                          show_values=False, 
                          peakwn=None, 
                          ytop=None,
                          style=None, 
                          show_initial_areas=False,
                          error_percent=0, 
                          label=None, 
                          initial_style=None,
                          initial_label=None, 
                          phase='olivine',
                          calibration='Bell', 
                          scale_water=3.,
                          shift_water=0., 
                          ignore_idx=[]):
        """
        Plots the area profile. 
        
        If axes is not None, it will plot the area profile to axes. 
        
        If axes is None (default), it will return three things: 
            the figure handle, the right-axes handle showing the absorbance,
            and the left-axes handle showing the rough estimate of the water.
            
        If show_water_ppm is set to False, the left-hand axis will not
        show up, and the third item returned will be None instead of an 
        axes handle.
        
        The water estimate defaults to assuming the phase='olivine' and 
        the calibration='Bell' for the work of Dave Bell, but you can use
        any phase and calibration accepted by the function 
        pynams.absorption_coefficients(). The scale_diffusion (default=3)
        is what the resulting water concentration estimated from the single
        profile is multiplied by to get a rough water estimate. These 
        estimates should be viewed with great skepticism, but it's a place to 
        start. George Rossman's 2006 discussion of how to quantify water in 
        a Reviews in Mineralogy volume is a good place to start for 
        background on these issues. You can shift that scale up or down
        using shift_water (default=0)
        
        Note that the right and left-hand axes are independent of each other,
        and changing the limits on one will not affect the other. To set
        the maximum area (right-hand axis upper limit) and scale both 
        axes simultaneously while making the plot, use the ytop keyword.
        
        Centered=True (False is default) puts 0 in the middle of the x-axis.
        
        If normalize_areas=True (False is default), areas/heights are 
        normalized to the maximum value. If normalize_positions = True
        (default=False), the positions are plotted relative to the length,
        so with resulting x-axis limits from 0 to 1.

        bestfitline=True draws a best-fit line through the data.
        
        Set peak_idx and whole_block for peak-specific and whole-block 
        profiles.
        
        The keyword ignore_idx allows you to pass in a list of indexes
        for points in your profile to exclude from the plot.
        """
        
        if len(self.positions_microns) < 1:
            print('Need positions_microns for profile')
            return

        # Get y data
        if wholeblock is False:
            if peakwn is None and peak_idx is None:
                try:
                    areas = self.areas
                except AttributeError:
                    areas = self.make_areas(peak=None)
            else:
                self.get_peak_info()

                if peak_idx is None:
                    peak_idx = np.where(self.peakpos==peakwn)[0][0]
                    print('peak at', peakwn, 'is index', peak_idx)
                else:
                    peakwn = self.peakpos[peak_idx]

            if heights_instead is True and peak_idx is not None:
                areas = self.peak_heights[peak_idx]
            elif peak_idx is not None:
                areas = self.peak_areas[peak_idx]
        else:
            if self.initial_profile is None:
                print('Need to specify an initial profile')
                return
            if peak_idx is None:
                if self.wb_areas is None:
                    self.make_wholeblock(peakfit=False, bulk=True)
                areas = self.wb_areas
            else:
                if self.peak_wb_areas is None:
                    self.make_wholeblock(peakfit=True, bulk=False)
    
                if heights_instead is False:
                    areas = self.peak_wb_areas[peak_idx]
                else:
                    areas = self.peak_wb_heights[peak_idx]                   
        if areas is False:
            return            
        
        if normalize_areas is True:
            areas_for_water = list(areas)
            areas = np.array(areas) / float(max(areas))
        else:
            areas_for_water = areas
            
        if np.shape(areas) != np.shape(self.positions_microns):
            print('Area and positions lists are not the same size!')
            print('area:', np.shape(self.areas))
            print('positions:', np.shape(self.positions_microns))
            return

        # Set length
        if self.length_microns is None:
            if self.sample is not None:
                leng = self.set_len()
            else:
                longest = max(self.positions_microns)
                leng = longest + 0.1*longest
        else:
            leng = self.length_microns

        # Set up plotting styles
        if style_bestfitline is None:
            style_bestfitline = styles.style_baseline
        if label is None:
            if style is not None:
                style['label'] = self.name
            else:
                label = self.name
        else:
            if style is not None:
                style['label'] = label

          
        # x data
        if normalize_positions is True:
            x = np.array(self.positions_microns) / leng
            if centered is True:
                x = x - 0.5
        else:
            if centered is True:
                x = np.array(self.positions_microns) - leng/2.0
            else:
                x = self.positions_microns            
                
        # ignore certain data points
        if ignore_idx is None:
            pass
        elif type(ignore_idx) != list:
            print('ignore_idx takes a list of integers: indices to ignore')
        elif len(ignore_idx) == 0:
            pass
        else:
             areas = [i for j, i in enumerate(areas) if j not in ignore_idx]
             areas_for_water = [i for j, i in enumerate(areas_for_water) if j not in ignore_idx]
             x = [i for j, i in enumerate(x) if j not in ignore_idx]
             
        # Use new or old figure axes
        if axes is None:
            if ytop is None:
                if peak_idx is None:
                    ytop = max(areas)+0.2*max(areas)
                else:
                    if heights_instead is False:
                        mpeak = max(self.peak_areas[peak_idx])
                    else:
                        mpeak = max(self.peak_heights[peak_idx])
                    ytop = mpeak + 0.2*mpeak
                       
            f, ax, ax_ppm = styles.plot_area_profile_outline(centered, 
                            peakwn, heights_instead=heights_instead, 
                            wholeblock=wholeblock, ytop=ytop,
                            show_water_ppm=show_water_ppm)

            if normalize_areas is True:
                if heights_instead is False:
                    ax.set_ylabel('Area normalized to max. area')
                else:
                    ax.set_ylabel('Height normalized to max. height')
            if normalize_positions is True:
                ax.set_xlabel('Position relative to total length')
                if centered is True:
                    ax.set_xlim(-0.5, 0.5)
                else:
                    ax.set_xlim(0., 1.)
            else:
                if centered is True:
                    ax.set_xlim(-leng/2.0, leng/2.0)
                else:
                    ax.set_xlim(0, leng)
        else:
            ax = axes
            show_water_ppm = False

        # plot
        yerror = np.array(areas)*error_percent/100.       
        if error_percent == 0:
            if style is None:
                ax.plot(x, areas, label=label)
            else:
                ax.plot(x, areas, **style)
        else:
            if style is not None:
                ax.errorbar(x, areas, yerr=yerror, **style)
            else:
                ax.errorbar(x, areas, yerr=yerror, label=label)

            
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

        # Title
        if peak_idx is None:
            tit = 'Bulk hydrogen'
        else:
            peakwn = self.peakpos[peak_idx]
            tit = ' '.join(('Peak at', str(peakwn) ,'/cm'))
        ax.set_title(tit)

        if show_water_ppm is True:
            abs_coeff = pynams.absorption_coefficients(phase=phase, 
                                                       calibration=calibration)
            if abs_coeff is False:
                show_water_ppm = False

        # The water concentration estimate on the right-hand axis
        if heights_instead is True:
            show_water_ppm = False
        if show_water_ppm is True:
            if shift_water == 0.:
                ax_ppm.set_ylabel(''.join(('ppm H2O in ', phase, ', ', 
                                           calibration, ' calibration *', 
                                           str(scale_water))))
            else:
                ax_ppm.set_ylabel(''.join(('ppm H2O in ', phase, ', ', 
                                           calibration, ' calibration *', 
                                           str(scale_water), ' + ',
                                           str(shift_water))))                
            ori = scale_water
            if normalize_areas is False:
                ppm_limits = np.array(ax.get_ylim()) * abs_coeff.n * ori
            else:
                ppm_limits = np.array([0, max(areas_for_water)*ytop])
                ppm_limits = ppm_limits * abs_coeff.n * ori
            ax_ppm.set_ylim(ppm_limits + shift_water)
        elif axes is None:
            ax_ppm.get_yaxis().set_visible(False)

        if axes is None:
            return f, ax, ax_ppm
        else:
            return


    def diffusion1D(self, log10D_m2s, time_seconds, wholeblock=False,
                    peak_idx=None, centered=True, heights_instead=False, 
                    init=1., fin=0., erf_or_sum='erf', infinity=100, 
                    points=100, symmetric=True, maximum_value=None):
        """
        Requires log10D_m2s and time_seconds
        
        Other keywords like pynams.diffusion.models.diffusion1D
        and profile.plot_diffusion()
        
        Returns x and y for diffusion modeling, without plotting
        """
        try:
            length = self.length_microns
        except AttributeError:
            length = max(self.positions_microns)
        
        if maximum_value is None:
            if peak_idx is None:
                try:
                    maximum_value = max(self.areas)
                except AttributeError:
                    maximum_value = max(self.make_areas())
            else:
                if heights_instead is True:
                    maximum_value = max(self.peak_heights[peak_idx])
                else:
                    maximum_value = max(self.peak_areas[peak_idx])

        fig, ax, x, y = models.diffusion1D(length, log10D_m2s, time_seconds, 
                                           init=init, fin=fin, 
                                           erf_or_sum=erf_or_sum, 
                                           show_plot=False, infinity=infinity, 
                                           points=points, centered=centered,
                                           symmetric=symmetric, 
                                           maximum_value=maximum_value)
        return x, y*maximum_value


    def plot_diffusion(self, log10D_m2s, time_seconds, 
                      axes=None, 
                      show_water_ppm=True,
                      wholeblock=False,
                      peak_idx=None, 
                      centered=True, 
                      heights_instead=False, 
                      style=styles.style_points, 
                      style_diffusion_line=styles.style_blue,
                      phase='olivine',
                      calibration='Bell',
                      scale_water=3,
                      shift_water=0.,
                      init=1., 
                      fin=0.,
                      erf_or_sum='erf', 
                      infinity=100, 
                      points=100, 
                      symmetric=True,
                      maximum_value=None,
                      ytop=None,
                      normalize_areas=False,
                      normalize_positions=False):

        """ 
        Plots area profile and with 1D diffusion profile on top.
    
        Returns figure handle and both left and right axes handles
        like profile.plot_area_profile()
        
        Includes most of the keywords from profile.plot_area_profile
        and diffusion.models.diffusion1D
        """
        fig, ax, ax_water = self.plot_area_profile(axes=axes,
                                       show_water_ppm=show_water_ppm,
                                       wholeblock=wholeblock,
                                       peak_idx=peak_idx,
                                       centered=centered,
                                       heights_instead=heights_instead,
                                       style=style,
                                       phase=phase,
                                       calibration=calibration,
                                       scale_water=scale_water,
                                       shift_water=shift_water,
                                       ytop=ytop,
                                       normalize_areas=normalize_areas,
                                       normalize_positions=normalize_positions)
        try:
            length = self.length_microns
        except AttributeError:
            length = max(self.positions_microns)
        
        if maximum_value is None:
            if peak_idx is None:
                try:
                    maximum_value = max(self.areas)
                except AttributeError:
                    maximum_value = max(self.make_areas())
            else:
                if heights_instead is True:
                    maximum_value = max(self.peak_heights[peak_idx])
                else:
                    maximum_value = max(self.peak_areas[peak_idx])

        models.diffusion1D(length, log10D_m2s, time_seconds, init=init, 
                           fin=fin, erf_or_sum=erf_or_sum, show_plot=True, 
                           style=style_diffusion_line, infinity=infinity, 
                           points=points, centered=centered, axes=ax, 
                           symmetric=symmetric, maximum_value=maximum_value)
       
        return fig, ax, ax_water
    

    def y_data_picker(self, wholeblock, heights_instead, peak_idx=None):
        """
        Pick out and return peak area or height data of interest.
        """
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
            if peak_idx is None:
                try:
                    y = self.areas
                except AttributeError:
                    y = self.make_areas()
            else:
                try:
                    y = self.peak_areas[peak_idx]
                except AttributeError:
                    print('Need peak areas')
                    return False
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
        """
        Saves a diffusivity in the profile attribute of interest.
        Consider using get_diffusivity() first
        """
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


    def scale_diffusion_picker(self, maximum_value=None, wholeblock=True,
                              heights_instead=False, peak_idx=None):
        """Returns value to scale diffusion to"""
        if maximum_value is not None:
            scale_diffusion = maximum_value
            
        elif wholeblock is True:
            if peak_idx is None:
                scale_diffusion = self.maximum_wb_area
            else:
                if heights_instead is False:
                    scale_diffusion = self.peak_maximum_areas_wb[peak_idx]
                else:
                    scale_diffusion = self.peak_maximum_heights_wb[peak_idx]
        
        else: 
            # bulk water                
            if peak_idx is None:
                if len(self.areas) < 1:
                    self.make_areas()
                scale_diffusion = np.max(self.areas)
        
            # peak-specific heights                
            elif heights_instead is True: 
                if self.peak_maximum_heights[peak_idx] != 0.0:
                    scale_diffusion = self.peak_maximum_heights[peak_idx]
                else:                
                    scale_diffusion = np.max(self.peak_heights[peak_idx])
                    self.peak_maximum_heights[peak_idx] = scale_diffusion
            # peak-specific areas
            else:
                if self.peak_maximum_areas[peak_idx] != 0.0:
                    scale_diffusion = self.peak_maximum_areas[peak_idx]
                else:
                    scale_diffusion = np.max(self.peak_areas[peak_idx])
                    self.peak_maximum_areas[peak_idx] = scale_diffusion
    
        return scale_diffusion


    def diffusion_residuals(self, time_seconds, log10D_m2s, wholeblock=False,
                            heights_instead=False, peak_idx=None,
                            initial_unit_value=1., final_unit_value=0.,
                            maximum_value=None):
        """
        Compares 1D diffusion curve with profile data.
        Returns a vector containing the residuals 
        and the variance, which is sqrt(sum of squares of the residuals)
        """
        if self.positions_microns is None:
            print('Need to set profile positions')
            return

        y = self.y_data_picker(wholeblock, heights_instead, peak_idx)
        scale_diffusion = self.scale_diffusion_picker(maximum_value, 
                                        wholeblock, heights_instead, peak_idx)
        
        x = self.positions_microns
        L = self.length_microns
        t = time_seconds

        # Need to work on this        
        init = scale_diffusion
            
        fin = final_unit_value
        
        params = models.params_setup1D(L, log10D_m2s, t, init, fin, 
                                          False, False, False)

        xdif, model = models.diffusion1D_params(params)
        resid = models.diffusion1D_params(params, x, y)
        RSS = np.sum(resid**2)        
        return resid, RSS
            

    def fitD(self, 
             time_seconds=None, 
             log10Dm2s=None, 
             initial_unit_value=1., 
             final_unit_value=0.,  
             starting_value=None, 
             vary_initial=False,
             vary_final=False,
             vary_time=True, 
             varyD=False, 
             peak_idx=None, 
             wholeblock=False, 
             centered=False,
             show_plot=True, 
             heights_instead=False,
             symmetric=True,
             points=200,
             ignore_idx=[]):
        """
        Fits a 1D diffusion curve to profile data.
        
        Prints and returns: initial, final, log10 diffusivity in m2/s, 
        time in minutes.

        Default fits both diffusivity and initial concentration and
        holds time (default=10 hours) and final unit value constant.
        Set vary_initial=False to hold initial constant. Set 
        vary_time and vary_final to True to fit those. Set varyD = False
        to hold diffusivity constant, and set the diffusivity with the
        keyword log10D_m2s (default is -12).
        
        Default initial is the maximum y-value, but you can change that
        with initial_unit_value or starting_value.
        
        Pass list of indexes of points you don't want included to ignore_idx
        """           
        if time_seconds is None:
            try:
                time_seconds = self.time_seconds
            except AttributeError:
                time_seconds = 36000.
            if self.time_seconds is None:
                time_seconds = 36000.

        if log10Dm2s is None:
            try:
                log10Dm2s = self.D_area
            except AttributeError:
                log10Dm2s = -12.

        if self.length_microns is None:
            print('Need to set profile attribute length_microns')
            return
            
        if self.positions_microns is None:
            print('Need to set profile positions')
            return
        
        x = self.positions_microns.copy()
        y_raw = self.y_data_picker(wholeblock, heights_instead, peak_idx).copy()

        # ignore certain data points
        if ignore_idx is None:
            pass
        elif type(ignore_idx) != list:
            print('ignore_idx takes a list')
        elif len(ignore_idx) == 0:
            pass
        else:
             x = [i for j, i in enumerate(x) if j not in ignore_idx]
             y_raw = [i for j, i in enumerate(y_raw) if j not in ignore_idx]

        scale_diffusion = max(y_raw)
        y = y_raw / scale_diffusion 

             
        if starting_value is None:
            init = initial_unit_value
        else:
            init = starting_value / max(y_raw)

        params = models.params_setup1D(microns=self.length_microns, 
                                       log10D_m2s=log10Dm2s, 
                                       time_seconds=time_seconds, 
                                       init=init,
                                       fin=final_unit_value,
                                       vinit=vary_initial, 
                                       vfin=vary_final,
                                       vTime=vary_time, 
                                       vD=varyD)

        dict_fitting = {'points' : points,
                        'symmetric' : symmetric,
                        'centered' : centered
                        }

        lmfit.minimize(models.diffusion1D_params, params, args=(x, y), 
                       kws=dict_fitting)
        best_D = ufloat(params['log10D_m2s'].value, 
                        params['log10D_m2s'].stderr)
        best_init = ufloat(params['initial_unit_value'].value, 
                         params['initial_unit_value'].stderr)
        best_final = ufloat(params['final_unit_value'].value, 
                         params['final_unit_value'].stderr)
        best_time = ufloat(params['time_seconds'].value,
                           params['time_seconds'].stderr)

        if wholeblock is True:
            if peak_idx is not None:
                if heights_instead is False:
                    self.D_peakarea_wb[peak_idx] = best_D.n
                    self.D_peakarea_wb_error[peak_idx] = best_D.s
                    self.peak_maximum_areas_wb[peak_idx] = best_init.n
                else:
                    self.D_height_wb[peak_idx] = best_D.n
                    self.D_height_wb_error[peak_idx] = best_D.s
                    self.peak_maximum_heights_wb[peak_idx] = best_init.n

            else:
                self.D_area_wb = best_D.n
                self.D_area_wb_error = best_D.s
                self.maximum_wb_area = best_init.n
        
        print('initial:', '{:.2f}'.format(best_init*scale_diffusion))
        print('final:', '{:.2f}'.format(best_final*scale_diffusion))
        print('log10D m2/s:', '{:.2f}'.format(best_D))
        print('minutes:', '{:.2f}'.format(best_time/60))
        if show_plot is True:
            fig, ax, ax2 = self.plot_diffusion(time_seconds=best_time.n,
                                   log10D_m2s=best_D.n,
                                   peak_idx=peak_idx,
                                   wholeblock=wholeblock,
                                   centered=centered,
                                   symmetric=symmetric,
                                   heights_instead=heights_instead, 
                                   maximum_value=best_init.n*scale_diffusion,
                                   fin=best_final.n)
        D = best_D
        init = best_init*scale_diffusion
        fin = best_final*scale_diffusion
        minutes = best_time / 60.
        return init, fin, D, minutes


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

        print('\n', self.name)
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


class TimeSeries(Profile):
    def __init__(self, sample=None, fnames=[], time_hours=[], folder='',
                 thickness_microns=None, style_base=styles.style_points):
        self.sample = sample
        self.thickness_microns = thickness_microns
        self.fnames = fnames
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
            thickness_microns = self.thickness_microns
            
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
             t, cc = models.diffusionThinSlab(log10D_m2s=D, 
                                        thickness_microns=thickness_microns, 
                                        max_time_hours=max_hours)
             ax.plot(t, cc, '-k', linewidth=1)
    
        # plot area data
        x = self.times_hours
        if y is None:
            if peak_idx is None:
                if len(self.areas) < 1:
                    self.make_areas()
                C = self.areas
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
        

def make_all_specta_lists(classname=Profile):
    for obj in gc.get_objects():
        if isinstance(obj, classname):
            obj.make_spectra()
    

def plotsetup_3x3minus2(yhi = 1, ylo = 0, xlo = 3000, xhi = 4000,
                        xtickgrid=250, ytickgrid=0.5):
    """Setup plot for spectra comparison e.g., Kunlun_peakcheck.py"""
    fig = plt.figure()
    fig.set_size_inches(6.5, 6.5)
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
        ax.set_xlim([xhi, xlo])
        ax.set_ylim([ylo, yhi])
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.setp(ax.get_xticklabels(), rotation=45)    

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
        if len(fin.areas) == 0:
            print('making area list for profile')
            fin.make_areas(show_plot=False)
    
        # What to normalize to? Priority given to self.wb_initial_profile, then
        # initial_profile passed in here, then initial area list passed in here.
        if fin.initial_profile is not None:
            
            init = fin.initial_profile
            if init.length_microns != fin.length_microns:
                print('initial and final lengths must be the same!')
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas) == 0:
                    print(profile.name)
                    print('making area list for profile')
                    profile.make_areas(show_plot=False)
            A0 = init.areas
            positions0 = init.positions_microns           
    
        elif initial_profile is not None:
            init = initial_profile
            if isinstance(init, Profile) is False:
                print('initial_profile argument must be a pynams Profile.')
                print('Consider using initial_area_list and positions instead')
                return False
            # Make sure area lists are populated
            for profile in [init, fin]:
                if len(profile.areas) == 0:
                    print('making area list for profile')
                    profile.make_areas(show_plot=False)
            A0 = init.areas
            positions0 = init.positions_microns
    
        elif initial_area_list is not None:
            if initial_area_positions_microns is None:
                print('Need initial_area_positions_microns for initial_area_list')
                return False
            A0 = initial_area_list
            positions0 = initial_area_positions_microns
            if len(fin.areas) == 0:
                print('making area list for final profile')
                fin.make_areas(show_plot=False)
        else:
            print('Need some information about initial state')
            return False
    
        # More initial checks
        if len(fin.areas) != len(fin.positions_microns):
            print('area and position lists do not match')
            print('length areas:', len(fin.areas))
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
        
        normalizeto = np.polyval(p, fin.areas)
        wb_areas = fin.areas / normalizeto
         
        # Save whole-block areas as part of profile
        fin.wb_areas = wb_areas    
    
    if show_plot is True:
        if fig_ax is None:
            f, ax = styles.plot_area_profile_outline()
        else:
            ax = fig_ax
        ax.set_ylim(0, top)
        ylabelstring = 'Final area / Initial area'
        if peakwn is not None:
            print('NOT READY FOR PEAKFITTING YET')
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