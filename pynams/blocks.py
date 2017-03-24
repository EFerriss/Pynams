"""
Code for grouping 3 orthogonal Profiles into a single Block object.
"""
from __future__ import print_function, division, absolute_import
import pynams.styles as styles
from pynams.diffusion import models
from pynams import Spectrum
import numpy as np
import matplotlib.pyplot as plt
import lmfit
from uncertainties import ufloat

class Block():
    def __init__(self,profiles=[], folder='', name='', get_peakfit=False, 
                 make_wb_areas=False, time_seconds=None, sample=None,
                 get_baselines=False, initialWB=None, celsius=None):
        """
        Sets up and checks new Block
        - Check that profiles list contains a list of three (3) profiles
        - Generate list of initial profiles
        - Generate list of profile directions
        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
        - Generate list of ray paths
        - Verify three ray path directions are compatible with directions list
        """
        if len(profiles) != 3:
            print('keyword profiles must be a list of 3 pynams profiles')
            print('If you only have 1 or 2, make a dummy profile')
            return
        
        self.profiles = profiles
        self.folder = folder
        self.name = name
        self.time_seconds = time_seconds
        self.sample = sample        
        self.peak_diffusivities = []
        self.peak_diffusivities_errors = []
        self.celsius = celsius
        
        d = []
        r = []
        ip = []
        L = []
        
        for prof in self.profiles:
            d.append(prof.direction)
            r.append(prof.raypath)
            L.append(prof.length_microns)
            prof.time_seconds = self.time_seconds

            if prof.sample is None and self.sample is not None:
                prof.sample = self.sample
            elif prof.sample is not None and self.sample is None:
                self.sample = prof.sample
            elif prof.sample != self.sample:
                print('Warning: profile sample does not match wb sample')
                
            if prof.positions_microns is None:
                print('profile', prof.profile_name)
                print('Needs positions_microns attribute')

            if prof.initial_profile is None:
                prof.initial_profile = prof
            ip.append(prof.initial_profile)
            
            if make_wb_areas is True:
                check = prof.make_wholeblock(peakfit=get_peakfit)
                if check is False:
                    return False

            if prof.spectra is None:
                prof.make_spectra()
                
            for spec in prof.spectra:
                if self.sample is not None:
                    if prof.raypath == 'a':
                        meanthick = np.mean(self.sample.length_a_microns)
                    elif prof.raypath == 'b':
                        meanthick = np.mean(self.sample.length_b_microns)
                    elif prof.raypath == 'c':
                        meanthick = np.mean(self.sample.length_c_microns)
                    try:
                        spec.thickness_microns = meanthick
                    except NameError:
                        print('Unable to determine thicknesses from sample')
                    
        self.directions = d
        self.raypaths = r
        self.initial_profiles = ip 
        self.lengths = L

        if get_peakfit is True:
            for prof in self.profiles:
                prof.get_peakfit()

        if get_baselines is True:        
            self.get_baselines()

        if initialWB is not None:
            for idx, prof in self.profiles:
                prof.initial_profile = initialWB.profiles[idx]
                
    def get_peakfits(self, peak_ending='-peakfit.CSV'):
        """Get peakfit information for all profiles"""
        for prof in self.profiles:
            prof.get_peakfits(peak_ending=peak_ending)

    def get_baselines(self, initial_too=True, folder=None, delim=',', 
                      baseline_ending='-baseline.CSV',
                      print_confirmation=False):
        """Get baselines for all spectra in whole block"""
        if self.initial_profiles is None:
            self.setupWB()
        for prof in self.profiles:
            for spectrum in prof.spectra:
                spectrum.get_baseline(baseline_ending=baseline_ending,
                                      folder=folder, delim=delim,
                                      print_confirmation=print_confirmation)
        if initial_too is True:
            for prof in self.initial_profiles:
                for spectrum in prof.spectra:
                    spectrum.get_baseline(baseline_ending=baseline_ending,
                                      folder=folder, delim=delim, 
                                      print_confirmation=print_confirmation)
            
    def plot_showbaselines(self):
        """Plot baselines for all spectra in the whole block"""
        for prof in self.profiles:
            for spec in prof.spectra:
                spec.plot_showbaseline()

    def plot_subtractbaselines(self):
        """
        Plot all spectra.plot_subtractbaseline() with default settings
        """
        for prof in self.profiles:
            for spec in prof.spectra:
                spec.plot_subtractbaseline()
            
    def make_areas(self, show_plot=False, printout_area=False, peak=None):
            """
            Make list of areas from all profiles, including whole-block areas.
            """
            for prof in self.profiles:
                prof.make_wholeblock()

    def average_spectra(self):
        """
        Create and return a single spectrum that is an average of all spectra 
        in all three profiles of the Block
        """
        absorbances_per_cm = []
        for prof in self.profiles:
            profile_average_spectrum = prof.average_spectra()
            absorb = profile_average_spectrum.abs_full_cm
            absorbances_per_cm.append(absorb)
        ave_abs = np.mean(absorbances_per_cm, axis=0)
        
        avespec = Spectrum(folder=None, fname=None)
        avespec.abs_full_cm = ave_abs
        avespec.wn_full = prof.spectra[0].wn_full
        avespec.abs_raw = ave_abs
        
        if self.name is not None:
            avespec.fname = (self.name + '\naverage across all profiles')
        else:
            avespec.fname = 'average across all profiles'
            
        return avespec


    def plot_spectra(self, profile_idx=None, show_baseline=True, 
                     show_initial_ave=True,
                     show_final_ave=True, plot_all=False, 
                     initial_and_final_together=False, style=styles.style_spectrum, 
                     stylei=styles.style_initial, wn=None):
        """Plot all spectra in all or specified profile in whole-block"""
        if profile_idx is None:
            proflist = ([self.profiles[0]] + 
                        [self.profiles[1]] + 
                        [self.profiles[2]])
        else:
            proflist = [self.profiles[profile_idx]]
            
        for prof in proflist:
            print(prof.profile_name)
            prof.plot_spectra(show_baseline=show_baseline, 
                  show_initial_ave=show_initial_ave, plot_all=plot_all,
                  show_final_ave=show_final_ave,
                  initial_and_final_together=initial_and_final_together,
                  style=style, stylei=stylei, wn=None)

    def plot_3panels_ave_spectra(self, peak_idx=None, peakwn=None, 
                                 top=1., high=4000, low=3000,
                                style=styles.style_spectrum_red, 
                                stylei=styles.style_initial, show_raypaths=False,
                                figsize=(6., 4), show_initial=True,
                                legloc=5, label='Final', figax3=None):
        """
        Three suplots showing average initial and final spectra in each
        direction
        """
        if self.initial_profiles is None:
            self.setupWB()

        if figax3 is None:
            f, ax3 = plt.subplots(1,3)
            f.set_size_inches(figsize)
            for idx, ax in enumerate(ax3[:3]):
                ax.set_xlim(high, low)
                ax.set_ylim(0., top)
                if show_raypaths is True:
                    raypathstring = ''.join(('profile || ', self.profiles[idx].direction, 
                                             '\nray path || ', self.profiles[idx].raypath))
                    ax.text(3500, top-top*0.22, raypathstring, 
                            backgroundcolor='w', horizontalalignment='center')
            
        else:
            ax3 = figax3
            
        for k in range(3):
            prof = self.profiles[k]
            avespec = prof.average_spectra()
            ax3[k].plot(avespec.wn_full, avespec.abs_full_cm, label=label, 
                        **style)
        
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
                    ax3[k].plot(initspec.wn_full, initspec.abs_full_cm, 
                            **stylei)

        plt.setp(ax3[1].get_yticklabels(), visible=False)
        plt.setp(ax3[2].get_yticklabels(), visible=False)
        ax3[1].set_xlabel('wavenumber (cm$^{-1}$)')
        ax3[0].set_ylabel('absorbance (cm$^{-1}$)')
        

        if show_initial is True:
            ax3[1].legend(loc=legloc)
        if figax3 is None:
            plt.tight_layout()
            plt.gcf().autofmt_xdate()
            tit = ' '.join(('Averaged profiles for', self.name))
            if peak_idx is not None:
                tit = str.join(tit, ', peak at ', str(peakpos[peak_idx]), '/cm')
            ax3[1].set_title(tit, zorder=100) # must come after tight_layout

        plt.subplots_adjust(top=0.85, bottom=0.25)
        if figax3 is None:
            return f, ax3

    def xy_picker(self, peak_idx=None, wholeblock=True, heights_instead=False,
                  centered=True, unit='microns'):
        """Picks out and returns appropriate x and y-data for 3D"""
        positions = []
        y = []
            
        for prof in self.profiles:
            positions.append(prof.positions_microns)

            # Bulk hydrogen            
            if peak_idx is None:
                
                # whole-block
                if wholeblock is True:
                    try:
                        y_to_add = prof.wb_areas
                    except AttributeError:
                        prof.make_wholeblock(peakfit=False, show_plot=False)
                        y_to_add = prof.wb_areas
                
                # absolute areas
                else:
                    if prof.areas is None:
                        prof.make_areas()
                    y_to_add = prof.areas


            # Peak-specific                
            else:
                for idx, spec in enumerate(prof.spectra):
                    if spec.peak_areas is None:
                        spec.get_peakareas()
                    prof.peak_heights[peak_idx][idx]=spec.peak_heights[peak_idx]
                    prof.peak_areas[peak_idx][idx]=spec.peak_areas[peak_idx]

                if wholeblock is True:
                    peak_wb, peakwn = prof.get_peak_wb_areas(peak_idx, 
                                           heights_instead=heights_instead)
                    y_to_add = peak_wb

                else:
                    if heights_instead is False:
                        y_to_add = prof.peak_areas[peak_idx]
                    else:
                        y_to_add = prof.peak_heights[peak_idx]
            y.append(y_to_add)
            
        
        if centered is True:
            a = np.mean(self.profiles[0].sample.length_a_microns) / 2.
            b = np.mean(self.profiles[1].sample.length_b_microns) / 2.
            c = np.mean(self.profiles[2].sample.length_c_microns) / 2.
            halflengths = [a, b, c]
            for idx in range(3):
                positions[idx] = positions[idx] - halflengths[idx]
            
        return positions, y

    def plot_areas_3panels(self, peak_idx=None, axes3=None, centered=False,
                           ytop=None, heights_instead=False, wholeblock=False,
                           xerror=0., yerror=None, pie=True,
                           label4legend=[None, None, None],
                           styles3=[styles.style_points]*3,
                           unit='microns',
                           show_line_at_1=False, 
                           show_errorbars=True):
        """
        Plots areas (default) or ratio of area to initial area 
        (wholeblock=True) for all 3 profiles.  
        
        Returns figure handle and a list of 3 axes handles (default) unless
        axes3 is not equal to None. 
        
        Keywords are similar to those for profile.plot_areas
        """
        self.make_areas()
        
        if peak_idx is not None:
            for prof in self.profiles:
                for spec in prof.spectra:
                    if spec.peak_areas is None:
                        spec.get_peakareas()    
            peakpos = self.profiles[0].peakpos

        positions, y = self.xy_picker(peak_idx=peak_idx, 
                                      wholeblock=wholeblock,
                                      heights_instead=heights_instead, 
                                      centered=centered, unit=unit)
        
        if peak_idx is not None:
            tit = ' '.join(('Peak at', str(peakpos[peak_idx]), '/cm'))
        else:
            tit = 'Bulk hydrogen'

        if ytop is None:
            z = []
            for ynum in y:
                if len(ynum) > 0:
                    z.append(max(ynum))
            ytop = max(z) + 0.1*max(z)

        if unit == 'microns':
            lengths = self.lengths
        elif unit == 'mm':
            lengths = np.array(self.lengths) / 1000.
        else:
            print('unit must be microns (default) or mm')
            return
                
        # Sent positions and areas to plotting command
        if axes3 is not None:
            styles.plot_3panels(positions, y, lengths, figaxis3=axes3,
                                styles3=styles3, ytop=ytop, 
                                wholeblock=wholeblock,
                                show_line_at_1=show_line_at_1,
                                heights_instead=heights_instead,
                                label4legend=label4legend,
                                use_errorbar=show_errorbars,
                                yerror=yerror, unit=unit,
                                xerror=xerror, centered=centered)
            axes3[1].set_title(tit)                                
        else:
            fig, ax = styles.plot_3panels(positions, y, lengths,
                                          styles3=styles3, ytop=ytop, 
                                          wholeblock=wholeblock,
                                          show_line_at_1=show_line_at_1,
                                          label4legend=label4legend,
                                          heights_instead=heights_instead,
                                          use_errorbar=show_errorbars,
                                          yerror=yerror, unit=unit,
                                          xerror=xerror, centered=centered)
            ax[1].set_title(tit)
            fig.set_size_inches(6.5, 3.)
            fig.autofmt_xdate()

        # add pie chart showing % of total height or area
        if pie is True:
            if peak_idx is None:
                pass
            else:
                ax_pie = fig.add_subplot(339)
                atot, htot = prof.get_area_total()
                
                if heights_instead is False:
                    a = spec.peak_areas[peak_idx] / atot
                    size = [a, 1. - a]
                    tit = '% total area'
        
                else:
                    h = spec.peak_heights[peak_idx] / htot
                    size = [h, 1. - h]
                    tit = '% sum of heights'
            
                colors = ['k', 'w']
                plt.pie(size, colors=colors, 
                        startangle=90,
                        radius=0.25, center=(0, 0), frame=False)
                ax_pie.axis('equal')
                ax_pie.set_title(tit)
            
            if axes3 is None:
                return fig, ax

    def make_composite_peak(self, peak_idx_list):
        """Make composite peaks for all spectra in whole block"""
        for prof in self.profiles:
            for spec in prof.spectra:
                spec.make_composite_peak(peak_idx_list)
            prof.get_peak_info()

    def print_spectra_names(self, show_initials=True):
        """Print out fnames of all spectra associated with the whole-block 
        instance."""
        if show_initials is True:
            if self.initial_profiles is None:
                self.setupWB()
                
            if self.initial_profiles is not None:
                print('--Initial profiles--')
                for prof in self.initial_profiles:
                    print(prof.profile_name)
                    spec_list = []
                    for spectrum in prof.spectra:
                        spec_list.append(spectrum.fname)
                    print(spec_list)
                    print(' ')
            else:
                print('No initial profiles given')

        if self.profiles is not None:

            print('--Final profiles--')
            for prof in self.profiles:
                print(prof.profile_name)
                spec_list = []
                for spectrum in prof.spectra:
                    spec_list.append(spectrum.fname)
                print(spec_list)
                print(' ')

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
        """
        Make spectra baselines for all spectra whole-block. 
        
        Keywords are similar to spectrum.make_baseline
        """  
        for prof in self.profiles:
            prof.make_baselines(raw_data=raw_data, wn_low=wn_low, 
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

    def save_baselines(self, initial_too=True, 
                       baseline_ending='-baseline.CSV'):
        """Make and save spectra baselines for all spectra."""
        for prof in self.profiles:
            for spectrum in prof.spectra:
                spectrum.save_baseline(baseline_ending=baseline_ending)

    def plot_peakfits(self, initial_too=False, profile_idx=None, legloc=1):
        """Whole block: Plot peakfits for all spectra in all profiles"""
        if profile_idx is None:
            for prof in self.profiles:
                prof.plot_peakfits(initial_too, legloc=legloc)
        else: 
            self.profiles[profile_idx].plot_peakfits(initial_too, 
                                                    legloc=legloc)

    def print_max_arearatio(self, peak_idx=None, heights_instead=False):
        """ Prints out the maximum whole-block area ratio observed 
        in any profile of the wholeblock for specified peak_idx"""
        if peak_idx is None:
            self.setupWB(False, True)
        else:
            self.get_peakfit()
            self.setupWB(True, False)

        a = []

        for prof in self.profiles:
            if peak_idx is None:
                maxval = max(prof.wb_areas)
            else:
                prof.make_wholeblock(peakfit=True)
                if heights_instead is False:
                    maxval = max(prof.peak_wb_areas[peak_idx])
                else:
                    maxval = max(prof.peak_wb_heights[peak_idx])
            a.append(maxval)
        print('\n', self.name)
        print(max(a))
        
        return max(a)

    def print_peakfits_ave(self, printall=True, print_max=False, 
                           print_aves=True):
        """ Prints out and returns average peak areas, peak heights, and 
        sum of the average peak areas summed over all profiles"""
        areas = []
        heights = []
        totalareas = []
        max_a = []
        max_h = []
        for prof in self.profiles:
            a, h, totala, ma, mh = prof.print_peakfits_ave(printout=printall)
            areas.append(a)
            heights.append(h)
            totalareas.append(totala)
            max_a.append(ma)
            max_h.append(mh)

            if print_max is True:
                print(ma)
                print(mh)
            
        if print_aves is True:
            asum = np.array(np.sum(areas, axis=0))
            hsum = np.array(np.sum(heights, axis=0))
            tasum = np.sum(totalareas)
    
            print('\n', self.name)
            print('peak positions (cm-1)')
            print(self.profiles[0].spectra[0].peakpos)
            print('\naverage peak areas summed over all profiles (cm-2)')
            print(asum)
            print('\naverage peak heights summed over all profiles (cm-1)')
            print(hsum)
            print('\naverage total area summed over all profiles (cm-1)')
            print(tasum)
        return asum, hsum, tasum

    def print_diffusivities(self, peak_idx=None, profile_idx=None,
                            show_plot=False, top=1.5):
        """Print diffusivities for each profile"""
        if show_plot is True:
            self.plot_3panels_ave_spectra(peak_idx=peak_idx, top=top)
                     
        if peak_idx is None and profile_idx is None:
            for prof in self.profiles:
                prof.get_diffusivities()
                prof.print_diffusivities()
            return
        
        if profile_idx is None:
            for prof in self.profiles:
                k = peak_idx
                print('\n', prof.profile_name)
                print('peak position and log10(diffusivity in m2/s)')
                print('bulk H :', prof.D_area_wb, '+/-', \
                    prof.D_area_wb_error)
                print(prof.peakpos[k], ':', prof.D_peakarea_wb[k],\
                    '+/-', prof.peak_D_area_wb_error[k])
            return
            
        prof = self.profiles[profile_idx]
        print('\n', prof.profile_name)
        print('peak positions, log10(D in m2/s), D errors')
        print('bulk H :  ', prof.D_area_wb, '+/-', \
                prof.D_area_wb_error)
        print(prof.peakpos)
        print(prof.D_peakarea_wb)
        print(prof.peak_D_area_wb_error)

    def save_diffusivities(self, folder=None, 
                           file_ending='-diffusivities.txt'):
        """Save diffusivities for all profiles in whole-block instance
        to files"""
        for prof in self.profiles:
            prof.save_diffusivities(folder, file_ending)
            
    def get_diffusivities(self, folder=None, 
                           file_ending='-diffusivities.txt'):
        """Gets diffusivities for all profiles in whole-block instance
        from previously saved files"""
        if folder is None:
            folder = self.folder
        for prof in self.profiles:
            prof.get_diffusivities(folder, file_ending)

    def diffusion_profiles(self, wholeblock_data=False,
                           wholeblock_diffusion=False,
                           peak_idx=None, 
                           time_seconds=None, 
                           log10D_m2s=[-12., -12., -12.], 
                           erf_or_sum='erf', 
                           points=50, 
                           heights_instead=False, 
                           init=1., 
                           fin=0, 
                           approximation1D=False):
        """
        Returns lmfit parameters, x-data, and y-data for 3-dimensional 
        diffusion in a block.
        
        Requires time in seconds either explicitly passed here or as 
        attributes of the WholeBlock object.

        Assumes 3D non-path-integrated (wholeblock_diffusion=False), 
        but wholeblock_diffusion can be set to True.
        """
        if self.lengths is None:
            self.setupWB(peakfit=False, make_wb_areas=False)
        if self.lengths is None:
            print('Need to setup self.lengths, which is in microns')
            return False

        if self.directions is None:           
            self.setupWB(peakfit=False, make_wb_areas=False)

        if self.initial_profiles is None:
            self.setupWB(peakfit=False, make_wb_areas=False)

        if self.raypaths is None:
            self.setupWB(peakfit=False, make_wb_areas=False)

        if time_seconds is None:
            if self.time_seconds is not None:
                time_seconds = self.time_seconds
            else:
                print('Need time information')
                return False, False, False
                
        # Pick which diffusivities to use
        if log10D_m2s is not None:
            D3 = models.D_checker(log10D_m2s)
        elif wholeblock_diffusion is True and peak_idx is None:
            D3 = self.D_area_wb
        else:
            D3 = []
            for prof in self.profiles:
                D = prof.D_picker(wholeblock_data, heights_instead, peak_idx)
                D3.append(D)
        if D3 is None or 0.0 in D3:
            print('D3:', D3)
            print('\nNeed diffusivities.')
            print('Input directly as diffusivities_log10D_m2s')
            print('or input bulk in profile.D_area_wb')
            print('or peak_diffusivities at specified peak_idx\n')
            return False

        L3 = self.lengths
        params = models.params_setup3D(L3, D3, time_seconds, 
                                       init, fin)

        if wholeblock_diffusion is True:        
            xdiff, ydiff = models.diffusion3Dwb_params(params, 
                                                   raypaths=self.raypaths, 
                                                   erf_or_sum=erf_or_sum,
                                                   show_plot=False)           
        else:
            v, ydiff, xdiff = models.diffusion3Dnpi_params(params, 
                                                           points=points, 
                                                           centered=False)
        
        if wholeblock_data is False:
            maxareas = []
            for prof in self.profiles:
                if len(prof.fnames) > 0:
                    try:
                        maxa = np.max(prof.areas)
                    except AttributeError:
                        prof.make_areas()
                        maxa = np.max(prof.areas)
                    maxareas.append(maxa)
            ydiff = np.array(ydiff) * np.max(maxareas)

        return params, xdiff, list(ydiff)
            
    def plot_diffusion(self, wholeblock_data=False, 
                       wholeblock_diffusion=False,
                       peak_idx=None, 
                       time_seconds=None, 
                       log10D_m2s=[-12., -12., -12.], 
                       erf_or_sum='erf', 
                       show_plot=True, 
                       show_data=True,
                       xaxis='centered',
                       show_slice=False, 
                       label4legend=[None, None, None],
                       axes3=None, 
                       points=50, 
                       top_spectra=1.0,
                       ytop=None, 
                       numformat='{:.1f}', 
                       heights_instead=False, 
                       init=1., 
                       centered=True,
                       fin=0, 
                       approximation1D=False, 
                       labelD=True,
                       show_errorbars=True, 
                       labelDy=None, 
                       labelDx=[None, None, None],
                       style_data=styles.style_points,
                       style_diffusion=styles.style_1,
                       show_line_at_init=True,
                       ):
        """
        Applies 3-dimensionsal diffusion equations using equations in 
        pynams.diffusion.models and plots non-path-integrated 3-dimensional
        diffusion curves (default) or path-integrated whole-block profiles
        described in Ferriss et al. 2015 (set wholeblock_diffusion=True).
        
        If show_data is True (the default), also plots the data, either
        directly as the measured areas (default) or as the ratio of the
        measured area to a best-fit line through the initial areas
        (set wholeblock_data=True)
        
        See the help documentation for Block.diffusion_profiles()
        for more details about the diffusion modeling. 
        
        If axes3 = a list of 3 axes handles, the data and diffusion curve
        are plotted there. Otherwise, the figure handle and a list of 
        the 3 axes handles are returned.

        The diffusivities are automatically labeled on the plot. 
            To remove the labels, set labelD=False
            To change the number of digits, play with numformat keyword.
            To change the y position across all axes, use labelDy.
            To change the x positions, pass *a list of x-values* into labelDx.
        
        Change the maximum y value with the top keyword.
        """        
        D3 = models.D_checker(log10D_m2s)
        hinstead = heights_instead
        approx = approximation1D
        params, xdiff, ydiff = self.diffusion_profiles(
                               wholeblock_diffusion=wholeblock_diffusion,
                               peak_idx=peak_idx,
                               time_seconds=time_seconds,
                               log10D_m2s=D3,
                               erf_or_sum=erf_or_sum,
                               points=points,
                               heights_instead=hinstead,
                               init=init, fin=fin,
                               approximation1D=approx)

        if params is False:
            return False, False

        if axes3 is None:
            fig, axes3 = self.plot_areas_3panels(peak_idx=peak_idx, ytop=ytop,
                                              wholeblock=wholeblock_data, 
                                              heights_instead=heights_instead,
                                              show_line_at_1=False,
                                              label4legend=label4legend,
                                              styles3=[style_data]*3,
                                              centered=centered,
                                              show_errorbars=show_errorbars)
        if centered is True:
            for idx_len in range(3):
                xdiff[idx_len] = xdiff[idx_len] - (self.lengths[idx_len]/2.)

        if wholeblock_data is False:
            maxareas = []
            for prof in self.profiles:
                if len(prof.fnames) > 0:
                    try:
                        maxa = np.max(prof.areas)
                    except AttributeError:
                        prof.make_areas()
                        maxa = np.max(prof.areas)
                    maxareas.append(maxa)
            init = np.max(maxareas)

        styles.plot_3panels(xdiff, np.array(ydiff), 
                            show_line_at_1=show_line_at_init, 
                            figaxis3=axes3, init=init, ytop=ytop,
                            styles3=[style_diffusion]*3, 
                            label4legend=label4legend,
                            centered=centered)

        # label diffusivities
        if labelD is True:
            for k in range(3):
                if centered is True:
                    labelDx[k] = 0.
                else:
                    labelDx[k] = self.lengths[k]/2.
                dlabel = ''.join(('logD ', str(numformat.format(D3[k])), 
                                  ' m$^2$/s'))
                if labelDy is None:
                    ytop = axes3[0].get_ylim()[1]
                    labelDy = ytop - ytop*0.2
                axes3[k].text(labelDx[k], labelDy, dlabel, 
                              horizontalalignment='center',
                              verticalalignment='center',
                              backgroundcolor='w')
   
        try:
            return fig, axes3
        except NameError:
            pass
       
    def fitD(self, peak_idx=None, init=1., fin=0.,
             guesses_log10D=[-13., -13., -13.], 
             heights_instead=False, wholeblock_data=True,
             vary_initials=False, vary_finals=False, 
             vary_diffusivities=[True, True, True],
             erf_or_sum='erf', show_plot=True, wholeblock_diffusion=True, 
             centered=True, show_initial_guess=True, style_initial=None,
             style_final={'color' : 'red'}, points=50, top=1.2):
        """
        Forward modeling to determine diffusivities in three dimensions 
        from blocks of data.
        """        
        # x and y are the data that we will fit to, centered for fitting
        x, y = self.xy_picker(peak_idx, wholeblock_data, heights_instead, 
                              centered=True)
                
        # for processing results
        bestD = [] 
        D3 = []
        e3 = []

        # set up fitting parameters in the reqquired format
        params = models.params_setup3D(microns3=self.lengths, 
                 log10D3=guesses_log10D, 
                 time_seconds=self.time_seconds, 
                 initial=init, final=fin,
                 vinit=vary_initials, vfin=vary_finals,
                 vD=vary_diffusivities)

        # other keywords needed for forward model
        dict_fitting = {'points' : points,
                        'erf_or_sum' : erf_or_sum} 

        if wholeblock_diffusion is True:
            # need raypaths and don't plot twice
            if self.raypaths is None:
                self.setupWB()
            dict_fitting['raypaths'] = self.raypaths
            dict_fitting['show_plot'] = False
            
            # run the minimizer
            lmfit.minimize(models.diffusion3Dwb_params, 
                           params, args=(x, y), 
                           kws=dict_fitting)

            resid = models.diffusion3Dwb_params(params, x, y, 
                                            raypaths=self.raypaths,
                                            erf_or_sum=erf_or_sum,
                                            show_plot=False)
            
        elif wholeblock_diffusion is False:
            lmfit.minimize(models.diffusion3Dnpi_params, 
                           params, args=(x, y), 
                           kws=dict_fitting)
     
            resid = models.diffusion3Dnpi(params, x, y)
        else:
            print('wholeblock_diffusion must be either True or False')
            return            

        # convert to ufloats because ufloats are fun
        bestD.append(ufloat(params['log10Dx'].value, 
                            params['log10Dx'].stderr))
        bestD.append(ufloat(params['log10Dy'].value, 
                            params['log10Dy'].stderr))
        bestD.append(ufloat(params['log10Dz'].value, 
                            params['log10Dz'].stderr))
        bestinit = (ufloat(params['initial_unit_value'].value, 
                             params['initial_unit_value'].stderr))

        # Plot and print results
        for k in range(3):
            D3.append(bestD[k].n)
            e3.append(bestD[k].s)

        if show_plot is True:
            if wholeblock_diffusion is True:
                self.plot_diffusion(init=init, 
                                    peak_idx=peak_idx,
                                    log10D_m2s=D3,
                                    heights_instead=heights_instead,
                                    centered=centered)
            else:
                 print('sorry, only plotting wholeblock right now')
                                             
        print('\ntime in hours:', params['time_seconds'].value / 3600.)
        print('\ninitial unit values:', bestinit)
        print('\nbestfit log10D in m2/s:')
        for D in bestD:
            print(D)
        print('residual sum of squares:', np.sum(np.array(resid)**2.))
        print(D3[0], e3[0], D3[1], e3[1], D3[2], e3[2])
                             
        # Store values in profile attributes        
        for k in range(3):
            self.profiles[k].D_saver(D3[k], e3[k], wholeblock_data, 
                            heights_instead, peak_idx)
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
