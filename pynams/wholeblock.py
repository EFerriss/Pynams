"""
Code for grouping 3 orthogonal Profiles into a single WholeBlock object.
"""
from __future__ import print_function, division, absolute_import
import pynams.styles as styles
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import xlsxwriter


class WholeBlock():
    def __init__(self,profiles=[], name='', peakfit=False, folder='',
                 make_wb_areas=False, time_seconds=None, worksheetname=None,
                 style_base = None, temperature_celsius=None,
                 diffusivities_log10_m2s=None, get_baselines=False,
                 diffusivity_errors=None, sample=None, D_area_wb=[-12., -12., -12.]):
        self.profiles = profiles
        self.folder = folder
        self.name = name
        self.time_seconds = time_seconds
        self.worksheetname = worksheetname
        self.style_base = style_base
        self.temperature_celsius = temperature_celsius
        self.sample = sample        
        self.D_area_wb = D_area_wb
        self.D_area_wb_error = [0., 0., 0.]
        self.peak_diffusivities = []
        self.peak_diffusivities_errors = []
        
        if len(self.profiles) > 0:
            self.setupWB(peakfit=peakfit, make_wb_areas=make_wb_areas,
                         get_baselines=get_baselines)
                
    def setupWB(self, peakfit=False, make_wb_areas=False, get_baselines=False):
        """Sets up and checks WholeBlock instance
        - Check that profiles list contains a list of three (3) profiles
        - Generate list of initial profiles
        - Generate list of profile directions
        - Verify three profile directions are orthogonal (['a', 'b', 'c'])
        - Generate list of ray paths
        - Verify three ray path directions are compatible with directions list
        """
        if len(self.profiles) != 3:
            print('For now, only a list of 3 profiles is allowed')
            return False
        d = []
        r = []
        ip = []
        L = []
        
        for prof in self.profiles:
#            if isinstance(prof, Profile) is False:
#                print 'Only profile objects allowed in profile list!'
#                return False
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
#                print 'Need to set initial_profile attribute for each profile.'
#                print 'Assuming these are initial profiles...'
                prof.initial_profile = prof
            ip.append(prof.initial_profile)
            
            if make_wb_areas is True:
                check = prof.make_wholeblock(peakfit=peakfit)
                if check is False:
                    return False

            if prof.spectra is None:
                prof.make_spectra()
                
            for spec in prof.spectra:
                if spec.thickness_microns is None and self.sample is not None:
                    if prof.raypath == 'a':
                        spec.thickness_microns = np.mean(self.sample.length_a_microns)
                    elif prof.raypath == 'b':
                        spec.thickness_microns = np.mean(self.sample.length_b_microns)
                    elif prof.raypath == 'c':
                        spec.thickness_microns = np.mean(self.sample.length_c_microns)
                    else:
                        print('Need wb.direction to assign thickness')

        self.directions = d
        self.raypaths = r
        self.initial_profiles = ip 
        self.lengths = L

        if peakfit is True:
            for prof in self.profiles:
                prof.get_peakfit()

        if get_baselines is True:        
            self.get_baselines()

        return True

    def get_peakfit(self, peak_ending='-peakfit.CSV', 
                    baseline_ending='-baseline.CSV'):
        """Get peakfit information for all profiles"""
        for prof in self.profiles:
            prof.get_peakfit(peak_ending=peak_ending,
                             baseline_ending=baseline_ending)

    def print_baseline_limits(self, initial_too=True):
        """Print out baseline wavenumber range for each profile"""
        for prof in self.profiles:
            print('\n', prof.profile_name)
            print(prof.spectra[0].base_low_wn, end=' ') 
            print(prof.spectra[0].base_high_wn)
            if initial_too is True:
                print(prof.initial_profile.profile_name)
                print(prof.initial_profile.spectra[0].base_low_wn, end=' ') 
                print(prof.initial_profile.spectra[0].base_high_wn)

    def get_baselines(self, initial_too=False, folder=None, delim=',', 
                      baseline_ending='-baseline.CSV'):
        """Get baselines for all spectra in whole block"""
        if self.initial_profiles is None:
            self.setupWB()
        for prof in self.profiles:
            for spectrum in prof.spectra:
                spectrum.get_baseline(baseline_ending=baseline_ending,
                                      folder=folder, delim=delim)
        if initial_too is True:
            for prof in self.initial_profiles:
                for spectrum in prof.spectra:
                    spectrum.get_baseline(baseline_ending=baseline_ending,
                                          folder=folder, delim=delim)
            
    def plot_showbaselines(self):
        """Plot baselines for all spectra in the whole block"""
        for prof in self.profiles:
            for spec in prof.spectra:
                spec.plot_showbaseline()


    def make_area_lists(self, polyorder=1, show_plot=False, set_class=None,
                       shiftline=None, printout_area=False, peak=None):
            """Make list of areas from all profiles"""
            self.areas = []
            for prof in self.profiles:
                a = prof.make_area_list(polyorder, show_plot, set_class,
                                        shiftline, printout_area, peak)
                self.areas.append(a)
           

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
        """Three suplots showing average initial and final spectra in each
        direction"""
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
                    if prof.wb_areas is None:
                        print('\nMaking whole block area ratios')
                        check = prof.make_wholeblock(peakfit=False, bulk=True,
                                                     show_plot=False)
                        if check is False:
                            return False        
                    y_to_add = prof.wb_areas
                
                # absolute areas
                else:
                    if prof.areas_list is None:
                        prof.make_area_list()
                    y_to_add = prof.areas_list


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

    def plot_areas_3panels(self, peak_idx=None, fig_ax3=None, centered=True,
                           top=None, wn=None, figsize=(6.5, 2.5), 
                           show_spectra=True, percent_error=0., 
                           xerror=0., yerror=None, pie=True,
                           label4legend=[None, None, None],
                           styles3=[styles.style_points]*3,
                           use_area_profile_styles=True, unit='microns',
                           heights_instead=False, wholeblock=True,
                           show_line_at_1=True, get_saved_baseline=True,
                           show_errorbars=True, peak_group=None):
        """Plot whole-block ratio of Area/Initial Area (default) 
        OR just areas (set wholeblock=False) on three panels"""
#        if get_saved_baseline is True:
#            self.get_baselines()
#
        if wholeblock is True:
            if ((self.directions is None) or (self.raypaths is None) or
                (self.initial_profiles is None) or (self.lengths is None)):
                    check = self.setupWB(make_wb_areas=False, peakfit=False)
                    if check is False:
                        print('Problem setting up whole block in setupWB')                    
                        return False
        
        if peak_idx is not None:
            for prof in self.profiles:
                for spec in prof.spectra:
                    if spec.peak_areas is None:
                        spec.get_peakareas()    

        peakpos = self.profiles[0].peakpos

        # concatenate positions and areas across three profiles to 
        # send in to the plotting function
        if peak_group is not None:
            positions, y_placeholder = self.xy_picker(peak_idx=peak_idx, 
                                                      wholeblock=wholeblock,
                                                      heights_instead=heights_instead, 
                                                      centered=centered, 
                                                      unit=unit)
            y = np.zeros_like(y_placeholder)
            tit = 'Sum of peaks'

            for peak_group_idx in peak_group:
                positions, y_add = self.xy_picker(peak_idx=peak_group_idx, 
                                                  wholeblock=wholeblock,
                                                  heights_instead=heights_instead, 
                                                  centered=centered, unit=unit)
                y = y + y_add
                tit = ' '.join((tit, str(peak_group_idx)))
            
        else:
            
            positions, y = self.xy_picker(peak_idx=peak_idx, wholeblock=wholeblock,
                                          heights_instead=heights_instead, 
                                          centered=centered, unit=unit)
            
            # Change title if peak-specific rather than bulk
            if peak_idx is not None:
                tit = ' '.join(('Peak at', str(peakpos[peak_idx]), '/cm'))
            else:
                tit = 'Bulk hydrogen'

        if top is None:
            z = [max(y[0]), max(y[1]), max(y[2])]
            top = max(z) + 0.1*max(z)

            
        if use_area_profile_styles is True and None in styles3:
            for k in range(3):
                styles3[k] = self.profiles[k].choose_marker_style()
                styles3[k]['markersize'] = 10
        
        if unit == 'microns':
            lengths = self.lengths
        elif unit == 'mm':
            lengths = np.array(self.lengths) / 1000.
        else:
            print('unit must be microns (default) or mm')
            return
                
        # Sent positions and areas to plotting command
        if fig_ax3 is not None:
            styles.plot_3panels(positions, y, lengths, figaxis3=fig_ax3,
                                styles3=styles3, top=top, wholeblock=wholeblock,
                                show_line_at_1=show_line_at_1,
                                heights_instead=heights_instead,
                                label4legend=label4legend,
                                use_errorbar=show_errorbars,
                                yerror=yerror, unit=unit,
                                percent_error=percent_error,
                                xerror=xerror, centered=centered)
            fig_ax3[1].set_title(tit)                                
        else:
            fig, ax = styles.plot_3panels(positions, y, lengths,
                                          styles3=styles3, top=top, 
                                          wholeblock=wholeblock,
                                          show_line_at_1=show_line_at_1,
                                          label4legend=label4legend,
                                          heights_instead=heights_instead,
                                          use_errorbar=show_errorbars,
                                          percent_error=percent_error,
                                          yerror=yerror, unit=unit,
                                          xerror=xerror, centered=centered)
            ax[1].set_title(tit)
            fig.set_size_inches(figsize)
            fig.autofmt_xdate()

            if pie is True:
                # add pie chart showing % of total height or area
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

    def make_profile_list(self, initial_too=False):
        """Return False or a list of profiles"""
        if initial_too is True:
            if self.initial_profiles is None:
                self.setupWB()
            if self.initial_profiles is None:
                print('No initial profiles')
                profile_list = self.profiles
            else:
                profile_list = self.initial_profiles + self.profiles
        else:
            profile_list = self.profiles
        return profile_list

    def make_baselines(self, initial_too=False, linetype='line', 
                       wn_high=3700., wn_low=3200., shiftline=None, 
                       show_fit_values=False, show_plot=False,
                       wn_mid=None): 
        """Make spectra baselines for all spectra in all profiles in 
        whole-block"""        
        profile_list = self.make_profile_list(initial_too)        
        for prof in profile_list:
            prof.make_baselines(linetype=linetype, shiftline=shiftline, 
                                show_fit_values=show_fit_values, 
                                show_plot=show_plot, wn_mid=wn_mid,
                                wn_high=wn_high, wn_low=wn_low) 

    def save_baselines(self, initial_too=True):
        """Make and save spectra baselines for all spectra."""
        profile_list = self.make_profile_list(initial_too)
        for prof in profile_list:
            for spectrum in prof.spectra:
                spectrum.save_baseline()

    def matlab(self, initial_too=False):
        """Print out a list of spectra names in a matlab-friendly way"""        
        print('\nFor use in FTIR_peakfit_loop.m\n')        
        
        if initial_too is True:
            if self.initial_profiles is None:
                self.setupWB()
            string = "{"
            for prof in self.initial_profiles:
                for spec in prof.spectra:
                    stringname = spec.fname
                    string = string + "'" + stringname + "' "
            string = string + "};"
            print(string, '\n')
        
        string = "{"
        for prof in self.profiles:
            for spec in prof.spectra:
                stringname = spec.fname
                string = string + "'" + stringname + "' "
        string = string + "}"
        print(string)

    def plot_areas(self, profile_index=None, peak_idx=None, peakwn=None, 
                   show_initials=False, show_finals=True, show_legend=True,
                   legloc=1, frame=False, top=None, bestfitlines=False,
                   heights_instead=False, together=False):
        """Plots profiles on one figure.
        Set initial_instead_of_final to True to see initial areas.
        Need to add legend and checks is wavenumber not in list."""

        # Which profiles to plot
        if show_initials is True:
            if self.initial_profiles is None:
                self.setupWB()
            if self.initial_profiles is None:
                print('Need initial profile')
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
            bstyleline[0] = dict(list(bstylelinebase.items()) + list(styles.style_Dx.items()))
            bstyleline[1] = dict(list(bstylelinebase.items()) + list(styles.style_Rx.items()))
            bstyleline[2] = dict(list(bstylelinebase.items()) + list(styles.style_Dy.items()))
            bstyleline[3] = dict(list(bstylelinebase.items()) + list(styles.style_Ry.items()))
            bstyleline[4] = dict(list(bstylelinebase.items()) + list(styles.style_Dz.items()))
            bstyleline[5] = dict(list(bstylelinebase.items()) + list(styles.style_Rz.items())) 
            
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

        # move y-axis upper limit to accomodate all data
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
                    
#        if peak_idx is None:
#            tit = 'Bulk hydrogen'
#        else:
#            tit = ' '.join(('Peak at',
#                            str(self.profiles[idx].peakpos[peak_idx]),'/cm'))
#        ax.text(0, top + 0.05*top, tit, horizontalalignment='center')
                   
        f.set_size_inches((6.5, 3.5))
        plt.subplots_adjust(top=0.9, bottom=0.35)
                   
        return f, ax

    def plot_peakfits(self, initial_too=False, profile_idx=None, legloc=1):
        """Whole block: Plot peakfits for all spectra in all profiles"""
        if profile_idx is None:
            for prof in self.profiles:
                prof.plot_peakfits(initial_too, legloc=legloc)
        else: 
            self.profiles[profile_idx].plot_peakfits(initial_too, 
                                                    legloc=legloc)


    def excelify(self, filename=None, exceldpi=150, top=1.0, 
                 bestfitlines_on_areas=True, include_peakfits=True):
        """Save all peak fit info AND peak figures to an excel spreadsheet """
        # file setup
        if filename is None and self.worksheetname is None:
            print('Need filename information here')
            print('or in worksheetname attribute of whole block instance')
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
        col4 = 28
        col_heights = 20
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

        peakpos = self.profiles[0].spectra[0].peakpos
        
        worksheet.write(0, 0, self.name)
        worksheet.write(1, 0, 'peak position (/cm)', wraptext)
        worksheet.write(1, 1, 'peak height (/cm)', wraptext)
        worksheet.write(1, 2, 'peak width (/cm)', wraptext)
        worksheet.write(1, 3, 'peak area (/cm2)', wraptext)
        worksheet.write(1, 5, 'Baselines', wraptext)        
        worksheet.write(1, 9, 'Baseline-subtracted', wraptext)
        worksheet.set_column(9, 9, width=11)
        
        if include_peakfits is True:
            # Bulk area and peak fits in 1st column for individual fits
            row = 2        
            pic_counter = 0
            for prof in self.profiles:
                row = row + 1
                worksheet.write(row, col1, prof.profile_name, boldtext)
                row = row + 1
    
                pos_idx = 0
                for spec in prof.spectra:                
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
                    f.savefig(newfig, dpi=exceldpi, format='png')
                    worksheet.insert_image(row, col2, newfig, 
                                           {'x_scale' : 0.5, 'y_scale' : 0.5})
    
                    # Peak fit figures in 2nd column
                    f, ax = spec.plot_peakfit()
                    newfig = 'peakfit{:d}.png'.format(pic_counter)
                    pic_counter = pic_counter + 1
                    f.savefig(newfig, dpi=exceldpi, format='png')
                    worksheet.insert_image(row, col2pt5, newfig, 
                                           {'x_scale' : 0.5, 'y_scale' : 0.5})
    
                    pos_idx = pos_idx + 1
                    row = row + 1
    
                    if spec.peakpos is None:
                        check = spec.get_peakfit()
                        if check is False:
                            print('trouble with getting peakfit')
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
        f.savefig('panel3.png', dpi=exceldpi, format='png')
        worksheet.insert_image(0, col3, 'panel3.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
            
        # area profiles, height profiles, whole-block profiles for both
        # heights are for peak-specific only
        worksheet.write(11, col3, 'Peak areas profiles', boldtext)
        worksheet.write(11, col4, 'Whole-block area profiles', boldtext)
        worksheet.write(25, col_heights, 'Peak height profiles', boldtext)
        worksheet.write(25, col_wb_heights, 'Whole-block height profiles', 
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
        worksheet.write(26, col_heights, errorstring3)
        worksheet.write(26, col_wb_heights, errorstring4)

        worksheet.write(13, col3, 'Errors typically plot in or near symbols')
        worksheet.write(13, col4, 'Errors typically plot in or near symbols')
        worksheet.write(27, col_heights, 
                                'Errors typically plot in or near symbols')
        worksheet.write(27, col_wb_heights, 
                                'Errors typically plot in or near symbols')
        
        f, ax = self.plot_areas_3panels(wholeblock=False)
        f.savefig('bulk_areas.png', dpi=exceldpi, format='png')
        worksheet.insert_image(15, col3, 'bulk_areas.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
        
        f, ax = self.plot_areas_3panels(wholeblock=True)
        f.savefig('wbbulk.png', dpi=exceldpi, format='png')
        worksheet.insert_image(15, col4, 'wbbulk.png', {'x_scale' : 0.5, 
                               'y_scale' : 0.5})
                       
        for peak_idx in range(len(self.profiles[0].spectra[0].peakpos)):
            f, ax = self.plot_areas_3panels(wholeblock=False, 
#                                            bestfitlines=bestfitlines_on_areas, 
                                            peak_idx=peak_idx)
            newfig = 'area{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi, format='png')
            worksheet.insert_image(29+(14*peak_idx), col3, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            

            f, ax = self.plot_areas_3panels(peak_idx=peak_idx, wholeblock=True)
            newfig = 'wbareas{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi, format='png')
            worksheet.insert_image(29+(14*peak_idx), col4, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})
                   
            f, ax = self.plot_areas_3panels(wholeblock=False,
#                                            bestfitlines=bestfitlines_on_areas, 
                                            peak_idx=peak_idx, 
                                            heights_instead=True)
            newfig = 'height{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi, format='png')
            worksheet.insert_image(29+(14*peak_idx), col_heights, newfig, 
                                   {'x_scale' : 0.5, 'y_scale' : 0.5})            

            f, ax = self.plot_areas_3panels(peak_idx=peak_idx, 
                                            heights_instead=True,
                                            wholeblock=True)
            newfig = 'wbheights{:d}.png'.format(peak_idx)
            f.savefig(newfig, dpi=exceldpi, format='png')
            worksheet.insert_image(29+(14*peak_idx), col_wb_heights, newfig, 
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
            spec = prof.spectra[0]
            ispec = iprof.spectra[0]
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
        for peak in peakpos:
            label1 = ''.join((str(peak), ' height (/cm)'))
            worksheet.write(1, col, label1, wraptext)
            col = col + 1
        for peak in peakpos:
            label2 = ''.join((str(peak), ' whole-block height ratio'))
            worksheet.write(1, col, label2, wraptext)
            col = col + 1
        worksheet.write(1, col, 'file label', wraptext)

        # values filling in the rows
        row = 2
        for prof in self.profiles:
            col = col5
  
            # label profile name with spaces on either side
            worksheet.write(row, col, prof.profile_name)            
            halflen = prof.set_len() / 2.            
            pos_idx = 0
            row = row + 1

            for spec in prof.spectra:
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

                # peak heights
                for k in range(len(peakpos)):
                    height = prof.peak_heights[k][pos_idx]
                    worksheet.write(row, col, height, 
                                    two_digits_after_decimal)
                    col = col + 1
                # peak whole-block heights
                for k in range(len(peakpos)):
                    wbh = prof.peak_wb_heights[k][pos_idx]
                    if np.isnan(wbh):  
                        worksheet.write(row, col, 'nan')
                    elif np.isinf(wbh): 
                        worksheet.write(row, col, 'inf')
                    else:
                        worksheet.write(row, col, wbh, 
                                    two_digits_after_decimal)
                    col = col + 1

                # filenames at the end of summary
                worksheet.write(row, col, spec.fname)

                row = row + 1
                pos_idx = pos_idx + 1

        workbook.close()        

    def print_peakfits(self, initial_too=False, excelfriendly=True):
        """Print out all peakfit information for each spectrum in 
        each profile"""
        if initial_too is True and self.initial_profiles is not None:
            proflist = self.profiles + self.initial_profiles
        else:
            proflist = self.profiles

        if excelfriendly is True:
            print('position height width area')
            
        for prof in proflist:
            prof.print_peakfits()

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

#            if np.isnan(maxval) is False:
            a.append(maxval)
#            else:
#                a.append(0)
        
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
#                print '\n', prof.profile_name
#                print 'max areas'
                print(ma)
#                print 'max_heights'
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


#        return asum, hsum, tasum, masum, mhsum


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
#            print('Need to setup self.lengths, which is in microns')
#            return False
#
#        if (wb_or_3Dnpi != 'npi') and (wb_or_3Dnpi != 'wb'):
#            print('wb_or_3Dnpi only takes "wb" or "npi"')
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
#                print('Need time information')
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
#            print('D3:', D3)
#            print('\nNeed diffusivities.')
#            print('Input directly as diffusivities_log10D_m2s')
#            print('or input bulk in profile.D_area_wb')
#            print('or peak_diffusivities at specified peak_idx\n')
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
#                for idx_len in range(3):
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
#                for k in range(3):
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
#            print('npi not working well right now, sorry')
#            lmfit.minimize(diffusion.diffusion3Dnpi_params, 
#                           params, args=(x, y), 
#                           kws=dict_fitting)
#     
#            resid = diffusion.diffusion3Dnpi(params, x, y)
#        else:
#            print('wb_or_3Dnpi can only be wb or npi')
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
#                 print('sorry, only plotting wholeblock right now')
#                                             
#        print('\ntime in hours:', params['time_seconds'].value / 3600.)
#        print('\ninitial unit values:', bestinit)
#        print('\nbestfit log10D in m2/s:')
#        for D in bestD:
#            print(D)
#        print('residual sum of squares:', np.sum(np.array(resid)**2.))
#        print(D3[0], e3[0], D3[1], e3[1], D3[2], e3[2])
#                             
#        # Store values in profile attributes        
#        for k in range(3):
#            self.profiles[k].D_saver(D3[k], e3[k], wholeblock, 
#                            heights_instead, peak_idx)
#        return bestD
    
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
