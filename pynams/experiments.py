# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 09:33:13 2015

@author: Ferriss

supporting functions for performing experiments on 
nominally anhydrous minerals. For now there's only one function.

pressure_design() takes dimensions and details for sample in 
high pressure piston cylinder experiments
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

plt.style.use('paper') # my personal style sheet

style_pressure_medium = {'hatch' : 'O', 'facecolor' : 'lightgrey'}
style_graphite = {'hatch' : '/', 'facecolor' : 'dimgrey'}
style_MgO = {'hatch' : '..', 'facecolor' : 'white'}
style_pyrophyllite = {'hatch' : 'xx', 'facecolor' : 'hotpink'}
style_capsule = {'facecolor' : 'orange', 'edgecolor' : 'k'}
style_buffer = {'facecolor' : 'w', 'hatch' : '*', 'edgecolor' : 'g'}

def pressure_design(capsule_material = 'copper',
                    pressure_medium_material='BaCO$_3$',
                    sleeve_material='pyrophyllite',
                    buffer_material='H$_2$O, Ni,\nNiO, SC ol.,\nSC enst.',
                    lid_shape = 'bevel', # or flat or suaged
                    h_graphite_button=1.5, 
                    h_pressure_medium=33.35, # h for height
                    h_graphite_cylinder=33.35,
                    h_sleeve=10.5,
                    h_sleeve_bottom=0.65,
                    h_capsule = 8.65,
                    h_lid=2., # total thickness, even if it has a lip
                    h_MgO_base=10.6,
                    h_MgO_wafer=1.5,
                    h_MgO_top=10.6,
                    od_pressure_medium=18.9, # od for outer diameter
                    od_graphite_cylinder=11.5,
                    id_graphite_cylinder=10.0,
                    id_sleeve=8.7,
                    id_capsule=6.7,
                    legend_on=True, 
                    figsize=(3., 3),):
    """Creates and returns figure and axis
    showing experimental setup for a high pressure 
    experiment used to hydrate olivine and potentially
    other nominally anhydrous minerals in the green 4-post piston cylinder 
    at Lamont. All dimensions are input in mm. h=height, d=diameter,
    od=outer diameter, id=inner diameter
    """
    # reset all style labels to None for legend 
    style_capsule['color'] = None
    style_buffer['label'] = buffer_material
    for style in [style_graphite, style_MgO, style_capsule]:
        style['label'] = None
    
    d_graphite_button = od_pressure_medium
    od_MgO_base = id_graphite_cylinder
    od_sleeve = id_graphite_cylinder

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_xlim(0, od_pressure_medium)
    ax.set_xlabel('(mm)')
    ax.set_ylabel('(mm)')
    h_guts = h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top
    highest_point = max(h_pressure_medium, h_graphite_cylinder, h_guts)
    ax.set_ylim(0., h_graphite_button + highest_point + 2.)
    plt.tick_params(axis='x', top='off')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    th_gc = (od_graphite_cylinder - id_graphite_cylinder) / 2.
    xgc = (od_pressure_medium - od_graphite_cylinder) / 2.

    style_pressure_medium['label'] = pressure_medium_material
    pressure_medium = patches.Rectangle((0., h_graphite_button), # (x,y)
                                        od_pressure_medium, # width
                                        h_pressure_medium, # height
                                        **style_pressure_medium)

    graphite_button = patches.Rectangle((0., 0.), d_graphite_button,
                                        h_graphite_button, **style_graphite)

    style_graphite['label'] = 'graphite'
    graphite_cylinder = patches.Rectangle((xgc, h_graphite_button), 
                             od_graphite_cylinder, h_graphite_cylinder, 
                             **style_graphite)

    the_guts = patches.Rectangle((xgc + th_gc, h_graphite_button), 
                                 id_graphite_cylinder,
                                 h_graphite_cylinder, facecolor='w')
                             
    MgO_base = patches.Rectangle((xgc+th_gc, h_graphite_button),
                                 od_MgO_base, h_MgO_base, **style_MgO)

    def make_capsule_shape(x, y, height, outerD, innerD, shape='regular'):
        thick = (outerD - innerD) / 2.
        if shape == 'regular':
            verts = [(x + thick*2 + innerD, y),
                     (x, y),
                     (x, y + height),
                     (x + thick, y + height),
                     (x + thick, y + innerD/2.),
                     (x + thick + innerD/2., y + thick),
                     (x + thick + innerD, y + innerD/2.),
                     (x + thick + innerD, y + height),
                     (x + thick*2 + innerD, y + height),
                     (0., 0.)]
            codes = [Path.MOVETO] + ([Path.LINETO] * 8) + [Path.CLOSEPOLY]
        elif shape == 'suaged':
            th_flap = thick/2.
            verts = [(x + thick*2 + innerD, y),
                     (x, y),
                     (x, y + height + thick),
                     (x + thick, y + height),
                     (x + thick, y + height + outerD/2.),
                     (x + thick, y + th_flap + outerD/2.),
                     (x + thick, y + innerD/2.),
                     (x + thick + innerD/2., y + thick),
                     (x + thick + innerD, y + innerD/2.),
                     (x + thick + innerD, y + height),
                     (x + thick*2 + innerD, y + height),
                     (0., 0.)]
            codes = [Path.MOVETO] + ([Path.LINETO] * 10) + [Path.CLOSEPOLY]
        path = Path(verts, codes)
        return path
    
    # sleeve around capsule
    sleeve_path = make_capsule_shape(x=xgc + th_gc, 
                                     y=h_graphite_button + h_MgO_base,
                                     height=h_sleeve, 
                                     outerD=od_sleeve, 
                                     innerD=id_sleeve)
    if sleeve_material == 'MgO':
        style_sleeve = style_MgO.copy()
    elif sleeve_material == 'pyrophyllite':
        style_sleeve = style_pyrophyllite.copy()
        style_sleeve['label'] = 'pyrophyllite'
    else: 
        print 'unknown sleeve material. Assuming pyrophyllite'
        style_sleeve = style_pyrophyllite.copy()
        style_sleeve['label'] = 'pyrophyllite'

    sleeve = patches.PathPatch(sleeve_path, **style_sleeve)

    # capsule
    th_sleeve = (od_sleeve - id_sleeve) / 2.
    if (lid_shape == 'bevel') or (lid_shape == 'flat'):
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule)
    elif lid_shape == 'suaged':
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule,
                       shape='suaged')  
    else:
        print 'valid entries for lid_shape are flat, bevel, and suaged'
        capsule_path = make_capsule_shape(x=xgc + th_gc + th_sleeve, 
                       y=h_graphite_button + h_MgO_base + h_sleeve_bottom,
                       height=h_capsule, outerD=id_sleeve, innerD=id_capsule)                          
                       
    if capsule_material == 'copper':
        style_capsule['label'] = 'copper'
        style_capsule['facecolor'] = 'orange'
    elif capsule_material == 'silver':
        style_capsule['label'] = 'silver'
        style_capsule['facecolor'] = 'silver'
    elif capsule_material == 'gold':
        style_capsule['label'] = 'gold'
        style_capsule['facecolor'] = 'gold'
    elif capsule_material == 'platinum':
        style_capsule['label'] = 'platinum'
        style_capsule['facecolor'] = 'lightblue'
    else:
        print 'unknown capsule material'
        style_capsule['label'] = 'capsule'        
    capsule = patches.PathPatch(capsule_path, **style_capsule)

    # MgO on top
    MgO_wafer = patches.Rectangle((xgc + th_gc, 
                                   h_graphite_button + h_MgO_base + h_sleeve),
                                   od_MgO_base, h_MgO_wafer, **style_MgO)

    style_MgO['label'] = 'MgO'
    MgO_top = patches.Rectangle((xgc + th_gc, 
                    h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer),
                    od_MgO_base, h_MgO_top, **style_MgO)

    thermocouple = patches.Rectangle((od_pressure_medium/2.-0.5, 
                    h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer),
                    1., h_MgO_top, facecolor='w')
    ax.plot([od_pressure_medium/2-0.15, od_pressure_medium/2.-0.15],
            [h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top + 2.,
             h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer], 
             color='r', linewidth=1)
    ax.plot([od_pressure_medium/2+0.15, od_pressure_medium/2.+0.15],
            [h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer + h_MgO_top + 2.,
             h_graphite_button + h_MgO_base + h_sleeve + h_MgO_wafer], 
             color='b', linewidth=1)

    # buffer
    th_capsule = (id_sleeve - id_capsule) / 2.
    buffer_inside = patches.Rectangle((xgc + th_gc + th_sleeve + th_capsule,
                    h_graphite_button + h_MgO_base + h_sleeve_bottom + th_capsule),
                    id_capsule, h_capsule - th_capsule, **style_buffer)

    # capsule lid
    del style_capsule['label'] # so it doesn't appear twice in the legend        
    if lid_shape == 'flat':
        lid = patches.Rectangle((xgc + th_gc + th_sleeve, 
                    h_graphite_button + h_MgO_base + h_sleeve_bottom + h_capsule),
                    id_sleeve, h_lid, **style_capsule)
    elif lid_shape == 'bevel':
        x = xgc + th_gc + th_sleeve
        y = h_graphite_button + h_MgO_base + h_sleeve_bottom + h_capsule
        th_lid = h_lid / 2.
        th_capsule = (id_sleeve - id_capsule) / 2.
        lid_verts = [(x + th_capsule, y), 
                     (x, y),
                     (x, y + th_lid), 
                     (x + id_sleeve, y + th_lid), 
                     (x + id_sleeve, y), 
                     (x + id_sleeve - th_capsule, y), 
                     (x + id_sleeve - th_capsule, y - th_lid), 
                     (x + th_capsule, y - th_lid), 
                     (0., 0.)]
        lid_codes = [Path.MOVETO] + ([Path.LINETO] * 7) + [Path.CLOSEPOLY]
        lid_path = Path(lid_verts, lid_codes)
        lid = patches.PathPatch(lid_path, **style_capsule)
    elif lid_shape == 'suaged':
        th_flap = th_capsule / 2.
        x = xgc + th_gc + th_sleeve + th_flap
        y = h_graphite_button + h_MgO_base + h_sleeve_bottom + h_capsule - th_flap
        th_lid = h_lid / 2.
        th_capsule = (id_sleeve - id_capsule) / 2.
        lid_verts = [(x + th_capsule - th_flap, y), 
                     (x, y),
                     (x, y + th_lid), 
                     (x + id_sleeve - th_capsule, y + th_lid), 
                     (x + id_sleeve - th_capsule, y), 
                     (x + id_sleeve - th_capsule - th_flap, y), 
                     (x + id_sleeve - th_capsule - th_flap, y - th_lid), 
                     (x + th_capsule - th_flap, y - th_lid), 
                     (0., 0.)]
        lid_codes = [Path.MOVETO] + ([Path.LINETO] * 7) + [Path.CLOSEPOLY]
        lid_path = Path(lid_verts, lid_codes)
        lid = patches.PathPatch(lid_path, **style_capsule)
    else:
        print 'valid entries for lid_shape are flat, bevel, and suaged'
        lid = patches.Rectangle((xgc + th_gc + th_sleeve, 
                    h_graphite_button + h_MgO_base + th_sleeve + h_capsule),
                    id_sleeve, h_lid, **style_capsule)


    ax.add_patch(pressure_medium)
    ax.add_patch(graphite_button)
    ax.add_patch(graphite_cylinder)
    ax.add_patch(the_guts)
    ax.add_patch(MgO_base)
    ax.add_patch(sleeve)
    ax.add_patch(buffer_inside)
    ax.add_patch(MgO_wafer)
    ax.add_patch(MgO_top)
    ax.add_patch(thermocouple)
    ax.add_patch(lid)
    ax.add_patch(capsule)
      
    fig.tight_layout()
    if legend_on is True:
        plt.subplots_adjust(right=0.55, left=0.17, bottom=0.15, top=0.9)
        ax.legend(bbox_to_anchor=(2.25, 0.8), frameon=False)
    else:
        plt.subplots_adjust(right=0.9, left=0.17, bottom=0.15, top=0.9)
    return fig, ax