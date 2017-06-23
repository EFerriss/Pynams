# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:42:23 2017

@author: Elizabeth Ferriss
"""
from bokeh.plotting import figure, output_file, show, ColumnDataSource, save
from bokeh.layouts import layout, widgetbox
from bokeh.models import HoverTool, Range1d, BoxZoomTool, SaveTool, ResetTool
from bokeh.models import PanTool
from bokeh.io import curdoc, output_notebook
import pandas as pd
import numpy as np
import pynams
from bokeh.models.widgets import RadioButtonGroup, RangeSlider
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.embed import file_html

olivine_path = ''.join((pynams.thisfolder, 'diffusion\\literaturevalues.csv'))
olivine = pd.read_csv(olivine_path)
olivine.fillna(0, inplace=True) # replace missing values with zero
olivine.loc[olivine['orientation'] == 'u', 'orientation'] = 'not oriented'
olivine.loc[olivine['name'] == 0, 'name'] = ''
olivine.loc[olivine['fO2Buffer'] == 0, 'fO2Buffer'] = ''

olivine["color"] = np.where(olivine["Author"] == 'Ferriss', "orange", "grey")
olivine["alpha"] = np.where(olivine["Author"] == 'Ferriss', 0.9, 0.25)

hover = HoverTool(tooltips=[
        ("||", "@orient"),
        ("Fe", "Fo# @fonum{0.0} @buffer"),
        ("Heat", "@hours hrs at @celsius C "),
        ("Type", "@exper, @mech"),
        ("Source", "@author @year @name"),
])

source = ColumnDataSource(data=dict(
        x = [],
        y = [],
        color = [],
        year = [],
        author = [],
        fonum = [],
        celsius = [],
        orient = [],
        exper = [],
        mech = [],
        alpha = [],
        hours = [],
        buffer = [],
        name = []
    ))

left, right, bottom, top = 6, 10, -16, -8
p = figure(plot_width=500, plot_height=500, 
           tools=[BoxZoomTool(), hover, PanTool(), SaveTool(), ResetTool()],
           title="H diffusion in olivine", 
           x_range=Range1d(left, right), y_range=Range1d(bottom, top))
p.circle('x', 'y', size=10, source=source, color='color', alpha='alpha')
p.xaxis.axis_label = "1e4 / Temperature (K)"
p.yaxis.axis_label = "log10 Diffusivity (m2/s)"

# widgets
widget_orient = RadioButtonGroup(
        labels=['|| a', '|| b', '|| c', 'not oriented', 'all'], active=4)
widget_mech = RadioButtonGroup(
        labels=['bulk H (pp)', 'bulk H (pv)', '[Si]', '[Ti]', '[tri]', '[Mg]', 'all'], active=6)
widget_fo = RangeSlider(title='Fo#', start=80, end=100, range=(80, 100))

def select_data():
    selected = olivine
    orient_val = widget_orient.active
    mech_val = widget_mech.active
    fomax_val = widget_fo.range[1]
    fomin_val = widget_fo.range[0]
    
    selected = selected[selected.Fo.values <= fomax_val]
    selected = selected[selected.Fo.values >= fomin_val]
    
    orient_labels = ['a', 'b', 'c', 'not']
    if (orient_val != 4):
        label = orient_labels[orient_val]
        selected = selected[selected.orientation.str.contains(label)==True]

    mechs = ['pp', 'pv', 'Si', 'Ti', 'tri', 'Mg']
    if (mech_val != 6):
        mech = mechs[mech_val]
        selected = selected[selected.mechanism.str.contains(mech)==True]
        
    return selected


def update():
    df = select_data()
    
    source.data=dict(
        x = 1e4 / (df['celsius'].values + 273.15),
        y = df['log10D'].values,
        color = df['color'],
        alpha = df['alpha'],
        year = df['Year'].values,
        author = df['Author'].values,
        fonum = df['Fo'].values,
        celsius = df['celsius'].values,
        orient = df['orientation'].values,
        exper = df['Experiment'].values,
        mech = df['mechanism'].values,
        hours = df['hours'].values,
        buffer = df['fO2Buffer'].values,
        name = df['name'].values
    )

update()

# set up callbacks
controls = [widget_orient, widget_mech, widget_fo]
for control in [widget_orient, widget_mech]:
    control.on_change('active', lambda attr, old, new: update())
widget_fo.on_change('range', lambda attr, old, new: update())

# layout
sizing_mode = 'fixed' 
inputs = widgetbox(*controls, sizing_mode=sizing_mode)
layout = layout([
            [inputs, p]
            ], sizing_mode=sizing_mode)

curdoc().add_root(layout)
curdoc().title = "Arrhenius Diagram"

output_file('ArrheniusDiagram.html', title='Bokeh Plot', 
            mode='cdn', root_dir=None)
save(layout)
#html = file_html(p, resources=(CDN('inline'), CDN),
#                 title="Arrhenius") # CDN,
                 
#html_file = olivine_path = ''.join((pynams.thisfolder, 
#                                    'diffusion\\ArrheniusDiagram.html'))
#with open('Arrhenius.html', 'w') as file:
#   file.write(html)
# %% Mf additions
#with open("Mark.html", 'w') as mf:
#    mf.write(html)
#

#    
#layout = row(
#    p,
#    widgetbox(amp_slider, freq_slider, phase_slider, offset_slider),
#)
#output_file("mf4.html", title="slider.py example")
#
#show(l)