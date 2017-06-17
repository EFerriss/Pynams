# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:42:23 2017

@author: Elizabeth Ferriss
"""

from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.layouts import layout, widgetbox
from bokeh.models import HoverTool
from bokeh.io import curdoc
import pandas as pd
import numpy as np
import pynams
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider, TextInput
from bokeh.client import push_session
from bokeh.driving import cosine
from bokeh.plotting import figure, curdoc

data_path = ''.join((pynams.thisfolder, 'diffusion\\literaturevalues.csv'))
data = pd.read_csv(data_path)
data.fillna(0, inplace=True) # replace missing values with zero
data.loc[data['orientation'] == 'u', 'orientation'] = 'not oriented'

data["color"] = np.where(data["Author"] == 'Ferriss', "orange", "grey")
data["alpha"] = np.where(data["Author"] == 'Ferriss', 0.9, 0.25)
#%%
hover = HoverTool(tooltips=[
        ("||", "@orient"),
        ("Fo#", "@fonum"),
        ("Temp.", "@celsius C"),
        ("Type", "@exper"),
        ("Source", "@author @year"),
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
        exper = []
    ))

p = figure(plot_width=400, plot_height=400, tools=[hover],
           title="H diffusion in olivine")
p.circle('x', 'y', size=10, source=source, color='color')
p.xaxis.axis_label = "1e4 / Temperature (K)"
p.yaxis.axis_label = "log10 Diffusivity (m2/s)"
    
# widgets
widget_orient = RadioButtonGroup(
        labels=['|| a', '|| b', '|| c', 'not oriented', 'all'], active=4)

def select_data():
    selected = data
    orient_val = widget_orient.active
    if (orient_val != 4):
        if (orient_val == 0):  
            selected = selected[selected.orientation.str.contains('a')==True]
        elif (orient_val == 1):
            selected = selected[selected.orientation.str.contains('b')==True]
        elif (orient_val == 2):    
            selected = selected[selected.orientation.str.contains('c')==True]
        elif (orient_val == 3):    
            selected = selected[selected.orientation.str.contains('not oriented')==True]
    return selected


def update():
    df = select_data()
    
    source.data=dict(
        x = 1e4 / (df['celsius'].values + 273.15),
        y = df['log10D'].values,
        color = df['color'],
        year = df['Year'].values,
        author = df['Author'].values,
        fonum = df['Fo'].values,
        celsius = df['celsius'].values,
        orient = df['orientation'].values,
        exper = df['Experiment'].values
    )

update()

# set up callbacks
widget_orient.on_change('active', lambda attr, old, new: update())

# layout
controls = [widget_orient]
sizing_mode = 'fixed' # 'scale_width' also looks nice
inputs = widgetbox(*controls, sizing_mode=sizing_mode)
l = layout([
        [inputs, p]
        ], sizing_mode=sizing_mode)

curdoc().add_root(l)
curdoc().title = "Arrhenius Diagram"