"""
make a figure of spring T vs fault depth and background geothermal gradient
"""

import numpy as np

import pandas as pd

import matplotlib.pyplot as pl
import matplotlib.cm


def simpleaxis(ax, removeh=True):
    """
    remove top and right axis from panel
    :param ax:
    :return:
    """

    ax.spines['top'].set_visible(False)
    if removeh == True:
        ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    return


__author__ = 'elco'

degree_symbol = unichr(176)

target_temperature = 50.0

cmap = matplotlib.cm.seismic

fn = '/home/elco/python_scripts/hydrotherm_escript/model_output/' \
     'model_params_and_results_48_runs.csv'

df = pd.read_csv(fn)

##################################################
fig, panel = pl.subplots(1, 1)

#x = df['bottom_temperature'] / df['total_depth']
x = df['thermal_gradient']
a = df['fault_bottoms']
y = np.array([-float(i.replace('[', '').replace(']', '')) for i in a])
z = df['max_surface_temperature'] - target_temperature

leg_sc = panel.scatter(x, y, c=z, cmap=cmap, s=40, vmin=-25.0, vmax=25.0)

cb = pl.colorbar(leg_sc, shrink=0.5)
cb.set_label('Model error\nspring temperature (%sC)' % degree_symbol)

panel.set_xlabel('Geothermal gradient (%sC/m)' % degree_symbol)
panel.set_ylabel('Depth of flow path (m)')

panel.yaxis.grid(True)
panel.xaxis.grid(False)
panel.set_xlim(x.min() * 0.9, x.max() * 1.1)

simpleaxis(panel)

fnf = fn[:-4] + '_fig.png'

fig.savefig(fnf, dpi=200)

print 'done'