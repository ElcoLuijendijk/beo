__author__ = 'elco'

"""
simple 2D model of advective heat flow

Elco Luijendijk, Goettingen University, 2015

"""

###############
# load modules
###############

import os
import pickle
import itertools

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.mlab

import useful_functions


def interpolate_data(xyz_array, Ti, dx, dy):

    xi = np.arange(xyz_array[:, 0].min(), xyz_array[:, 0].max() + dx, dx)
    yi = np.arange(xyz_array[:, 1].min(), xyz_array[:, 1].max() + dy, dy)
    xg, yg = np.meshgrid(xi, yi)
    xgf, ygf = xg.flatten(), yg.flatten()
    #zgf = scipy.interpolate.griddata(xyz_array, Ti, np.vstack((xgf, ygf)).T,
    #                                 method='linear')
    zg = matplotlib.mlab.griddata(xyz_array[:, 0], xyz_array[:, 1], Ti, xi, yi, interp='linear')

    return xg, yg, zg


degree_symbol = unichr(176)
day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6

cnt_int = 10.0

dx = 10.0
dy = 10.0

# add extra row with zoom-in of model results
add_zoom_panel = False

# extent for zoom panel
xlim = [1500, 3500]
ylim = [-2000, 100]

# timesteps to select for output
#timeslices = [2, 20, 100]
timeslices = [0]

# input files
files = ['model_output/T_field_model_run_14_(-3000.0, 0.03).pck',
         'model_output/T_field_model_run_32_(-6000.0, 0.03).pck']

labels = ['flow path depth = 3000 m',
          'flow path depth = 6000 m',
          ]

xsize = 8
golden_ratio = (1.0 + np.sqrt(5))/2.0
ysize = xsize / golden_ratio
if add_zoom_panel is True:
    ysize = xsize

fig = pl.figure(figsize=(xsize, ysize))

#fig.subplots_adjust(wspace=0.05, hspace=0.05)

import matplotlib.gridspec as gridspec
nrows = 4
hr = [25, 75, 12, 3]
hspace = 0.1

if add_zoom_panel is True:
    nrows += 1
    hr = [25, 75, 75, 9, 3]
    hspace = 0.3

ncols = len(files)
gs = gridspec.GridSpec(nrows, ncols,
                       height_ratios=hr)

panels = [fig.add_subplot(gs[1, i]) for i in range(ncols)]
tpanels = [fig.add_subplot(gs[0, i]) for i in range(ncols)]

if add_zoom_panel is True:
    zpanels = [fig.add_subplot(gs[2, i]) for i in range(ncols)]
else:
    zpanels = [None] * len(panels)

cpanel = fig.add_subplot(gs[-1, :2])
cpanel.set_xticks([])
cpanel.set_yticks([])

#rpanels = [tp.twinx() for tp in tpanels]

gs.update(wspace=0.2, hspace=hspace)

for i, panel, tpanel, zpanel, fn in zip(itertools.count(), panels, tpanels,
                                        zpanels, files):

    #fn = 'model_output/T_field_duration_500.pck'
    print 'loading data from %s' % fn
    fin = open(fn, 'r')
    output_data = pickle.load(fin)
    fin.close()

    [runtimes, xyz_array, T_steady, T_array, xyz_element_array,
     qh_array, qv_array,
     fault_fluxes, durations, xzs, Tzs, AHe_data] = output_data

    timeslice = -1
    print 'making a figure of model run %s' % fn
    print 'at timeslices ', timeslices
    print 'total number of saved timeslices = %i' % len(T_array)

    xmin, xmax = xyz_array[:, 0].min(), xyz_array[:, 0].max()
    ymin, ymax = xyz_array[:, 1].min(), xyz_array[:, 1].max()

    nt = len(T_array)

    # temperature field contours and values:
    T_slice = T_array[timeslice]
    Ta = T_slice - T_steady

    vmin = Ta.min()
    vmax = Ta.max()

    xg, yg, zg = interpolate_data(xyz_array, Ta, dx, dy)
    cnts = np.arange(vmin, vmax+cnt_int, cnt_int)
    leg_cn = panel.contourf(xg, yg, zg, cnts)

    if add_zoom_panel is True:
        leg_cn = zpanel.contourf(xg, yg, zg, cnts)


    print 'adding arrows'
    qhi = qh_array[timeslice]
    qvi = qv_array[timeslice]
    xq, yq, qhg = interpolate_data(xyz_element_array, qhi * year, dx, dy)
    xq, yq, qvg = interpolate_data(xyz_element_array, qvi * year, dx, dy)

    # set arrow scale
    va = (qhg ** 2 + qvg ** 2) ** 0.5
    scale = np.abs(va).max() * 20.0
    print 'quiver scale = %0.2e' % scale
    # select only nodes with non-zero q
    ind = va != 0.0
    # and sort with increasing x coord:
    a = np.argsort(xg[ind])

    thin = 40
    leg_q = panel.quiver(xg[ind][a][::thin], yg[ind][a][::thin],
                         qhg[ind][a][::thin], qvg[ind][a][::thin],
                         angles='xy', scale=scale, headwidth=5)

    if add_zoom_panel is True:
        thin2 = 10
        leg_q = zpanel.quiver(xg[ind][a][::thin2], yg[ind][a][::thin2],
                              qhg[ind][a][::thin2], qvg[ind][a][::thin2],
                              angles='xy', scale=scale, headwidth=5)

    lss = ['-', ':']
    i = 0
    leg_st, = tpanel.plot(xzs[i], Tzs[i][timeslice],
                          color='black', ls=lss[i])

    #
    #if AHe_data is not None:
    #    Ahe_ages_all, xs_Ahe_all = AHe_data
    #    for i in range(n_depths):
    #        for rp, timeslice in zip(rpanels, timeslices):
    #            leg_ahe, = rp.plot(xs_Ahe_all[i],
    #                               Ahe_ages_all[i][timeslice] / My,
    #                               color='blue', ls=lss[i])

panels[0].set_ylabel('Elevation (m)')
if add_zoom_panel is True:
    zpanels[0].set_ylabel('Elevation (m)')

for p in panels:

    # horizontal line at surface
    p.axhline(y=0, color='black')
    p.yaxis.grid(True)

    p.set_xlim(0, xyz_array[:, 0].max())
    p.set_ylim(xyz_array[:, 1].min(), xyz_array[:, 1].max())
    p.set_xticks(p.get_xticks()[:-1])
    if add_zoom_panel is False:
        p.set_xlabel('Distance (m)')
    useful_functions.simpleaxis(p)

if add_zoom_panel is True:
    for zp in zpanels:

        # horizontal line at surface
        zp.axhline(y=0, color='black')
        zp.yaxis.grid(True)
        zp.set_xlim(xlim[0], xlim[1])
        zp.set_ylim(ylim[0], ylim[1])
        zp.set_xlabel('Distance (m)')
        useful_functions.simpleaxis(p)

for p, tp, zp in zip(panels[1:], tpanels[1:], zpanels[1:]):
    p.set_yticklabels([])
    tp.set_yticklabels([])
    if add_zoom_panel is True:
        zp.set_yticklabels([])

for tp in tpanels[:]:
    tp.set_xticklabels([])
    tp.set_ylim(0, Tzs[0].max() * 1.1)
    tp.yaxis.grid(True)

for tp in tpanels[:]:
    tp.spines['top'].set_visible(False)
    tp.get_xaxis().tick_bottom()

#for rp in rpanels[:-1]:
#    rp.set_yticklabels([])

#for rp in rpanels:
#    rp.set_xlim(xlim)
#    if AHe_data is not None:
#        rp.set_ylim(0, Ahe_ages_all[-1].max() / My * 1.1)
#    rp.spines['top'].set_visible(False)
#    rp.get_xaxis().tick_bottom()

    #useful_functions.simpleaxis(rp, removeh=False)

#for tp, dt in zip(tpanels, runtimes[timeslices]):
#    dti = dt
    #dti = dt - runtimes[10]
    #title = '%0.0f years' % (dti / year)
    #tp.set_title(title)

for panel in panels:
    #panel.set_yticks(panel.get_yticks()[::2])
    pass
for tpanel in tpanels:
    tpanel.set_yticks(tpanel.get_yticks()[::2])
if add_zoom_panel is True:
    for zpanel in zpanels:
        zpanel.set_yticks(zpanel.get_yticks()[::2])
#for rpanel in rpanels:
#    rpanel.set_yticks(rpanel.get_yticks()[::2])

for panel, tpanel, zpanel in zip(panels, tpanels, zpanels):
    panel.set_xticks(panel.get_xticks()[::2])
    tpanel.set_xticks(tpanel.get_xticks()[::2])
    if add_zoom_panel is True:
        #zpanel.set_xticks(panel.get_xticks()[::2])
        pass

for tp, label in zip(tpanels, labels):
    tp.set_title(label)

for p in panels:
    p.set_ylim(xyz_array[:, 1].min(), xyz_array[:, 1].max())

tpanels[0].set_ylabel('Surface\ntemperature (%sC)' % degree_symbol)
#rpanels[-1].set_ylabel('U-Th/He age (My)')
#rpanels[-1].yaxis.label.set_color('blue')
#rpanels[-1].tick_params(axis='y', colors='blue')

cb = fig.colorbar(leg_cn, cax=cpanel, orientation='horizontal')
cb.set_label('Temperature change (%sC)' % degree_symbol)

fn_fig = fn[:-4] + '_T_field'
#if steady_state is True:
#    fn_fig += '_ss'
fn_fig += '.png'
print 'saving %s' % fn_fig
fig.savefig(fn_fig, dpi=150)

fn_local1 = os.path.split(fn_fig)[-1]
fn_local2 = os.path.join('model_output', fn_local1)

print 'saving %s' % fn_fig
fig.savefig(fn_local2, dpi=150)

pl.clf()

print 'done'
