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
import pdb

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.mlab

import model_parameters.figure_params as fp


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


def interpolate_data(xyz_array, Ti, dx, dy):

    xi = np.arange(xyz_array[:, 0].min(), xyz_array[:, 0].max() + dx, dx)
    yi = np.arange(xyz_array[:, 1].min(), xyz_array[:, 1].max() + dy, dy)
    xg, yg = np.meshgrid(xi, yi)
    xgf, ygf = xg.flatten(), yg.flatten()
    #zgf = scipy.interpolate.griddata(xyz_array, Ti, np.vstack((xgf, ygf)).T,
    #                                 method='linear')

    zg = matplotlib.mlab.griddata(xyz_array[:, 0], xyz_array[:, 1], Ti,
                                  xi, yi,
                                  interp='linear')

    return xg, yg, zg


degree_symbol = unichr(176)
day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6

# timesteps to select for output
#timeslices = [2, 20, 100]


#
#xlim = [1500, 3500]
#ylim = [-2000, 0]

# read model output files
#model_output_folder = '/home/elco/model_files/hydrotherm_escript/'
result_dir = 'model_output'
files = os.listdir(result_dir)

files = [os.path.join(result_dir, f) for f in files if f[-4:] == '.pck']

files.sort(key=os.path.getmtime)

files = files[::-1]

print 'output files, from newest to oldest:'
for i, file in enumerate(files):
    print i, file

print 'enter a number to select a file, ' \
      'or press enter to make a figure of all files'

a = raw_input()
if a != '':
    files = [files[int(a)]]


for fn in files:

    #fn = 'model_output/T_field_duration_500.pck'
    fin = open(fn, 'r')
    output_data = pickle.load(fin)
    fin.close()

    go = True
    try:
        #[runtimes, xyz_array, T_array, xyz_element_array, qh_array, qv_array,
        # fault_fluxes, durations, xzs, Tzs, AHe_data] = output_data
        [runtimes_all, runtimes, xyz_array, surface_levels,
         T_init_array, T_array, boiling_temp_array,
         xyz_array_exc, exceed_boiling_temp_array,
         xyz_element_array,
         qh_array, qv_array,
         fault_fluxes, durations,
         xzs, Tzs, x_surface, T_surface,
         Ahe_ages_all, xs_Ahe_all, Ahe_depths,
         AHe_ages_surface, AHe_xcoords_surface,
         AHe_ages_samples_surface, AHe_data_file] \
            = output_data

    except ValueError:
        msg = 'error, could not read file %s' % fn
        msg += 'probably the versions of beo.py and make_figure.py do not match'
        raise ValueError(msg)
    # T_array, t_array, dx, dy, fault_mid, xi, yi, nt_heating,
    # subsurface_height, q_advective, duration_heating

    if go is True:

        print 'making a figure of model run %s' % fn
        print 'at timeslices ', fp.timeslices
        print 'total number of saved timeslices = %i' % len(T_array)

        xmin, xmax = xyz_array[:, 0].min(), xyz_array[:, 0].max()
        ymin, ymax = xyz_array[:, 1].min(), xyz_array[:, 1].max()

        nt = len(T_array)
        vmin = T_array.min()
        vmax = T_array.max()

        # temperature field contours and values:
        #fp.dx = 10.0
        #fp.dy = 10.0

        Tas = [T_array[ti] for ti in fp.timeslices]

        #fig, panels = pl.subplots(1, 3, figsize=(8, 6), sharey=True)
        fig = pl.figure(figsize=(fp.xsize, fp.ysize))

        legs = []
        labels =[]
        #fig.subplots_adjust(wspace=0.05, hspace=0.05)

        import matplotlib.gridspec as gridspec
        nrows = 4
        ncols = len(fp.timeslices)
        gs = gridspec.GridSpec(nrows, ncols, height_ratios=[25, 75, 17, 3])
        #ax = fig.add_subplot(gs[0, 0])

        panels = [fig.add_subplot(gs[1, i]) for i in range(ncols)]
        tpanels = [fig.add_subplot(gs[0, i]) for i in range(ncols)]

        gs.update(wspace=0.05, hspace=0.2)

        rpanels = [tp.twinx() for tp in tpanels]

        cpanel = fig.add_subplot(gs[-1, :-1])
        cpanel.set_xticks([])
        cpanel.set_yticks([])

        leg_panel = fig.add_subplot(gs[-1:, -1], frameon=False)
        leg_panel.set_xticks([])
        leg_panel.set_yticks([])

        kwargs = dict(edgecolor='None', vmin=vmin, vmax=vmax)

        #cnt_int = 10.0
        if fp.cnt_range is None:
            cnts = np.arange(vmin, vmax+fp.cnt_int, fp.cnt_int)
        else:
            cnts = np.arange(fp.cnt_range[0], fp.cnt_range[-1] + fp.cnt_int, fp.cnt_int)

        for p, Ta in zip(panels, Tas):
            xg, yg, zg = interpolate_data(xyz_array, Ta, fp.dx, fp.dy)
            leg_cn = p.contourf(xg, yg, zg, cnts)
            #p.scatter(xyz_array[:, 0], xyz_array[:, 1], s=0.1, color='gray')
            #c=Ta, **kwargs)

        if fp.show_vapour is True:
            print 'showing location of water vapour'
            for p, bt in zip(panels, exceed_boiling_temp_array[fp.timeslices]):
                #ind = bt >= Ta
                #nind = ind == False
                ind = bt > 0.0
                nind = ind == False
                leg_vp = p.scatter(xyz_array_exc[:, 0][ind],
                                   xyz_array_exc[:, 1][ind],
                                   s=1, color='gray', alpha=0.3)
                #leg_vp = panels[-1].scatter(xyz_array[:, 0][nind],
                # xyz_array[:, 1][nind], s=1, color='gray')

            legs.append(leg_vp)
            labels.append('water vapour')

        for p, qhi, qvi in zip(panels, qh_array, qv_array):
            print 'adding arrows'
            xq, yq, qhg = interpolate_data(xyz_element_array, qhi * year,
                                           fp.dx, fp.dy)
            xq, yq, qvg = interpolate_data(xyz_element_array, qvi * year,
                                           fp.dx, fp.dy)

            # set arrow scale
            va = (qhg ** 2 + qvg ** 2) ** 0.5
            scale = np.abs(va).max() * fp.scale_multiplier
            print 'quiver scale = %0.2e' % scale

            # select only nodes with non-zero q
            ind = va != 0.0
            # and sort with increasing x coord:
            a = np.argsort(xg[ind])

            leg_q = p.quiver(xg[ind][a][::fp.skip_arrows], yg[ind][a][::fp.skip_arrows],
                             qhg[ind][a][::fp.skip_arrows], qvg[ind][a][::fp.skip_arrows],
                             angles='xy', scale=scale, headwidth=5, pivot='tip',
                             alpha=fp.arrow_transparency)

            #legs.append(leg_q)
            #labels.append('flow arrows')

        lss = ['-'] * 50
        #colors = ['darkblue', 'green']
        # surface temperatures
        n_depths = len(Tzs)

        #for i in range(n_depths):
        for tp, timeslice in zip(tpanels, fp.timeslices):
            leg_st, = tp.plot(x_surface[timeslice],
                              T_surface[timeslice],
                              color='black', ls=lss[i])

        legs.append(leg_st)
        labels.append('surface temperature')

        for p, timeslice in zip(panels, fp.timeslices):
            p.axhline(y=surface_levels[timeslice], color='black', lw=0.25)

        #
        if Ahe_ages_all is not None:
            #Ahe_ages_all, xs_Ahe_all = AHe_data
            #for i in range(n_depths):
            for rp, timeslice in zip(rpanels, fp.timeslices):

                leg_ahe, = rp.plot(AHe_xcoords_surface[timeslice],
                                   AHe_ages_surface[timeslice] / My,
                                   color=fp.AHe_color, ls=lss[i])

            legs.append(leg_ahe)
            labels.append('modeled AHe ages')

        #AHe_ages_samples_surface, AHe_data_file
        if AHe_data_file is not None:

            if len(np.unique(AHe_data_file['profile'].values)) > 1:
                print 'found multiple profiles in AHe data file'
                print 'select profile number to show in figure:'
                profile_no = int(raw_input())

                AHe_data_file = AHe_data_file[AHe_data_file['profile'] == profile_no]

            x = AHe_data_file['distance'].values
            y = AHe_data_file['AHe_age_uncorr'].values
            yerr = AHe_data_file['AHe_age_uncorr_2se'].values

            for rp, timeslice in zip(rpanels, fp.timeslices):
                leg_ahe_samples = rp.errorbar(x, y, yerr=yerr,
                                              marker='o',
                                              ms=fp.marker_size,
                                              markeredgecolor='black',
                                              color=fp.AHe_color,
                                              linestyle='None')

            legs.append(leg_ahe_samples)
            labels.append('measured AHe ages')

        #for p in panels[:]:
        panels[0].set_ylabel('Elevation (m)')
        #bpanels[0].set_ylabel('Elevation (m)')

        for p in panels:

            # horizontal line at surface
            #p.axhline(y=0, color='black')

            p.yaxis.grid(True)

            #p.set_xlim(0, xyz_array[:, 0].max())

            if fp.xlim[0] > xmin:
                xmin = fp.xlim[0]
            if fp.xlim[1] < xmax:
                xmax = fp.xlim[1]
            if fp.ylim[0] > ymin:
                ymin = fp.ylim[0]
            if fp.ylim[1] < ymax:
                ymax = fp.ylim[1]

            p.set_xlim(xmin, xmax)
            p.set_ylim(ymin, ymax)
            p.set_xticks(p.get_xticks()[:-1])
            p.set_xlabel('Distance (m)')

        for p, tp in zip(panels[1:], tpanels[1:]):
            p.set_yticklabels([])
            tp.set_yticklabels([])

        for tp in tpanels[:]:
            tp.set_xticklabels([])
            tp.set_ylim(0, T_surface[-1].max() * 1.1)
            tp.yaxis.grid(True)
            tp.set_xlim(xmin, xmax)

        for tp in tpanels[:]:
            tp.spines['top'].set_visible(False)
            tp.get_xaxis().tick_bottom()
        #useful_functions.simpleaxis(tpanels[-1], removeh=False)

        for rp in rpanels[:-1]:
            rp.set_yticklabels([])

        for rp in rpanels:
            rp.set_xlim(xmin, xmax)
            if Ahe_ages_all is not None:
                rp.set_ylim(0, Ahe_ages_all[0].max() / My * 1.25)
            rp.spines['top'].set_visible(False)
            rp.get_xaxis().tick_bottom()

            #useful_functions.simpleaxis(rp, removeh=False)

        for tp, dt in zip(tpanels, runtimes[fp.timeslices]):
            dti = dt
            #dti = dt - runtimes[10]
            title = '%0.0f years' % (dti / year)
            tp.set_title(title)

        for panel in panels:
            panel.set_yticks(panel.get_yticks()[::2])
        for tpanel in tpanels:
            tpanel.set_yticks(tpanel.get_yticks()[::2])
        for rpanel in rpanels:
            rpanel.set_yticks(rpanel.get_yticks()[::2])

        for panel, tpanel in zip(panels, tpanels):
            panel.set_xticks(panel.get_xticks()[::2])
            tpanel.set_xticks(tpanel.get_xticks()[::2])

        tpanels[0].set_ylabel('Temperature (%s C)' % degree_symbol)
        rpanels[-1].set_ylabel('AHe age (My)')
        rpanels[-1].yaxis.label.set_color(fp.AHe_color)
        rpanels[-1].tick_params(axis='y', colors=fp.AHe_color)

        for p in panels:
            p.set_xlim(xmin, xmax)
            p.set_ylim(ymin, ymax)

        cb = fig.colorbar(leg_cn, cax=cpanel, orientation='horizontal')
        cb.set_label('Temperature (%s C)' % degree_symbol)

        if fp.add_legend is True:
            leg = leg_panel.legend(legs, labels, frameon=False,
                                   fontsize=fp.legend_font_size,
                                   loc='upper center')
            #leg = fig.legend(legs, labels, frameon=False,
            #                 loc='lower right', fontsize=fp.legend_font_size)

        fn_fig = fn[:-4] + '_T_field.%s' % fp.fig_format

        fn_local1 = os.path.split(fn_fig)[-1]
        fn_local2 = os.path.join('model_output', fn_local1)

        print 'saving %s' % fn_fig
        fig.savefig(fn_local2, dpi=fp.figure_resolution)

        pl.clf()

        print 'done with model scenario'

print 'done'
