__author__ = 'Elco Luijendijk'

"""
Make a figure modeled temperatures and apatite (U-Th)/He ages.
Input consists of a python pickle (.pck) file generated by Beo

"""

###############
# load modules
###############

import matplotlib
matplotlib.use('Agg')

import string
import os
import pickle
import pdb
import itertools
import argparse

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.mlab
from matplotlib import ticker

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

print '-' * 50

parser = argparse.ArgumentParser(description='make figures of model temperatures and AHe ages for model runs Beo')

parser.add_argument(dest='output_files', metavar='output files', default=None, nargs='?',
                    help='one or more Beo output files (.pck)')

parser.add_argument('-t', dest='timeslices', default=None, help='one or more timeslices (years) to display in '
                                                                'the figure',
                    nargs='+', type=float)

parser.add_argument('-x', dest='xbounds', default=None, help='min. and max. x bound of figure as two numbers',
                    nargs='+', type=float)

parser.add_argument('-y', dest='ybounds', default=None, help='min. and max. y bound of figure as two numbers',
                    nargs='+', type=float)

parser.add_argument('-c', dest='combine_figs', help='combine modeled AHe ages or borehole temperatures in one figure',
                    action="store_true")

parser.add_argument('-m', dest='show_mesh', help='show mesh',
                    action="store_true")

parser.add_argument('-v', dest='do_not_show_vapour', help='do not show water vapour',
                    action="store_true")

parser.print_help()

print '\nnote, more options for figure layout are available in the figure_params.py file in ' \
      'the model_parameters directory'

print '-' * 50

args = parser.parse_args()

degree_symbol = unichr(176)
day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6


#
#xlim = [1500, 3500]
#ylim = [-2000, 0]

# read model output files
#model_output_folder = '/home/elco/model_files/hydrotherm_escript/'

result_dir = 'model_output'

if args.output_files is not None:

    files = args.output_files

else:
    # list files in model output directory
    files = os.listdir(result_dir)

    files = [os.path.join(result_dir, f) for f in files if f[-4:] == '.pck']

    files.sort(key=os.path.getmtime)

    files = files[::-1]

    print 'output files, from newest to oldest:'
    for i, fn in enumerate(files):
        print i, fn

    print 'enter a number to select a file, two numbers separated by - for a series of figures, ' \
          'or press enter to make a figure of all files'

    a = raw_input()

    if '-' in a:
        b = a.split('-')
        files = [files[i] for i in range(int(b[0]), int(b[1]) + 1)]
    elif ',' in a:
        b = a.split(',')
        files = [files[int(i)] for i in b]
    elif a != '':
        files = [files[int(a)]]

if args.do_not_show_vapour is True:
    fp.show_vapour = False

if args.show_mesh is True:
    fp.show_mesh = True

for fn in files:

    #fn = 'model_output/T_field_duration_500.pck'
    print 'reading model output datafile %s' % fn
    fin = open(fn, 'r')
    output_data = pickle.load(fin)
    fin.close()

    go = True
    try:
        if len(output_data) == 25:
            print 'reading output data, previous version of code, ' \
                  'without corrected AHe ages'
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

            Ahe_ages_all_corr = None
            AHe_ages_surface_corr = None
            AHe_ages_samples_surface_corr = None
        else:
            print 'reading output data, new version of code with corrected ' \
                  'AHe data'
            [runtimes_all, runtimes, xyz_array,
             surface_levels, x_loc_fault, z_loc_fault,
             T_init_array, T_array, boiling_temp_array,
             xyz_array_exc, exceed_boiling_temp_array, xyz_element_array,
             qh_array, qv_array, fault_fluxes,
             durations, xzs, Tzs,
             x_surface, T_surface,
             Ahe_ages_all, Ahe_ages_all_corr, xs_Ahe_all, Ahe_depths,
             AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
             AHe_ages_samples_surface, AHe_ages_samples_surface_corr, AHe_data_file,
             borehole_xlocs, borehole_zlocs, borehole_depths,
             borehole_temp_measured, borehole_temps_modeled] \
                = output_data

#4            [runtimes, runtimes[output_steps], xyz_array, surface_levels,
#1             T_init_array,
#2             T_array[output_steps], boiling_temp_array[output_steps],
#2             xyz_array_exc, exceed_boiling_temp_array[output_steps],
#1             xyz_element_array,
#2             qh_array[output_steps], qv_array[output_steps],
#2             fault_fluxes, durations,
#4             xzs, Tzs_cropped, x_surface, T_surface,
#4             AHe_ages_cropped, AHe_ages_corr_cropped, xs_Ahe_all, target_depths,
#3             AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
#2             AHe_ages_samples_surface, AHe_data_file]

        # 3  runtimes, runtimes[output_steps], xyz_array,
        # 1 surface_levels[output_steps],
        # 2 x_loc_fault, z_loc_fault,
        # 1 T_init_array,
        # 2 T_array, boiling_temp_array[output_steps],
        # 2 xyz_array_exc, exceed_boiling_temp_array[output_steps],
        # 1 xyz_element_array,
        # 2 qh_array[output_steps], qv_array[output_steps],
        # 2 fault_fluxes, durations,
        # 4 xzs, Tzs_cropped, x_surface, T_surface,
        # 4 AHe_ages_cropped, AHe_ages_corr_cropped, xs_Ahe_all, target_depths,
        # 3 AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
        # 3 AHe_ages_samples_surface, AHe_ages_samples_surface_corr, AHe_data_file,
        # 2 borehole_xlocs, borehole_zlocs,
        # 3 borehole_depths, borehole_temp_measured, borehole_temps_modeled]

    except ValueError:
        msg = 'error, could not read file %s' % fn
        msg += 'probably the versions of beo.py and make_figure.py do not match'
        raise ValueError(msg)
    # T_array, t_array, dx, dy, fault_mid, xi, yi, nt_heating,
    # subsurface_height, q_advective, duration_heating

    if args.timeslices is None:
        print 'saved timeslices for this model run:'
        for i, runtime in enumerate(runtimes):
            print '%i\t%0.2f yr' % (i, runtime / year)

        print '\nplease enter the timeslices you want to include. For several timesteps enter numbers divided by commas'

        key_inp = raw_input()
        if len(key_inp) == 1:
            fp.timeslices = [int(a)]
        else:
            inp_list = key_inp.split(',')
            fp.timeslices = [int(ts) for ts in inp_list]

    else:
        #times_raw = [float(ti) for ti in args.timeslices]

        timeslices = []
        rta = np.array(runtimes) / year
        for ti in args.timeslices:
            ind = np.argmin(np.abs(rta - ti))
            timeslices.append(ind)

        fp.timeslices = timeslices

    print 'timeslices to show: ', fp.timeslices

    xmin, xmax = xyz_array[:, 0].min(), xyz_array[:, 0].max()
    ymin, ymax = xyz_array[:, 1].min(), xyz_array[:, 1].max()

    print 'model dimensions:'
    print 'x coordinates: %0.1f to %0.1f ' % (xmin, xmax)
    print 'y coordinates: %0.1f to %0.1f ' % (ymin, ymax)

    if args.xbounds is None:
        print 'select figure bounds in x direction. Enter two numbers divided by a comma'
        print 'press enter to show full model domain'

        key_inp = raw_input()

        if key_inp == '':
            fp.xlim = None
        else:
            fp.xlim = np.array(key_inp.split(',')).astype(float)
    else:
        fp.xlim = args.xbounds

    if args.ybounds is None:
        print 'select figure bounds in y direction, press enter to show full model domain'
        key_inp = raw_input()
        if key_inp == '':
            fp.ylim = None
        else:
            fp.ylim = np.array(key_inp.split(',')).astype(float)
    else:
        fp.ylim = args.ybounds

    if borehole_xlocs is not None:
        print 'add panel with borehole temperatures (y/n) ?'
        if 'y' in raw_input():
            fp.add_temperature_panel = True
        else:
            fp.add_temperature_panel = False

    #
    show_surface_only = True
    if len(Ahe_ages_all) > 1:
        print 'found multiple depth levels with modeled AHe ages'
        print 'show only the surface level (y/n)?'

        if 'y' in raw_input():
            show_surface_only = True
        else:
            show_surface_only = False

    if go is True:

        # check if timeslices are really present
        nt = len(T_array)
        if nt == 1:
            fp.timeslices = [0]

        if type(fp.timeslices) is int:
            nts = fp.timeslices
            fp.timeslices = np.linspace(0, nt-1, nts).astype(int)

        fnew = []
        for f in fp.timeslices:
            if f >= nt:
                f = nt - 1
            fnew.append(f)
            fp.timeslices = fnew

        print 'making a figure of model run %s' % fn
        print 'at timeslices ', fp.timeslices
        print 'total number of saved timeslices = %i' % len(T_array)

        print 'x and y bounds of figure:'
        print 'x = ', fp.xlim
        print 'y = ', fp.ylim

        print 'using resolution for interpolation temperature field of %0.1f x %0.1f m' % (fp.dx, fp.dy)
        print 'this can be adjusted by modifying dx and dy in the figure_params.py file'

        nt = len(T_array)
        vmin = T_array.min()
        vmax = T_array.max()

        # temperature field contours and values:
        #fp.dx = 10.0
        #fp.dy = 10.0

        Tas = [T_array[ti] for ti in fp.timeslices]

        #pl.locator_params(nbins=3)
        #pl.locator_params(axis='y', nbins=3)
        #pl.locator_params(axis='x', nbins=3)

        #fig, panels = pl.subplots(1, 3, figsize=(8, 6), sharey=True)
        fig = pl.figure(figsize=(fp.xsize, fp.ysize))

        legs = []
        labels =[]
        #fig.subplots_adjust(wspace=0.05, hspace=0.05)

        import matplotlib.gridspec as gridspec
        nrows = 4
        ncols = len(Tas)
        height_ratios = fp.height_ratios
        width_ratios = [2] * ncols

        if fp.add_temperature_panel is True:
            ncols += 1
            width_ratios.append(2 * fp.relative_size_temp_panel)

        gs = gridspec.GridSpec(nrows, ncols,
                               height_ratios=height_ratios,
                               width_ratios=width_ratios)
        #ax = fig.add_subplot(gs[0, 0])

        if fp.add_temperature_panel is True:
            panels = [fig.add_subplot(gs[1, i]) for i in range(ncols - 1)]
            tpanels = [fig.add_subplot(gs[0, i]) for i in range(ncols - 1)]
            temp_panel = fig.add_subplot(gs[1, -1])
        else:
            panels = [fig.add_subplot(gs[1, i]) for i in range(ncols)]
            tpanels = [fig.add_subplot(gs[0, i]) for i in range(ncols)]

        gs.update(wspace=0.05, hspace=0.2)

        if Ahe_ages_all is not None:
            rpanels = [tp.twinx() for tp in tpanels]

        cpanel = fig.add_subplot(gs[-1, 0:fp.ncols_colorbar])
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

        print 'interpolating temperatures to regular grid'
        for p, Ta in zip(panels, Tas):
            xg, yg, zg = interpolate_data(xyz_array, Ta, fp.dx, fp.dy)
            leg_cn = p.contourf(xg, yg, zg, cnts, cmap=fp.cmap)
            #p.scatter(xyz_array[:, 0], xyz_array[:, 1], s=0.1, color='gray')
            #c=Ta, **kwargs)

        if fp.show_vapour is True and exceed_boiling_temp_array is not None:
            print 'showing location of water vapour'
            for p, bt in zip(panels, exceed_boiling_temp_array[fp.timeslices]):
                #ind = bt >= Ta
                #nind = ind == False
                ind = bt > 0.0
                nind = ind == False
                leg_vp = p.scatter(xyz_array_exc[:, 0][ind],
                                   xyz_array_exc[:, 1][ind],
                                   s=1, color='gray', alpha=0.3, zorder=201)
                #leg_vp = panels[-1].scatter(xyz_array[:, 0][nind],
                # xyz_array[:, 1][nind], s=1, color='gray')

            legs.append(leg_vp)
            labels.append('water vapour')

        if fp.show_mesh is True:
            for p in panels:
                p.scatter(xyz_array[:, 0], xyz_array[:, 1], s=0.25, color='black')

        for p, qhi, qvi in zip(panels, qh_array[fp.timeslices], qv_array[fp.timeslices]):
            print 'adding arrows'
            #xq, yq, qhg = interpolate_data(xyz_element_array[::fp.skip_arrows],
            #                               qhi[::fp.skip_arrows] * year,
            #                               fp.dx, fp.dy)
            #xq, yq, qvg = interpolate_data(xyz_element_array[::fp.skip_arrows],
            #                               qvi[::fp.skip_arrows] * year,
            #                               fp.dx, fp.dy)

            # set arrow scale
            #va = (qhg ** 2 + qvg ** 2) ** 0.5
            va = ((qhi * year)**2 + (qvi * year)**2) ** 0.5
            scale = np.abs(va).max() * fp.scale_multiplier
            print 'quiver scale = %0.2e' % scale

            # select only nodes with non-zero q
            ind = va != 0.0
            # and sort with increasing x coord:
            #a = np.argsort(xg[ind])

            #leg_q = p.quiver(xg[ind][a], yg[ind][a],
            #                 qhg[ind][a], qvg[ind][a],
            #                 angles='xy', scale=scale, headwidth=5, pivot='tip',
            #                 alpha=fp.arrow_transparency)
            leg_q = p.quiver(xyz_element_array[:, 0][ind][::fp.skip_arrows],
                             xyz_element_array[:, 1][ind][::fp.skip_arrows],
                             qhi[ind][::fp.skip_arrows] * year,
                             qvi[ind][::fp.skip_arrows] * year,
                             angles='xy', scale=scale, headwidth=5, pivot='tip',
                             alpha=fp.arrow_transparency)

            #if p == panels[-1]:
            #    p.quiverkey(leg_q, 0.85, 0.1, 'flow direction')

        if fp.add_temperature_panel is True:
            # show borehole locations
            for p, timeslice in zip(panels, fp.timeslices):

                for borehole_xloc, borehole_zloc, borehole_depth in \
                        zip(borehole_xlocs, borehole_zlocs, borehole_depths):

                    leg_bh, = p.plot((borehole_xloc[timeslice], borehole_xloc[timeslice]),
                                     (borehole_zloc[timeslice],
                                      borehole_zloc[timeslice] - borehole_depth.max()),
                                     color='gray', lw=1.5)

            legs.append(leg_bh)
            labels.append('borehole location')

        for p, timeslice in zip(panels, fp.timeslices):
            leg_surf = p.axhline(y=surface_levels[timeslice], color='black',
                                 lw=0.5)

        legs.append(leg_surf)
        labels.append('land surface')

        lss = fp.linestyles
        #colors = ['darkblue', 'green']
        # surface temperatures
        n_depths = len(Tzs)

        #for i in range(n_depths):
        for j, tp, timeslice in zip(itertools.count(), tpanels, fp.timeslices):
            leg_st, = tp.plot(x_surface[timeslice],
                              T_surface[timeslice],
                              color=fp.colors[j], ls=lss[j])

        legs.append(leg_st)
        labels.append('modelled temperature')

        if Ahe_ages_all is not None:
            #Ahe_ages_all, xs_Ahe_all = AHe_data
            #for i in range(n_depths):

            if args.combine_figs is True:
                n_repeats = len(files)
            else:
                n_repeats = 1

            AHe_label = 'modelled AHe ages'

            for n_repeat in range(n_repeats):

                # read additional data
                if args.combine_figs is True and n_repeat > 0:

                    fn = files[n_repeat]

                    print 'reading model output datafile %s for additional AHe model results' % fn
                    fin = open(fn, 'r')
                    output_data_add = pickle.load(fin)
                    fin.close()

                    [Ahe_ages_all, Ahe_ages_all_corr, xs_Ahe_all, Ahe_depths,
                     AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
                     AHe_ages_samples_surface, AHe_ages_samples_surface_corr, AHe_data_file] = output_data_add[20:30]
                    #AHe_xcoords_surface_all, AHe_ages_surface_corr_all, AHe_xcoords_surface_all, AHe_ages_surface_all

                    runtimes_new = output_data_add[1]

                    print 'new timeslices to show: '
                    print runtimes_new[fp.timeslices] / year

                    if len(Ahe_ages_all) > 1:
                        print 'found multiple depth levels with modeled AHe ages'
                        print 'show only the surface level (y/n)?'

                        if 'y' in raw_input():
                            show_surface_only = True
                        else:
                            show_surface_only = False

                if args.combine_figs is True:
                    print 'enter label for modeled AHe ages:'
                    AHe_label = raw_input()

                if show_surface_only is True:

                    for p, rp, timeslice in zip(panels, rpanels, fp.timeslices):

                        if (fp.show_corrected_AHe_ages is True
                            and AHe_ages_surface_corr is not None):

                            print 'plotting corrected AHe ages'
                            leg_ahe, = rp.plot(AHe_xcoords_surface[timeslice],
                                               AHe_ages_surface_corr[timeslice] / My,
                                               color=fp.AHe_color[n_repeat], ls=fp.AHe_linestyle)
                        else:
                            print 'plotting uncorrected AHe ages'
                            leg_ahe, = rp.plot(AHe_xcoords_surface[timeslice],
                                               AHe_ages_surface[timeslice] / My,
                                               color=fp.AHe_color[n_repeat], ls=fp.AHe_linestyle)

                        if timeslice == fp.timeslices[0]:
                            legs.append(leg_ahe)
                            labels.append(AHe_label)
                else:

                    for p, rp, timeslice in zip(panels, rpanels, fp.timeslices):

                        for AHe_age_i, AHe_age_corr_i, ahe_x, ahe_depth in zip(Ahe_ages_all,  Ahe_ages_all_corr,
                                                                               xs_Ahe_all, Ahe_depths):

                            print 'show AHe age for depth slice %0.1f m? (y/n)' % ahe_depth

                            if 'y' in raw_input():
                                if (fp.show_corrected_AHe_ages is True
                                        and AHe_ages_surface_corr is not None):

                                    print 'plotting corrected AHe ages'
                                    leg_ahe, = rp.plot(ahe_x,
                                                       AHe_age_corr_i[timeslice] / My,
                                                       ls=fp.AHe_linestyle)
                                else:
                                    print 'plotting uncorrected AHe ages'
                                    leg_ahe, = rp.plot(ahe_x,
                                                       AHe_age_i[timeslice] / My,
                                                       ls=fp.AHe_linestyle)

                                p.axhline(ahe_depth, ls=fp.AHe_linestyle, color=leg_ahe.get_color())

                                if timeslice == fp.timeslices[0]:
                                    legs.append(leg_ahe)
                                    labels.append('modelled AHe ages, %0.0f m' % ahe_depth)

        #AHe_ages_samples_surface, AHe_data_file
        if Ahe_ages_all is not None and AHe_data_file is not None:

            if len(np.unique(AHe_data_file['profile'].values)) > 1:
                print 'found multiple profiles in AHe data file'
                print 'select profile number to show in figure:'
                profile_no = int(raw_input())

                AHe_data_file = AHe_data_file[AHe_data_file['profile'] == profile_no]

            if fp.show_AHe_data is True:
                xr = AHe_data_file['distance_to_fault'].values
                y = AHe_data_file['AHe_age_corr'].values
                yerr = AHe_data_file['AHe_age_corr_2se'].values

                for rp, timeslice in zip(rpanels, fp.timeslices):
                    x = x_loc_fault[timeslice] + xr

                    if fp.show_average_AHe_ages is False:
                        leg_ahe_samples = rp.errorbar(x, y, yerr=yerr,
                                                      marker='o',
                                                      ms=fp.marker_size,
                                                      markeredgecolor='black',
                                                      color=fp.AHe_color[0],
                                                      linestyle='None')
                        #legs.append(leg_ahe_samples)
                        #labels.append('measured AHe ages')

                    else:
                        leg_ahe_samples_single = rp.scatter(x, y,
                                                     marker='o',
                                                     s=4,
                                                     color='gray')

                        sample_names_adj = [s.split(' ')[0] for s in AHe_data_file['sample'].values]
                        AHe_data_file['sample_names_adj'] = sample_names_adj
                        sample_ids = np.unique(AHe_data_file['sample_names_adj'])


                        avg_ages = []
                        avg_age_ses = []
                        avg_dists = []
                        for sample_id in sample_ids:
                            ind = AHe_data_file['sample_names_adj'] == sample_id
                            avg_age_i = np.mean(AHe_data_file.loc[ind, 'AHe_age_corr'])
                            avg_age_se_i = np.mean(AHe_data_file.loc[ind, 'AHe_age_corr_2se'])
                            avg_dist_i = np.mean(AHe_data_file.loc[ind, 'distance_to_fault'])

                            avg_ages.append(avg_age_i)
                            avg_age_ses.append(avg_age_se_i)
                            avg_dists.append(avg_dist_i)

                        x = x_loc_fault[timeslice] + np.array(avg_dists)
                        leg_ahe_samples = rp.errorbar(x, avg_ages, yerr=avg_age_ses,
                                                      marker='o',
                                                      ms=fp.marker_size,
                                                      markerfacecolor='None',
                                                      markeredgecolor=fp.AHe_color[0],
                                                      color=fp.AHe_color[0],
                                                      linestyle='None',
                                                      lw=1.0)

                if fp.show_average_AHe_ages is False:
                    legs.append(leg_ahe_samples)
                    labels.append('measured AHe ages')

                else:

                    legs.append(leg_ahe_samples_single)
                    labels.append('measured single grain AHe ages')
                    legs.append(leg_ahe_samples)
                    labels.append('average AHe ages samples')

        if fp.add_temperature_panel is True:

            for borehole_depth, borehole_temp_modeled, borehole_temp_measured_i in \
                    zip(borehole_depths, borehole_temps_modeled, borehole_temp_measured):

                leg_bh_temp_meas, = temp_panel.plot(borehole_temp_measured_i, -borehole_depth,
                                color='gray', lw=1.5, ls='--')

                for j, timeslice in enumerate(fp.timeslices):

                    leg_bh_temp_mod, = temp_panel.plot(borehole_temp_modeled[timeslice],
                                                       -borehole_depth,
                                                       lw=1.0,
                                                       ls=lss[j],
                                                       color=fp.colors[j])

                legs += [leg_bh_temp_meas]
                labels += ['measured temperature']

        panels[0].set_ylabel('Elevation (m)')

        for p in panels:

            p.yaxis.grid(False)

            if fp.xlim is not None:
                if fp.xlim[0] > xmin:
                    xmin = fp.xlim[0]
                if fp.xlim[1] < xmax:
                    xmax = fp.xlim[1]

            if fp.ylim is not None:
                if fp.ylim[0] > ymin:
                    ymin = fp.ylim[0]
                if fp.ylim[1] < ymax:
                    ymax = fp.ylim[1]

            p.set_xlim(xmin, xmax)
            p.set_ylim(ymin, ymax)

            p.set_xlabel('Distance (m)')

            #p.set_xticks(p.get_xticks()[:-1])

        for p, tp in zip(panels[1:], tpanels[1:]):
            p.set_yticklabels([])
            tp.set_yticklabels([])

        for tp in tpanels[:]:
            tp.set_xticklabels([])
            maxT = np.max(np.hstack(T_surface)) * 1.1
            tp.set_ylim(0, maxT)
            tp.yaxis.grid(False)
            tp.set_xlim(xmin, xmax)
            #tp.set_xticks(tp.get_xticks()[:-1])

        for tp in tpanels[:]:
            tp.spines['top'].set_visible(False)
            tp.get_xaxis().tick_bottom()

            if Ahe_ages_all is None:
                tp.spines['right'].set_visible(False)
                tp.get_yaxis().tick_left()

        if Ahe_ages_all is not None:

            for rp in rpanels[:-1]:
                rp.set_yticklabels([])

            for rp in rpanels:
                rp.set_xlim(xmin, xmax)
                if Ahe_ages_all is not None:
                    rp.set_ylim(0, Ahe_ages_all[0].max() / My * 1.5)
                rp.spines['top'].set_visible(False)
                rp.get_xaxis().tick_bottom()

        titles = []

        for dti in runtimes[fp.timeslices]:
            titles += ['%0.0f years' % (dti / year)]

        for tp, title in zip(tpanels, titles):
            tp.set_title(title)

        if fp.add_temperature_panel is True:
            temp_panel.set_yticklabels([])
            temp_panel.set_ylim(ymin, ymax)
            temp_panel.spines['top'].set_visible(False)
            temp_panel.get_xaxis().tick_bottom()
            temp_panel.spines['right'].set_visible(False)
            temp_panel.get_yaxis().tick_left()
            temp_panel.set_xlabel('Borehole\ntemperature (%sC)' % degree_symbol)

        #for panel in panels:
        #    panel.set_yticks(panel.get_yticks()[::2])

        #for tpanel in tpanels:
        #    tpanel.set_yticks(tpanel.get_yticks()[::2])

        #if Ahe_ages_all is not None:
        #    for rpanel in rpanels:
        #        rpanel.set_yticks(rpanel.get_yticks()[::2])

        tpanels[0].set_ylabel('Surface\ntemperature (%sC)' % degree_symbol)

        if Ahe_ages_all is not None:
            rpanels[-1].set_ylabel('AHe age (My)')
            rpanels[-1].yaxis.label.set_color(fp.AHe_color[0])
            rpanels[-1].tick_params(axis='y', colors=fp.AHe_color[0])

        for p, tp in zip(panels, tpanels):
            p.set_xlim(xmin, xmax)
            tp.set_xlim(xmin, xmax)
            p.set_ylim(ymin, ymax)

        if fp.add_temperature_panel is True:
            temp_panel.set_ylim(ymin, ymax)

        cb = fig.colorbar(leg_cn, cax=cpanel, orientation='horizontal')
        tick_locator = ticker.MaxNLocator(nbins=fp.bins_colorbar)
        cb.locator = tick_locator
        cb.update_ticks()
        cb.set_label('Temperature (%s C)' % degree_symbol)

        if fp.add_legend is True:
            leg = leg_panel.legend(legs, labels, frameon=False,
                                   fontsize=fp.legend_font_size,
                                   loc='upper right',
                                   ncol=fp.n_columns_legend)
            #leg = fig.legend(legs, labels, frameon=False,
            #                 loc='lower right', fontsize=fp.legend_font_size)

        fn_fig = fn[:-4] + '_T_field.%s' % fp.fig_format

        fn_local1 = os.path.split(fn_fig)[-1]
        fn_local2 = os.path.join('model_output', fn_local1)

        #gs.tight_layout(fig)

        all_panels = tpanels + panels
        if fp.add_temperature_panel is True:
            all_panels.append(temp_panel)

        for i, p in enumerate(all_panels):
            p.text(0.01, 1.01, string.ascii_lowercase[i],
                   weight='bold', transform=p.transAxes, ha='left', va='bottom', fontsize='medium')

        for p in all_panels:

            locy = ticker.MaxNLocator(nbins=fp.bins_yaxis)  # this locator puts ticks at regular intervals
            locx = ticker.MaxNLocator(nbins=fp.bins_xaxis)  # this locator puts ticks at regular intervals

            p.xaxis.set_major_locator(locx)
            p.yaxis.set_major_locator(locy)

            #temp_panel.set_yticks(temp_panel.get_yticks()[::2])

        print 'saving %s' % fn_fig
        fig.savefig(fn_local2, dpi=fp.figure_resolution)

        pl.clf()

        print 'done with model scenario'

print 'done'
