"""
run multiple models sequentially with different parameter values

"""

__author__ = 'Elco Luijendijk'

import sys
import os
import pickle
import itertools
import inspect
import pdb
import datetime
import scipy.interpolate

import numpy as np
import pandas as pd

from model_parameters.model_parameters import ModelParams
import model_parameters.parameter_ranges as pr
import beo


def interpolate_data(xyz, data, target_x, y_steps=300):

    """

    """

    xi = np.array([xyz[:, 0].min(), target_x, xyz[:, 0].max()])
    yi = np.linspace(xyz[:, 1].min(), xyz[:, 1].max(), y_steps)
    nx = len(xi)
    ny = len(yi)

    xg, yg = np.meshgrid(xi, yi)
    xgf, ygf = xg.flatten(), yg.flatten()

    # interpolate u to grid
    zgf = scipy.interpolate.griddata(xyz, data, np.vstack((xgf, ygf)).T,
                                     method='linear')

    # make a 2d grid again
    zg = np.resize(zgf, (ny, nx))

    return yg[:, 1], zg[:, 1]


def coefficient_of_determination(y, f):

    """
    calculate coefficient of determination

    Parameters
    ----------
    y : numpy array
        observed values
    f : numpy array
        calculated values

    Returns
    -------
    R2 : float
        coefficient of determination

    """

    ind = (np.isnan(y)==False) & (np.isnan(f)==False) & \
        (y != -np.inf) & (y != np.inf) &\
        (f != -np.inf) & (f != np.inf)

    yc = y[ind]
    fc = f[ind]

    # residual sum of squares
    SS_res = np.sum((yc - fc)**2)
    # total sum of squares
    SS_tot = np.sum((yc - np.mean(yc))**2)

    R2 = 1.0 - (SS_res / SS_tot)

    return R2


mp = ModelParams

day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6

scriptdir = os.path.realpath(sys.path[0])
output_folder = os.path.join(scriptdir, mp.output_folder)

# create list with param values for each model run
scenario_param_names_raw = dir(pr)
scenario_param_names = [m for m in scenario_param_names_raw
                        if '__' not in m and '_s' in m]

scenario_parameter_list = [getattr(pr, p)
                           for p in scenario_param_names]

# construct list with all parameter combinations
if pr.parameter_combinations is True:
    scenario_parameter_combinations = \
        list(itertools.product(*scenario_parameter_list))
else:
    nscens = np.sum(np.array([len(sp) for sp in scenario_parameter_list
                              if sp is not None]))
    nparams = len(scenario_parameter_list)
    scenario_parameter_combinations = []

    if pr.initial_base_run is True:
        scenario_parameter_combinations.append([None] * nparams)

    for j, sl in enumerate(scenario_parameter_list):
        if sl[0] is not None:
            sc = [None] * nparams
            for sli in sl:
                sci = list(sc)
                sci[j] = sli
                scenario_parameter_combinations.append(sci)

param_list = scenario_parameter_combinations

# read default model parameter file
Parameters = mp()


# get attributes
attributes = inspect.getmembers(
    Parameters, lambda attribute: not (inspect.isroutine(attribute)))
attribute_names = [attribute[0] for attribute in attributes
                   if not (attribute[0].startswith('__') and
                           attribute[0].endswith('__'))]

# set up pandas dataframe to store model input params
n_model_runs = len(param_list)
n_ts = np.sum(np.array(mp.N_outputs))
#if mp.exhumation_rate > 0:
#    n_ts = mp.exhumation_steps
n_rows = n_model_runs * n_ts

ind = np.arange(n_rows)
columns = ['model_run', 'model_error', 'timestep', 'runtime_yr'] + attribute_names
columns += ['surface_elevation',
            'max_surface_temperature', 'T_change_avg']

df = pd.DataFrame(index=ind, columns=columns)

if mp.analyse_borehole_temp is True:
    print 'loading temperature data'
    dft = pd.read_csv(mp.temperature_file)

today = datetime.datetime.now()
today_str = '%i-%i-%i' % (today.day, today.month,
                          today.year)

fn = 'model_params_and_results_%i_runs_%s_%s.csv' \
     % (len(param_list), mp.output_fn_adj, today_str)
fn_path_csv = os.path.join(output_folder, fn)

AHe_ages_surface_all = []
AHe_ages_surface_corr_all = []
AHe_xcoords_surface_all = []
AHe_data_file = None

if mp.model_AHe_samples is True:
    dfhs = pd.read_csv(mp.AHe_data_file)
    AHe_data_file = dfhs

for model_run, param_set in enumerate(param_list):

    print '-' * 20

    # reload default params
    Parameters = mp()

    # update default parameters in Parameter class
    for scenario_param_name, scenario_parameter in \
            zip(scenario_param_names, param_set):

        if scenario_parameter is not None:
            # find model parameter name to adjust
            model_param_name = scenario_param_name[:-2]

            print 'updating parameter %s from %s to %s' \
                  % (model_param_name,
                     str(getattr(Parameters, model_param_name)),
                     str(scenario_parameter))

            # update model parameter
            setattr(Parameters, model_param_name, scenario_parameter)

    print '-' * 20

    # store input parameters in dataframe
    attributes = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]
    #for a in attribute_dict:
    #    if a[0] in df.columns:
    #        if type(a[1]) is list or type(a[1]) is np.ndarray:
    #            df.loc[model_run, a[0]] = str(a[1])
    #        else:
    #            df.loc[model_run, a[0]] = a[1]

    print 'running single model'
    try:
        output = beo.model_run(Parameters)
    except Exception, msg:
        print '!' * 10
        print 'error running model run %i' % model_run
        print msg
        print '!' * 10

        for j in range(n_ts):

            output_number = model_run * n_ts + j

            #k = output_steps[j]

            for a in attribute_dict:
                if a[0] in df.columns:
                    if type(a[1]) is list or type(a[1]) is np.ndarray:
                        df.loc[output_number, a[0]] = str(a[1])
                    else:
                        df.loc[output_number, a[0]] = a[1]

            df.loc[output_number, 'model_error'] = str(msg)

        continue
        
    (runtimes, xyz_array, surface_levels, x_flt, z_flt,
     T_init_array, T_array, boiling_temp_array,
     xyz_array_exc, exceed_boiling_temp_array,
     xyz_element_array, qh_array, qv_array,
     fault_fluxes, durations,
     xzs, Tzs,
     Ahe_ages_all, Ahe_ages_corr_all, xs_Ahe_all,
     target_depths,
     AHe_ages_samples_all) = output

    # crop output to only the output timesteps, to limit filesize
    output_steps = []
    for duration, N_output in zip(mp.durations, mp.N_outputs):
        nt = int(duration / mp.dt)

        output_steps_i = list(np.linspace(0, nt-1, N_output).astype(int))
        output_steps += output_steps_i

    # select data for output steps only
    output_steps = np.array(output_steps)

    #if mp.exhumation_rate != 0:
    #    print 'exhumation, making sure output steps are equal to steps where ' \
    #          'a new surface level is reached'
    #    output_steps = [i for i, s in enumerate(surface_levels) if s in target_depths]

    n_ts = len(output_steps)

    Tzs_cropped = [Tzi[output_steps] for Tzi in Tzs]

    if Ahe_ages_all is not None:
        AHe_ages_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_all]
        AHe_ages_corr_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_corr_all]
    else:
        AHe_ages_cropped = None
        AHe_ages_corr_cropped = None

    if mp.calculate_he_ages is True and mp.model_AHe_samples is True:
        AHe_ages_samples_cropped = [AHe_i[output_steps]
                                    for AHe_i in AHe_ages_samples_all]
    else:
        AHe_ages_samples_cropped = None

    N_output_steps = len(output_steps)

    T_array = T_array[output_steps]

    #
    T_surface = []
    x_surface = []

    # add surface AHe data to output
    AHe_ages_surface = []
    AHe_ages_surface_corr = []
    AHe_xcoords_surface = []
    AHe_ages_samples_surface = []
    borehole_xlocs = None
    borehole_zlocs = None
    borehole_depths = None
    borehole_temp_measured = []
    borehole_temps_modeled = []

    if mp.analyse_borehole_temp is True:

        print 'extracting borehole temperature data'

        borehole_xlocs = np.zeros((len(mp.borehole_names), n_ts))
        borehole_zlocs = np.zeros_like(borehole_xlocs)
        borehole_depths = []
        borehole_temp_measured = []

        for borehole_number, borehole, xloc_raw in zip(itertools.count(),
                                                       mp.borehole_names,
                                                       mp.borehole_xs):

            ind = dft['borehole'] == borehole
            dfti = dft.loc[ind]

            borehole_depth = dfti['depth'].values
            T_obs = dfti['temperature'].values

            borehole_depths.append(borehole_depth)
            borehole_temp_measured.append(T_obs)

            borehole_temp_modeled = np.zeros((n_ts, len(T_obs)))

            for j in range(n_ts):

                z_surface = surface_levels[output_steps[j]]

                #ind_ts = np.where()
                surface_ind = np.where(z_flt < z_surface)[0][0]
                x_flt_step = x_flt[surface_ind]
                xloc = x_flt_step + xloc_raw

                borehole_xlocs[borehole_number, j] = xloc
                borehole_zlocs[borehole_number, j] = z_flt[surface_ind]

                output_number = model_run * n_ts + j


                col_name = 'modeled_xcoord_borehole_%s_run_%i_timestep_%i' \
                           % (borehole, model_run, output_steps[j])
                dft[col_name] = xloc

                z_obs = surface_levels[output_steps[j]] - borehole_depth

                x_buffer = 200.0
                pre_select = np.abs(xyz_array[:, 0] - xloc) < x_buffer

                zi, Ti = interpolate_data(xyz_array[pre_select],
                                          T_array[j][pre_select],
                                          xloc, y_steps=300)

                T_mod = np.interp(z_obs, zi, Ti)

                borehole_temp_modeled[j, :] = T_mod

                col_name = 'modeled_T_%s_run_%i_timestep_%i' \
                           % (borehole, model_run, output_steps[j])
                dft[col_name] = T_mod

                T_error = T_mod - T_obs
                df.loc[output_number, 'ME_temperature_%s' % borehole] = \
                    T_error.mean()
                df.loc[output_number, 'MAE_temperature_%s' % borehole] = \
                    np.mean(np.abs(T_error))

                RMSE_T = np.sqrt(np.mean(T_error ** 2))
                df.loc[output_number, 'RMSE_temperature_%s' % borehole] = \
                    RMSE_T

                R2_T = coefficient_of_determination(T_obs, T_mod)
                df.loc[output_number, 'R2_temperature_%s' % borehole] = R2_T

            #borehole_xlocs.append(borehole_xloc)
            borehole_temps_modeled.append( borehole_temp_modeled)

    for j in range(n_ts):

        output_number = model_run * n_ts + j

        #k = output_steps[j]

        for a in attribute_dict:
            if a[0] in df.columns:
                if type(a[1]) is list or type(a[1]) is np.ndarray:
                    df.loc[output_number, a[0]] = str(a[1])
                else:
                    df.loc[output_number, a[0]] = a[1]

        # store model results in dataframe
        df.loc[output_number, 'model_run'] = model_run
        df.loc[output_number, 'timestep'] = output_steps[j]
        df.loc[output_number, 'output_timestep'] = j
        df.loc[output_number, 'runtime_yr'] = runtimes[output_steps[j]] / year
        df.loc[output_number, 'surface_elevation'] = \
            surface_levels[output_steps[j]]

        T_change = T_array[j] - T_init_array
        df.loc[output_number, 'T_change_avg'] = T_change.mean()

        n_depths = len(target_depths)

        for i in range(n_depths):
            df.loc[output_number, 'max_temperature_layer_%i' % i] = Tzs_cropped[i][j].max()

        # add output T at surface
        surface_elev = surface_levels[output_steps[j]]

        if surface_elev in target_depths:
            surface_ind = np.where(target_depths == surface_elev)[0][0]
            T_surface_i = Tzs_cropped[surface_ind][j]
            x_coords_i = xzs[surface_ind]

        else:
            # interpolate T from nearest surfaces
            diff = target_depths - surface_elev
            ind_low = np.where(diff < 0)[0][-1]
            ind_high = np.where(diff > 0)[0][0]

            #fraction = np.abs(diff[ind_low]) / (target_depths[ind_high] - target_depths[ind_low])

            T_down = Tzs_cropped[ind_low][j]
            T_up = Tzs_cropped[ind_high][j]

            if len(T_down) != len(T_up):
                print 'warning, trying to interpolate two layers with unequal number of nodes'

            #T_surface_i = ((1.0-fraction) * Tzs_cropped[ind_low][j] + fraction * Tzs_cropped[ind_high][j])
            #x_coords_i = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]
            T_surface_i = Tzs_cropped[ind_low][j]
            x_coords_i = xzs[ind_low]

        T_surface.append(T_surface_i)
        x_surface.append(x_coords_i)

        df.loc[output_number, 'max_surface_temperature'] = T_surface_i.max()

        # calculate partial resetting and full resetting distance in the AHe data
        if Ahe_ages_all is not None:

            n_depths = len(AHe_ages_cropped)
            nt_output = AHe_ages_cropped[0].shape[0]

            for i in range(n_depths):
                ages_raw = AHe_ages_cropped[i][j] / My
                ages_raw_corr = AHe_ages_corr_cropped[i][j] / My
                x_step = 1.0
                x_coords_int = np.arange(xzs[i].min(), xzs[i].max() + x_step, x_step)
                ages = np.interp(x_coords_int, xzs[i], ages_raw)
                dev_age = ages / ages.max()

                min_age = np.min(ages)
                max_age = np.max(ages)
                ind_min_age = np.argmin(ages)

                col_name = 'elevation_layer%i' % i
                df.loc[output_number, col_name] = target_depths[i]

                x_min_age = x_coords_int[ind_min_age]
                col_name = 'lowest_age_layer%i' % i
                df.loc[output_number, col_name] = min_age
                col_name = 'highest_age_layer%i' % i
                df.loc[output_number, col_name] = max_age
                col_name = 'x_lowest_age_layer%i' % i
                df.loc[output_number, col_name] = x_min_age

                if dev_age.min() < mp.partial_reset_limit:
                    ind_partial = np.where(dev_age < mp.partial_reset_limit)[0]
                    x_partial_min = x_coords_int[ind_partial[0]]
                    x_partial_max = x_coords_int[ind_partial[-1]]

                    col_name = 'x_min_partial_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_partial_min
                    col_name = 'x_max_partial_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_partial_max
                else:
                    col_name = 'x_min_partial_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan
                    col_name = 'x_max_partial_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan

                if ages.min() < mp.reset_limit:
                    ind_full = np.where(ages < mp.reset_limit)[0]
                    x_full_min = x_coords_int[ind_full[0]]
                    x_full_max = x_coords_int[ind_full[-1]]
                    col_name = 'x_min_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_full_min
                    col_name = 'x_max_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_full_max
                else:
                    col_name = 'x_min_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan
                    col_name = 'x_max_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan

            # figure out which depth is currently at the surface
            # and calculate the partial and full reset widths for these
            if surface_elev in target_depths:
                surface_ind = np.where(target_depths == surface_elev)[0][0]
                ages_raw = AHe_ages_cropped[surface_ind][j] / My
                ages_raw_corr = AHe_ages_corr_cropped[surface_ind][j] / My
                x_coords = xzs[surface_ind]

            else:
                # interpolate AHe age from nearest surfaces
                diff = target_depths - surface_elev
                ind_low = np.where(diff < 0)[0][-1]
                ind_high = np.where(diff > 0)[0][0]

                fraction = np.abs(diff[ind_low]) \
                           / (target_depths[ind_high] - target_depths[ind_low])

                #ages_raw = ((1.0-fraction) * AHe_ages_cropped[ind_low][j]
                #            + fraction * AHe_ages_cropped[ind_high][j]) / My
                #x_coords = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]
                ages_raw = AHe_ages_cropped[ind_low][j] / My
                ages_raw_corr = AHe_ages_corr_cropped[ind_low][j] / My
                x_coords = xzs[ind_low]

            # add surface AHe data to output
            AHe_ages_surface.append(ages_raw * My)
            AHe_ages_surface_corr.append(ages_raw_corr * My)
            AHe_xcoords_surface.append(x_coords)

            x_coords_int = np.arange(x_coords.min(), x_coords.max() + x_step, x_step)

            if len(ages_raw) != len(x_coords):
                print 'warning, lenght of AHe ages array and x coordinates ' \
                      'array are not equal. Skipping calculating partial/full ' \
                      'reset zone for surface layer at timestep %i of %i' \
                      % (j, n_ts)
                continue

            ages = np.interp(x_coords_int, x_coords, ages_raw)
            dev_age = ages / ages.max()

            min_age = np.min(ages)
            max_age = np.max(ages)
            ind_min_age = np.argmin(ages)

            x_min_age = x_coords_int[ind_min_age]
            col_name = 'lowest_age_surface'
            df.loc[output_number, col_name] = min_age
            col_name = 'highest_age_surface'
            df.loc[output_number, col_name] = max_age
            col_name = 'x_lowest_age_surface'
            df.loc[output_number, col_name] = x_min_age

            if dev_age.min() < mp.partial_reset_limit:
                ind_partial = np.where(dev_age < mp.partial_reset_limit)[0]
                x_partial_min = x_coords_int[ind_partial[0]]
                x_partial_max = x_coords_int[ind_partial[-1]]

                col_name = 'x_min_partial_reset_surface'
                df.loc[output_number, col_name] = x_partial_min
                col_name = 'x_max_partial_reset_surface'
                df.loc[output_number, col_name] = x_partial_max
            else:
                col_name = 'x_min_partial_reset_surface'
                df.loc[output_number, col_name] = np.nan
                col_name = 'x_max_partial_reset_surface'
                df.loc[output_number, col_name] = np.nan

            if ages.min() < mp.reset_limit:
                ind_full = np.where(ages < mp.reset_limit)[0]
                x_full_min = x_coords_int[ind_full[0]]
                x_full_max = x_coords_int[ind_full[-1]]
                col_name = 'x_min_full_reset_surface'
                df.loc[output_number, col_name] = x_full_min
                col_name = 'x_max_full_reset_surface'
                df.loc[output_number, col_name] = x_full_max
            else:
                col_name = 'x_min_full_reset_surface'
                df.loc[output_number, col_name] = np.nan
                col_name = 'x_max_full_reset_surface'
                df.loc[output_number, col_name] = np.nan

            AHe_ages_surface_all.append(AHe_ages_surface)
            AHe_ages_surface_corr_all.append(AHe_ages_surface_corr)
            AHe_xcoords_surface_all.append(AHe_xcoords_surface)

        AHe_ages_samples_surface = []
        if mp.calculate_he_ages and mp.model_AHe_samples is True:

            for i in range(N_output_steps):
                surface_elev = surface_levels[i]

                if surface_elev in target_depths:
                    surface_ind = np.where(target_depths == surface_elev)[0][0]
                    ages_raw = AHe_ages_samples_cropped[surface_ind][i]
                    #x_coords = xzs[surface_ind]

                else:
                    # interpolate AHe age from nearest surfaces
                    diff = target_depths - surface_elev
                    ind_low = np.where(diff < 0)[0][-1]
                    ind_high = np.where(diff > 0)[0][0]

                    fraction = np.abs(diff[ind_low]) / (target_depths[ind_high]
                                                        - target_depths[ind_low])

                    ages_raw = ((1.0-fraction) * AHe_ages_samples_cropped[ind_low][i]
                                + fraction * AHe_ages_samples_cropped[ind_high][i])

                    #x_coords = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]

                # add surface AHe data to output
                AHe_ages_samples_surface.append(ages_raw)

    # analyze model-data fit of AHe surface samples
    if mp.model_AHe_samples is True:

        print 'analyzing fit of modeled and measured AHe ages'

        n_grains = len(AHe_ages_samples_surface[0])

        me_ahe = np.zeros(N_output)
        mae_ahe = np.zeros(N_output)
        mswd_ahe = np.zeros(N_output)

        for timestep in range(N_output):

            output_number2 = model_run * n_ts + timestep
            profile_loc = dfhs['profile'] == mp.profile_number
            dfhs2 = dfhs.loc[profile_loc]

            AHe_data_file = dfhs2

            diff = dfhs2['AHe_age_uncorr'].values - AHe_ages_samples_surface[timestep] / My

            me_ahe[timestep] = np.mean(diff)
            mae_ahe[timestep] = np.mean(np.abs(diff))
            mswd_ahe[timestep] = np.sum((diff / (0.5 * dfhs2['AHe_age_uncorr_2se'])) ** 2) / (n_grains - 1)

            df.loc[output_number2, 'mean_error_AHe_samples'] = me_ahe[timestep]
            df.loc[output_number2, 'mean_abs_error_AHe_samples'] = mae_ahe[timestep]
            df.loc[output_number2, 'mswd_AHe_samples'] = mswd_ahe[timestep]

    # option to save corrected ages for figure output


    # gather and save model output
    output_selected = \
        [runtimes, runtimes[output_steps], xyz_array,
         surface_levels[output_steps],
         T_init_array,
         T_array, boiling_temp_array[output_steps],
         xyz_array_exc, exceed_boiling_temp_array[output_steps],
         xyz_element_array,
         qh_array[output_steps], qv_array[output_steps],
         fault_fluxes, durations,
         xzs, Tzs_cropped, x_surface, T_surface,
         AHe_ages_cropped, AHe_ages_corr_cropped, xs_Ahe_all, target_depths,
         AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
         AHe_ages_samples_surface, AHe_data_file,
         borehole_xlocs, borehole_zlocs,
         borehole_depths, borehole_temp_measured, borehole_temps_modeled]

#    borehole_xlocs = None
#    borehole_zlocs = None
#    borehole_depths = None
#    borehole_temp_measured = []
#    borehole_temps_modeled = []


    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month,
                              today.year)

    fn = 'results_model_run_%i_%s_%s_%s.pck' \
         % (model_run, str(param_set), mp.output_fn_adj, today_str)
    fn_path = os.path.join(output_folder, fn)

    print 'saving model results as %s' % fn_path
    fout = open(fn_path, 'w')
    pickle.dump(output_selected, fout)
    fout.close()

    print '-' * 30
    print 'saving summary of parameters and model results as %s' % fn_path_csv
    df.to_csv(fn_path_csv, index_label='row', encoding='utf-8')

    if mp.analyse_borehole_temp is True:
        fn_new = os.path.split(mp.temperature_file)[-1].split('.')[:-1]
        fn_new = ''.join(fn_new)
        fn_new += '_modeled_%i_runs_%s_%s.csv' % (len(param_list),
                                                  mp.output_fn_adj,
                                                  today_str)
        print 'saving modeled temperatures for boreholes to %s' % fn_new
        dft.to_csv(os.path.join(output_folder, fn_new))

    if mp.save_AHe_ages is True and AHe_ages_surface_all != []:

        nxs = np.max(np.array([AHe_ii.shape[0]
                               for AHe_i in AHe_ages_surface_all
                               for AHe_ii in AHe_i]))
        nts = len(AHe_ages_surface)

        cols = []
        for j in range(n_model_runs):
            for i in range(nts):
                cols += ['x_run_%i_ts%i' % (j, i),
                         'AHe_age_uncorr_run_%i_ts%i' % (j, i)]

        dfh = pd.DataFrame(columns=cols, index=np.arange(nxs))

        for j in range(n_model_runs):
            for i in range(nts):
                a = len(dfh)
                b = len(AHe_xcoords_surface_all[j][i])

                l = np.min((a, b))
                try:
                    dfh.loc[:l-1, 'x_run_%i_ts%i' % (j, i)] = \
                        AHe_xcoords_surface_all[j][i][:l]
                    dfh.loc[:l-1, 'AHe_age_uncorr_run_%i_ts%i' % (j, i)] = \
                        AHe_ages_surface_all[j][i][:l] / My
                except Exception, msg:
                    print 'error, something went wrong with saving AHe data ' \
                          'to a .csv file for model run %i and timestep %i' \
                          % (j, i)
                    print msg
                    print 'continuing with next timestep'

        fnh = 'AHe_surface_modeled_%i_runs_%s_%s.csv' % (n_model_runs,
                                                         mp.output_fn_adj,
                                                         today_str)

        print 'saving modeled AHe ages at the surface to %s' % fnh
        dfh.to_csv(os.path.join(output_folder, fnh), index_label='row')

    if mp.model_AHe_samples is True:

        # save AHe ages
        for timestep in range(N_output):

            # and store in datafile
            col_name = 'modeled_AHe_age_uncorr_run_%i_timestep_%i' \
                       % (model_run, output_steps[timestep])

            profile_loc = dfhs['profile'] == mp.profile_number

            if True in profile_loc:
                try:
                    dfhs.loc[profile_loc, col_name] = AHe_ages_samples_surface[timestep] / My
                except Exception, msg:
                    print 'error, something went wrong with saving AHe ' \
                          'sample data to a .csv file for model run %i ' \
                          'and timestep %i' % (j, i)
                    print msg
                    print 'continuing with next timestep'

        output_fn1 = os.path.split(mp.AHe_data_file)[-1]
        output_fn2 = output_fn1[:-4] + '_modeled_%i_runs_%s_%s.csv' \
                                       % (n_model_runs,
                                          mp.output_fn_adj,
                                          today_str)
        output_fn = os.path.join(output_folder, output_fn2)
        print 'saving modeled AHe ages samples to %s' % output_fn
        dfhs.to_csv(output_fn, index=False)

print '-' * 30
print 'done with all model runs'
