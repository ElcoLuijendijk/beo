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

import numpy as np
import pandas as pd

from model_parameters.model_parameters import ModelParams
import model_parameters.parameter_ranges as pr
import beo

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
n_rows = n_model_runs * n_ts

ind = np.arange(n_rows)
columns = ['model_run', 'timestep', 'runtime_yr'] + attribute_names
columns += ['surface_elevation',
            'max_surface_temperature', 'T_change_avg']

df = pd.DataFrame(index=ind, columns=columns)

for model_run, param_set in enumerate(param_list):

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
    output = beo.model_run(Parameters)

    (runtimes, xyz_array, surface_levels,
     T_init_array, T_array, boiling_temp_array,
     xyz_array_exc, exceed_boiling_temp_array,
     xyz_element_array, qh_array, qv_array,
     fault_fluxes, durations,
     xzs, Tzs,
     Ahe_ages_all, xs_Ahe_all,
     target_depths) = output

    # crop output to only the output timesteps, to limit filesize
    output_steps = []
    for duration, N_output in zip(mp.durations, mp.N_outputs):
        nt = int(duration / mp.dt)

        output_steps_i = list(np.linspace(0, nt-1, N_output).astype(int))
        output_steps += output_steps_i

    # select data for output steps only
    output_steps = np.array(output_steps)

    Tzs_cropped = [Tzi[output_steps] for Tzi in Tzs]
    AHe_ages_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_all]

    T_array = T_array[output_steps]

    #
    T_surface = []
    x_surface = []

    # add surface AHe data to output
    AHe_ages_surface = []
    AHe_xcoords_surface = []

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
            # interpolate AHe age from nearest surfaces
            diff = target_depths - surface_elev
            ind_low = np.where(diff < 0)[0][-1]
            ind_high = np.where(diff > 0)[0][0]

            fraction = np.abs(diff[ind_low]) / (target_depths[ind_high] - target_depths[ind_low])

            T_surface_i = ((1.0-fraction) * Tzs_cropped[ind_low][j] + fraction * Tzs_cropped[ind_high][j])

            x_coords_i = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]

        T_surface.append(T_surface_i)
        x_surface.append(x_coords_i)

        df.loc[output_number, 'max_surface_temperature'] = T_surface_i.max()

        # calculate partial resetting and full resetting distance in the AHe data
        if Ahe_ages_all is not None:

            n_depths = len(AHe_ages_cropped)
            nt_output = AHe_ages_cropped[0].shape[0]

            for i in range(n_depths):
                ages_raw = AHe_ages_cropped[i][j] / My
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
                surface_ind = np.where(target_depths == surface_elev)[0]
                ages_raw = AHe_ages_cropped[surface_ind][j] / My
                x_coords = xzs[surface_ind]

            else:
                # interpolate AHe age from nearest surfaces
                diff = target_depths - surface_elev
                ind_low = np.where(diff < 0)[0][-1]
                ind_high = np.where(diff > 0)[0][0]

                fraction = np.abs(diff[ind_low]) \
                           / (target_depths[ind_high] - target_depths[ind_low])

                ages_raw = ((1.0-fraction) * AHe_ages_cropped[ind_low][j]
                            + fraction * AHe_ages_cropped[ind_high][j]) / My

                x_coords = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]

            # add surface AHe data to output
            AHe_ages_surface.append(ages_raw * My)
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
         AHe_ages_cropped, xs_Ahe_all, target_depths,
         AHe_ages_surface, AHe_xcoords_surface]

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month,
                              today.year)

    fn = 'results_model_run_%i_%s_%s.pck' \
         % (model_run, str(param_set), today_str)
    fn_path = os.path.join(output_folder, fn)

    print 'saving model results as %s' % fn_path
    fout = open(fn_path, 'w')
    pickle.dump(output_selected, fout)
    fout.close()

today = datetime.datetime.now()
today_str = '%i-%i-%i' % (today.day, today.month,
                          today.year)

fn = 'model_params_and_results_%i_runs_%s.csv' \
     % (len(param_list), today_str)
fn_path = os.path.join(output_folder, fn)

print '-' * 30
print 'saving summary of parameters and model results as %s' % fn_path

df.to_csv(fn_path, index_label='row', encoding='utf-8')

print 'done with all model runs'
