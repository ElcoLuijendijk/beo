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

import numpy as np
import pandas as pd

import model_parameters.model_parameters as mp
import model_parameters.parameter_ranges as pr
import beo

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
Parameters = mp

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
columns = ['timestep', 'runtime_yr'] + attribute_names
columns += ['surface_elevation',
            'max_surface_temperature', 'T_change_avg']

df = pd.DataFrame(index=ind, columns=columns)

for model_run, param_set in enumerate(param_list):

    Parameters = mp

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
    for a in attribute_dict:
        if a[0] in df.columns:
            if type(a[1]) is list or type(a[1]) is np.ndarray:
                df.loc[model_run, a[0]] = str(a[1])
            else:
                df.loc[model_run, a[0]] = a[1]

    print 'running single model'
    output = beo.model_run(mp)

    (runtimes, xyz_array, surface_levels,
     T_init_array, T_array, xyz_element_array,
     qh_array, qv_array,
     fault_fluxes, durations, xzs, Tzs, Ahe_ages_all, xs_Ahe_all,
     z_Ahe) = output

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
    output_selected = \
        [runtimes, runtimes[output_steps], xyz_array,
         surface_levels[output_steps],
         T_init_array,
         T_array[output_steps], xyz_element_array,
         qh_array[output_steps], qv_array[output_steps],
         fault_fluxes, durations, xzs, Tzs_cropped,
         AHe_ages_cropped, xs_Ahe_all, z_Ahe]

    T_array = T_array[output_steps]

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

        df.loc[output_number, 'max_surface_temperature'] = Tzs[0][j].max()
        T_change = T_array[j] - T_init_array
        df.loc[output_number, 'T_change_avg'] = T_change.mean()

        # calculate partial resetting and full resetting distance in the AHe data
        if Ahe_ages_all is not None:

            n_depths = len(AHe_ages_cropped)
            nt_output = AHe_ages_cropped[0].shape[0]

            for i in range(n_depths):
                ages = AHe_ages_cropped[i][j] / My
                dev_age = ages / ages.max()

                min_age = np.min(ages)
                max_age = np.max(ages)
                ind_min_age = np.argmin(ages)

                col_name = 'elevation_layer%i' % i
                df.loc[output_number, col_name] = z_Ahe[i]

                x_min_age = xzs[i][ind_min_age]
                col_name = 'lowest_age_layer%i' % i
                df.loc[output_number, col_name] = min_age
                col_name = 'highest_age_layer%i' % i
                df.loc[output_number, col_name] = max_age
                col_name = 'x_lowest_age_layer%i' % i
                df.loc[output_number, col_name] = x_min_age

                if dev_age.min() < mp.partial_reset_limit:
                    ind_partial = np.where(dev_age < mp.partial_reset_limit)[0]
                    x_partial_min = xzs[i][ind_partial[0]]
                    x_partial_max = xzs[i][ind_partial[-1]]

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
                    x_full_min = xzs[i][ind_full[0]]
                    x_full_max = xzs[i][ind_full[-1]]
                    col_name = 'x_min_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_full_min
                    col_name = 'x_max_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = x_full_max
                else:
                    col_name = 'x_min_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan
                    col_name = 'x_max_full_reset_layer%i' % i
                    df.loc[output_number, col_name] = np.nan

    fn = 'T_field_model_run_%i_%s.pck' \
         % (model_run, str(param_set))
    fn_path = os.path.join(output_folder, fn)

    print 'saving model results as %s' % fn_path
    fout = open(fn_path, 'w')
    pickle.dump(output_selected, fout)
    fout.close()

fn = 'model_params_and_results_%i_runs.csv' \
     % len(param_list)
fn_path = os.path.join(output_folder, fn)

print '-' * 30
print 'saving summary of parameters and model results as %s' % fn_path

df.to_csv(fn_path, index_label='model_run', encoding='utf-8')

print 'done with all model runs'
