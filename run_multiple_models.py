"""
run multiple model sequentially with different parameter values

"""

__author__ = 'elcopone'

import os
import pickle
import itertools
import inspect

import numpy as np
import pandas as pd

import model_parameters.model_parameters as mp
import model_parameters.parameter_ranges as pr
import hydrotherm_escript

day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6

# create list with param values for each model run
param_list = \
    list(itertools.product(pr.fault_bottoms,
                           pr.thermal_gradients))

# read default model parameter file
Parameters = mp

# get attributes
attributes = inspect.getmembers(
    Parameters, lambda attribute: not (inspect.isroutine(attribute)))
attribute_names = [attribute[0] for attribute in attributes
                   if not (attribute[0].startswith('__') and
                           attribute[0].endswith('__'))]

# set up pandas dataframe to store model input params
ind = np.arange(len(param_list))
columns = attribute_names
columns += ['max_surface_temperature', 'T_change_avg']

df = pd.DataFrame(index=ind, columns=columns)


for model_run, param_set in enumerate(param_list):

    fault_bottom, thermal_gradient = param_set
    print 'updated parameters ', param_set

    # update parameters in param file
    mp.fault_bottoms[0] = fault_bottom
    mp.thermal_gradient = thermal_gradient

    # store input parameters in dataframe
    Parameters = mp
    attributes = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]
    for a in attribute_dict:
        if a[0] in df.columns:
            if type(a[1]) is list:
                df.loc[model_run, a[0]] = str(a[1])
            else:
                df.loc[model_run, a[0]] = a[1]

    print 'running single model'
    output = hydrotherm_escript.model_run(mp)

    runtimes, xyz_array, T_init_array, T_array, xyz_element_array, qh_array, qv_array, \
          fault_fluxes, durations, xzs, Tzs, Ahe_ages_all, xs_Ahe_all = output

#    runtimes, xyz_array, T_init_array, T_array, xyz_element_array, qh_array, qv_array, \
#          fault_fluxes, durations, xzs, Tzs, Ahe_ages_all, xs_Ahe_all = output

    # calculate partial resetting and full resetting distance in the AHe data
    if Ahe_ages_all is not None:

        n_depths = len(Ahe_ages_all)
        nt_output = Ahe_ages_all[0].shape[0]

        for i in range(n_depths):
            for j in range(nt_output):
                ages = Ahe_ages_all[i][j] / My
                dev_age = ages / ages.max()

    # store model results in dataframe
    df.loc[model_run, 'max_surface_temperature'] = Tzs[0].max()
    T_change = T_array[-1] - T_init_array
    df.loc[model_run, 'T_change_avg'] = T_change.mean()

    #output_folder = os.path.join(folder, 'model_output')
    output_folder = mp.output_folder
    fn = 'T_field_model_run_%i_%s.pck' \
         % (model_run, str(param_set))
    fn_path = os.path.join(output_folder, fn)

    print 'saving model results as %s' % fn_path

    fout = open(fn_path, 'w')
    pickle.dump(output, fout)
    fout.close()

fn = 'model_params_and_results_%i_runs.csv' \
     % len(param_list)
fn_path = os.path.join(output_folder, fn)

print 'saving summary of parameters and model results as %s' % fn_path

df.to_csv(fn_path, index_label='model_run', encoding='utf-8')

print 'done with all model runs'
