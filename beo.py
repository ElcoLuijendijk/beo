"""
run multiple models sequentially with different parameter values

"""

__author__ = 'Elco Luijendijk'

import sys
import os
import re
import glob
import pickle
import difflib
import itertools
import inspect
import datetime
import time
#import imp

import pdb

import scipy.interpolate

import numpy as np
import pandas as pd

import beo_core

# escript is optional, Beo can also be run with the fipy backend
try:
    import esys.escript as es
except ImportError:
    es = None


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


def find_referenced_parameter_names(source_dir):

    """
    Scan the Beo source code for every parameter name that is read from the
    model parameter object, either directly (mp.name) or through getattr
    (getattr(mp, 'name', default)).

    Used by warn_about_possible_parameter_typos to check the parameter names
    defined in a model parameter file against the names that the code
    actually reads.
    """

    pattern = re.compile(r"\bmp\.(\w+)\b|getattr\(\s*mp\s*,\s*['\"](\w+)['\"]")

    source_files = glob.glob(os.path.join(source_dir, '*.py'))
    source_files += glob.glob(os.path.join(source_dir, 'lib', '*.py'))

    referenced = set()
    for source_file in source_files:
        with open(source_file) as fin:
            text = fin.read()
        for match in pattern.finditer(text):
            referenced.add(match.group(1) or match.group(2))

    return referenced


def warn_about_possible_parameter_typos(parameter_names, source_dir):

    """
    Warn about parameter names that are defined in the model parameter file
    but are not read anywhere in the Beo source code, if the name is close to
    one that is read.

    A required parameter that is misspelled in the model parameter file
    raises an AttributeError as soon as the code tries to read it, which is
    easy to notice. A parameter that is only read with
    getattr(mp, name, default) instead silently falls back to the default
    value if the name is misspelled, without the model run failing or
    printing anything. This check is meant to catch that case, and cannot
    catch a misspelled parameter name for which no closely named parameter
    exists.
    """

    referenced = find_referenced_parameter_names(source_dir)

    for name in parameter_names:

        if name in referenced:
            continue

        matches = difflib.get_close_matches(name, referenced, n=1, cutoff=0.8)

        if matches:
            print('warning, the model parameter file defines "%s", which is '
                  'not read anywhere in the Beo source code. did you mean '
                  '"%s"?' % (name, matches[0]))


day = 24.0 * 60.0 * 60.0
year = 365.25 * day
My = year * 1e6

scriptdir = os.path.realpath(sys.path[0])

if len(sys.argv) > 1:

    inp_file_loc = os.path.join(scriptdir, sys.argv[-1])

    print('model input files: ', inp_file_loc)

    try:
        #model_parameters = imp.load_source('model_parameters', inp_file_loc)

        import importlib.util
        spec = importlib.util.spec_from_file_location('model_parameters', inp_file_loc)
        model_parameters = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(model_parameters)
    except IOError:
        msg = 'cannot find parameter file %s' % inp_file_loc
        raise IOError(msg)

    ModelParams = model_parameters.ModelParams
    pr = model_parameters.ParameterRanges

else:

    print('running model input data from file ' \
          'model_parameters/model_parameters.py')

    from model_parameters.model_parameters import ModelParams
    from model_parameters.model_parameters import ParameterRanges as pr

mp = ModelParams

output_folder = os.path.join(scriptdir, mp.output_folder)

if os.path.exists(output_folder) is False:
    os.mkdir(output_folder)
    print('created new directory for model output: %s' % output_folder)

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

warn_about_possible_parameter_typos(attribute_names, scriptdir)

# set up pandas dataframe to store model input params
n_model_runs = len(param_list)

if getattr(mp, 'dt_growth_factor', None) is not None and getattr(mp, 'stored_timesteps', None) is None:
    msg = 'error, stored_timesteps must be set when dt_growth_factor is used.\n'
    msg += 'the regular dt_stored output interval assumes a constant dt, which does\n'
    msg += 'not hold while the timestep is still growing, so an explicit list of\n'
    msg += 'output times is required instead'
    raise ValueError(msg)

if getattr(mp, 'stored_timesteps', None) is not None:
    if getattr(mp, 'dt_growth_factor', None) is not None:
        if beo_core.beo_fipy is None:
            msg = 'error, dt_growth_factor requires the fipy backend, but fipy and the '
            msg += 'gmsh python module could not be imported.\n'
            msg += 'install these with pip install fipy gmsh'
            raise ImportError(msg)

        # replay the same growing, then constant, timestep schedule that
        # the fipy backend will use, see build_growing_dt_schedule and
        # OutputTimeSelector in lib/beo_fipy.py, so that the exact
        # cumulative model time at every solve timestep is known here as
        # well, without having to run the model first
        all_dts = []
        for duration in mp.durations:
            all_dts += beo_core.beo_fipy.build_growing_dt_schedule(
                mp.dt, mp.dt_growth_factor, mp.dt_max, duration)
        time_grid = np.cumsum(np.array(all_dts))
    else:
        # a requested time that is closer to the previous one than dt is
        # stored at the same timestep as the previous one, see
        # OutputTimeSelector in lib/beo_fipy.py. snapping every requested
        # time up to the nearest multiple of dt the same way here gives the
        # exact number of distinct timesteps that will be stored, so that
        # the dataframe is preallocated with exactly the right number of
        # rows
        nt_total = int(np.sum(mp.durations) / mp.dt) + 1
        time_grid = mp.dt * np.arange(1, nt_total + 1)

    snapped_times = []
    for target in mp.stored_timesteps:
        reachable = time_grid[time_grid >= target]
        if len(reachable) > 0:
            snapped_times.append(reachable[0])
    n_ts = len(np.unique(snapped_times))
else:
    n_ts = np.sum(np.array(mp.N_outputs))

n_rows = n_model_runs * n_ts

ind = np.arange(n_rows)
columns = ['model_run', 'model_error', 'timestep', 'runtime_yr', 'computational_time'] + attribute_names
columns += ['surface_elevation',
            'max_surface_temperature', 'T_change_avg']

df = pd.DataFrame(index=ind, columns=columns)

if mp.analyse_borehole_temp is True:
    print('loading temperature data')
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

    print('-' * 20)

    print(f"model run {model_run}")

    if pr.initial_base_run is True and model_run == 0:
        print("using base-case parameter values for the first run")

    # reload default params
    Parameters = mp()

    # update default parameters in Parameter class
    for scenario_param_name, scenario_parameter in \
            zip(scenario_param_names, param_set):

        if scenario_parameter is not None:
            # find model parameter name to adjust
            model_param_name = scenario_param_name[:-2]

            print('updating parameter %s from %s to %s' \
                  % (model_param_name,
                     str(getattr(Parameters, model_param_name)),
                     str(scenario_parameter)))

            # update model parameter
            setattr(Parameters, model_param_name, scenario_parameter)

    print('-' * 20)

    # store input parameters in dataframe
    attributes = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]

    print('running single model')

    if mp.steady_state is True and mp.add_exhumation is True:
        msg = 'Error, both steady-state and exhumation are set to True. Please change your model parameters file'
        raise ValueError(msg)

    if (getattr(mp, 'stored_timesteps', None) is not None
            and mp.repeat_timeslices is not None and mp.repeat_timeslices > 1):
        msg = 'error, stored_timesteps and repeat_timeslices cannot be used together. '
        msg += 'stored_timesteps is a single list of absolute times covering the entire '
        msg += 'model run, which cannot be repeated in the same way as dt_stored'
        raise ValueError(msg)

    # check if timeslices should be repeated and if yes, extend durations, fault fluxes
    # and aquifer fluxes arrays
    if mp.repeat_timeslices is not None and mp.repeat_timeslices > 1:

        try:
            assert (type(mp.durations) is list
                    and type(mp.fault_fluxes) is list
                    and type(mp.aquifer_fluxes) is list
                    and type(mp.N_outputs) is list)
        except AssertionError:
            msg = 'error, make sure the parameters durations, fault_fluxes, ' \
                  'aquifer_fluxes and N_outputs are lists'
            raise ValueError(msg)

        mp.durations = mp.durations * mp.repeat_timeslices
        mp.fault_fluxes = mp.fault_fluxes * mp.repeat_timeslices
        mp.aquifer_fluxes = mp.aquifer_fluxes * mp.repeat_timeslices
        mp.N_outputs = mp.N_outputs * mp.repeat_timeslices

    no_exceptions = True
    if no_exceptions is True:
        start_time = time.time()
        output = beo_core.model_run(Parameters)
        comp_time_model_run = time.time() - start_time

        for j in range(n_ts):

            output_number = model_run * n_ts + j
            df.loc[output_number, 'computational_time'] = comp_time_model_run

    else:
        try:
            start_time = time.time()
            output = beo_core.model_run(Parameters)

            comp_time_model_run = time.time() - start_time

            for j in range(n_ts):

                output_number = model_run * n_ts + j
                df.loc[output_number, 'computational_time'] = comp_time_model_run

        except Exception as msg:
            print('!' * 10)
            print('error running model run %i' % model_run)
            print(msg)
            print('!' * 10)

            for j in range(n_ts):

                output_number = model_run * n_ts + j

                for a in attribute_dict:
                    if a[0] in df.columns:
                        if type(a[1]) is list or type(a[1]) is np.ndarray:
                            df.loc[output_number, a[0]] = str(a[1])
                        else:
                            df.loc[output_number, a[0]] = a[1]

                df.loc[output_number, 'model_error'] = str(msg)

            continue

    (runtimes, xyz_array, surface_levels, x_flt, z_flt,
     Ts, q_vectors, T_init_array, T_array, boiling_temp_array,
     xyz_array_exc, exceed_boiling_temp_array,
     xyz_element_array, qh_array, qv_array,
     fault_fluxes, durations,
     xzs, Tzs, Tzs_diff,
     Ahe_ages_surface_all, Ahe_ages_surface_corr_all, xs_Ahe_surface_all,
     target_depths,
     AHe_ages_surface_samples_all, AHe_ages_surface_samples_all_corr,
     z_borehole, AHe_ages_borehole, AHe_ages_borehole_corrected) = output

    if getattr(mp, 'stored_timesteps', None) is not None:
        # all stored timesteps are already exactly the times requested in
        # stored_timesteps, so none of the further selection below is needed
        output_steps = np.arange(len(runtimes))

        print('selecting output steps: ', output_steps)
        print('generating time output at steps: ')
        print(np.array(runtimes)[output_steps] / year)

    else:
        output_steps = [0]

        for duration, N_output in zip(mp.durations, mp.N_outputs):
            nt = int(duration / mp.dt_stored)
            print('timesteps = %i' % nt)

            output_steps_i = list(np.linspace(0, nt, N_output + 1).astype(int) + output_steps[-1])[1:]
            output_steps += output_steps_i

        # select data for output steps only
        output_steps = np.array(output_steps)

        print('selecting output steps: ', output_steps)

        times_test = np.arange(0, np.sum(mp.durations) + mp.dt_stored, mp.dt_stored)

        print('generating time output at steps: ')
        print(times_test[output_steps] / year)

    if mp.steady_state is False:
        n_ts = len(output_steps)
    else:
        n_ts = 2
        output_steps = np.array([0, -1])

    Tzs_cropped = [Tzi[output_steps] for Tzi in Tzs]
    Tzs_diff_cropped = [Tzi[output_steps] for Tzi in Tzs_diff]

    if Ahe_ages_surface_all is not None:
        AHe_ages_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_surface_all]
        AHe_ages_corr_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_surface_corr_all]
    else:
        AHe_ages_cropped = None
        AHe_ages_corr_cropped = None

    if AHe_ages_surface_samples_all is not None:
        AHe_ages_samples_cropped = [AHe_i[output_steps]
                                    for AHe_i in AHe_ages_surface_samples_all]
        AHe_ages_samples_corr_cropped = [AHe_i[output_steps]
                                         for AHe_i in AHe_ages_surface_samples_all_corr]
    else:
        AHe_ages_samples_cropped = None
        AHe_ages_samples_corr_cropped = None

    N_output_steps = len(output_steps)

    T_array = T_array[output_steps]

    save_VTK_file = hasattr(mp, "save_VTK_file") is True and mp.save_VTK_file is True
    backend = getattr(mp, 'backend', 'escript')

    if save_VTK_file is True and (es is None or backend == 'fipy'):
        print('no output to VTK file, saving VTK files is only supported for the escript backend')

    elif save_VTK_file is True:
        VTK_dir = 'VTK_files_model_run_%i_%s_%s_%s' \
                  % (model_run, str(param_set), mp.output_fn_adj, today_str)
        VTK_dir_full = os.path.join(output_folder, VTK_dir)
        print('saving VTK file of temperatures and flux in directory %s' % VTK_dir_full)

        if os.path.exists(VTK_dir_full) is False:
            os.mkdir(VTK_dir_full)

        VTK_data = es.DataManager(formats=[es.DataManager.VTK], work_dir=VTK_dir_full)

        for output_step in output_steps:

            VTK_data.addData(temperature=Ts[output_step], q=q_vectors[output_step])
            VTK_data.setTime(runtimes[output_step] / year)
            VTK_data.export()

    else:
        print('no output to VTK file. add save_VTK_file=True to input file to change this')

    #
    T_surface = []
    x_surface = []

    # add surface AHe data to output
    AHe_ages_surface = []
    AHe_ages_surface_corr = []
    AHe_xcoords_surface = []
    AHe_ages_samples_surface = []
    AHe_ages_samples_surface_corr = []
    borehole_xlocs = None
    borehole_zlocs = None
    borehole_depths = None
    borehole_temp_measured = []
    borehole_temps_modeled = []

    # find x coord of first fault, for relative positioning of borehole temperature data and AHe samples
    x_loc_fault = np.zeros(n_ts)
    z_loc_fault = np.zeros(n_ts)

    for j in range(n_ts):
        z_surface = surface_levels[output_steps[j]]
        surface_ind = np.where(z_flt <= z_surface)[0][0]
        x_flt_step = x_flt[surface_ind]
        x_loc_fault[j] = x_flt_step
        z_loc_fault[j] = z_flt[surface_ind]

    if mp.analyse_borehole_temp is True:

        print('extracting borehole temperature data')

        borehole_xlocs = np.zeros((len(mp.borehole_names), n_ts))
        borehole_zlocs = np.zeros_like(borehole_xlocs)
        borehole_depths = []
        borehole_temp_measured = []

        for borehole_number, borehole, xloc_raw in zip(itertools.count(),
                                                       mp.borehole_names,
                                                       mp.borehole_xs):

            if borehole not in dft['borehole'].values:
                msg = 'error, could not find borehole %s in the borehole temperature data file %s' \
                      % (borehole, mp.temperature_file)
                raise IndexError(msg)

            ind = dft['borehole'] == borehole
            dfti = dft.loc[ind]

            borehole_depth = dfti['depth'].values
            T_obs = dfti['temperature'].values

            borehole_depths.append(borehole_depth)
            borehole_temp_measured.append(T_obs)

            borehole_temp_modeled = np.zeros((n_ts, len(T_obs)))

            for j in range(n_ts):

                # find location of borehole
                z_surface = surface_levels[output_steps[j]]
                xloc = x_loc_fault[j] + xloc_raw

                borehole_xlocs[borehole_number, j] = xloc
                borehole_zlocs[borehole_number, j] = z_loc_fault[j]

                output_number = model_run * n_ts + j

                if mp.report_borehole_xcoords is True:
                    col_name = 'modeled_xcoord_borehole_%s_run_%i_timestep_%i' \
                               % (borehole, model_run, output_steps[j])
                    dft[col_name] = xloc

                z_obs = surface_levels[output_steps[j]] - borehole_depth

                x_buffer = 200.0
                pre_select = np.abs(xyz_array[:, 0] - xloc) < x_buffer

                zi, Ti = beo_core.interpolate_data(xyz_array[pre_select], T_array[j][pre_select], xloc, y_steps=300)

                T_mod = np.interp(z_obs, zi, Ti)

                borehole_temp_modeled[j, :] = T_mod

                col_name = 'modeled_T_%s_run_%i_timestep_%i' \
                           % (borehole, model_run, output_steps[j])
                dft.loc[ind, col_name] = T_mod

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

            borehole_temps_modeled.append(borehole_temp_modeled)

    for j in range(n_ts):

        output_number = model_run * n_ts + j

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
        df.loc[output_number, 'T_change_max'] = T_change.max()
        df.loc[output_number, 'T_change_min'] = T_change.min()

        n_depths = len(target_depths)

        for i in range(n_depths):
            df.loc[output_number, 'max_temperature_layer_%i' % i] = Tzs_cropped[i][j].max()

            for Tc in mp.T_change_report:

                if Tzs_diff_cropped[i][j].max() >= Tc:
                    ind = Tzs_diff_cropped[i][j] >= Tc
                    xc = xzs[i][ind]
                    Tx_min = xc.min()
                    Tx_max = xc.max()
                    Tx = Tx_max - Tx_min

                else:
                    Tx = 0

                df.loc[output_number, 'area_T_change_exc_%0.1f_layer_%i' % (Tc, i)] = Tx

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

            T_surface_i = Tzs_cropped[ind_low][j]
            x_coords_i = xzs[ind_low]

        T_surface.append(T_surface_i)
        x_surface.append(x_coords_i)

        df.loc[output_number, 'max_surface_temperature'] = T_surface_i.max()

        # calculate partial resetting and full resetting distance in the AHe data
        if Ahe_ages_surface_all is not None:

            n_depths = len(AHe_ages_cropped)
            nt_output = AHe_ages_cropped[0].shape[0]

            x_step = 1.0

            for i in range(n_depths):
                ages_raw = AHe_ages_cropped[i][j] / My
                ages_raw_corr = AHe_ages_corr_cropped[i][j] / My
                x_coords_int = np.arange(xzs[i].min(), xzs[i].max() + x_step, x_step)
                ages = np.interp(x_coords_int, xzs[i], ages_raw)
                dev_age = ages / ages.max()

                min_age = np.min(ages)
                max_age = np.max(ages)
                ind_min_age = np.argmin(ages)

                save_data_all_layers = True

                if save_data_all_layers is True:

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

                fraction = np.abs(diff[ind_low]) / (target_depths[ind_high] - target_depths[ind_low])

                ages_raw = AHe_ages_cropped[ind_low][j] / My
                ages_raw_corr = AHe_ages_corr_cropped[ind_low][j] / My
                x_coords = xzs[ind_low]

            # add surface AHe data to output
            AHe_ages_surface.append(ages_raw * My)
            AHe_ages_surface_corr.append(ages_raw_corr * My)
            AHe_xcoords_surface.append(x_coords)

            x_coords_int = np.arange(x_coords.min(), x_coords.max() + x_step, x_step)

            if len(ages_raw) != len(x_coords):
                print('warning, lenght of AHe ages array and x coordinates ' \
                      'array are not equal. Skipping calculating partial/full ' \
                      'reset zone for surface layer at timestep %i of %i' \
                      % (j, n_ts))
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
        AHe_ages_samples_surface_corr = []

        if Ahe_ages_surface_all is not None and mp.calculate_he_ages and mp.model_AHe_samples is True:

            for i in range(N_output_steps):
                surface_elev = surface_levels[output_steps[i]]

                if surface_elev in target_depths:
                    surface_ind = np.where(target_depths == surface_elev)[0][0]
                    ages_raw = AHe_ages_samples_cropped[surface_ind][i]
                    ages_raw_corr = AHe_ages_samples_corr_cropped[surface_ind][i]

                else:
                    # interpolate AHe age from nearest surfaces
                    diff = target_depths - surface_elev
                    ind_low = np.where(diff < 0)[0][-1]
                    ind_high = np.where(diff > 0)[0][0]

                    fraction = np.abs(diff[ind_low]) / (target_depths[ind_high]
                                                        - target_depths[ind_low])

                    ages_raw = AHe_ages_samples_cropped[ind_low][i]
                    ages_raw_corr = AHe_ages_samples_corr_cropped[ind_low][i]

                # add surface AHe data to output
                AHe_ages_samples_surface.append(ages_raw)
                AHe_ages_samples_surface_corr.append(ages_raw_corr)

    # analyze model-data fit of AHe surface samples
    if Ahe_ages_surface_all is not None and mp.model_AHe_samples is True:

        print('analyzing fit of modeled and measured AHe ages')

        n_grains = len(AHe_ages_samples_surface[0])

        me_ahe = np.zeros(N_output_steps)
        mae_ahe = np.zeros(N_output_steps)
        mswd_ahe = np.zeros(N_output_steps)

        me_ahe_corr = np.zeros(N_output_steps)
        mae_ahe_corr = np.zeros(N_output_steps)
        mswd_ahe_corr = np.zeros(N_output_steps)

        AHe_data_file = dfhs

        for timestep in range(N_output_steps):

            output_number2 = model_run * n_ts + timestep

            if mp.profile_number in dfhs['profile'].values:
                profiles = [mp.profile_number]
            else:
                profiles = np.unique(dfhs['profile'])

            print('available profiles in AHe data file: ',  np.unique(dfhs['profile']))
            print('selected profiles: ', profiles)

            for profile in profiles:
                profile_loc = dfhs['profile'] == profile
                dfhs2 = dfhs.loc[profile_loc]

                if mp.profile_number in dfhs['profile'].values:
                    diff = dfhs2['AHe_age_uncorr'].values - AHe_ages_samples_surface[timestep] / My
                    diff_corr = dfhs2['AHe_age_corr'].values - AHe_ages_samples_surface_corr[timestep] / My
                else:
                    diff = dfhs2['AHe_age_uncorr'].values - AHe_ages_samples_surface[timestep][profile_loc] / My
                    diff_corr = dfhs2['AHe_age_corr'].values - AHe_ages_samples_surface_corr[timestep][profile_loc] / My

                me_ahe[timestep] = np.mean(diff)
                mae_ahe[timestep] = np.mean(np.abs(diff))
                mswd_ahe[timestep] = np.sum((diff / (0.5 * dfhs2['AHe_age_uncorr_2se'])) ** 2) / (n_grains - 1)

                me_ahe_corr[timestep] = np.mean(diff_corr)
                mae_ahe_corr[timestep] = np.mean(np.abs(diff_corr))
                mswd_ahe_corr[timestep] = \
                    np.sum((diff_corr / (0.5 * dfhs2['AHe_age_uncorr_2se'])) ** 2) / (n_grains - 1)

                df.loc[output_number2, 'mean_error_AHe_samples_profile%s' % (str(profile))] = me_ahe[timestep]
                df.loc[output_number2, 'mean_abs_error_AHe_samples_profile%s' % (str(profile))] = mae_ahe[timestep]
                df.loc[output_number2, 'mswd_AHe_samples_profile%s' % (str(profile))] = mswd_ahe[timestep]

                df.loc[output_number2, 'mean_error_AHe_samples_corrected_profile%s' % (str(profile))] = \
                    me_ahe_corr[timestep]
                df.loc[output_number2, 'mean_abs_error_AHe_samples_corrected_profile%s' % (str(profile))] = \
                    mae_ahe_corr[timestep]
                df.loc[output_number2, 'mswd_AHe_samples_corrected_profile%s' % (str(profile))] = \
                    mswd_ahe_corr[timestep]

    # option to save corrected ages for figure output

    #
    if exceed_boiling_temp_array is not None:
        exceed_boiling_temp_array_cropped = exceed_boiling_temp_array[output_steps]
        boiling_temp_array_cropped = boiling_temp_array[output_steps]
    else:
        exceed_boiling_temp_array_cropped = None
        boiling_temp_array_cropped = None

    if AHe_ages_borehole is not None:
        AHe_ages_borehole_cropped = AHe_ages_borehole[output_steps]
        AHe_ages_borehole_corrected_cropped = AHe_ages_borehole_corrected[output_steps]
    else:
        AHe_ages_borehole_cropped = None
        AHe_ages_borehole_corrected_cropped = None

    # gather and save model output
    output_selected = \
        [runtimes, runtimes[output_steps], xyz_array,
         surface_levels[output_steps], x_loc_fault, z_loc_fault,
         T_init_array, T_array, boiling_temp_array_cropped,
         xyz_array_exc, exceed_boiling_temp_array_cropped, xyz_element_array,
         qh_array[output_steps], qv_array[output_steps], fault_fluxes,
         durations, xzs, Tzs_cropped,
         x_surface, T_surface, AHe_ages_cropped,
         AHe_ages_corr_cropped, xs_Ahe_surface_all, target_depths,
         AHe_ages_surface, AHe_ages_surface_corr, AHe_xcoords_surface,
         AHe_ages_samples_surface, AHe_ages_samples_surface_corr, AHe_data_file,
         borehole_xlocs, borehole_zlocs, borehole_depths,
         borehole_temp_measured, borehole_temps_modeled,
         z_borehole, AHe_ages_borehole_cropped, AHe_ages_borehole_corrected_cropped]

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month,
                              today.year)

    fn = 'results_model_run_%i_%s_%s_%s.pck' \
         % (model_run, str(param_set), mp.output_fn_adj, today_str)
    fn_path = os.path.join(output_folder, fn)

    print('saving model results as %s' % fn_path)
    try:
        fout = open(fn_path, 'wb')
        pickle.dump(output_selected, fout)
        fout.close()
    except IOError as err:
        print(f"{type(err)}:{err}")
        print(f"check file path for file {fn_path}")
    except TypeError as err:
        print(f"{type(err)}:{err}")
        print(f"could not save model results pickle file {fn_path}")    
    except BaseException as err:
        print(f"Unexpected {err}, {type(err)}")
        raise
    
    print('-' * 30)
    print('saving summary of parameters and model results as %s' % fn_path_csv)
    df.to_csv(fn_path_csv, index_label='row', encoding='utf-8')

    try:
        # save temperatures histories at depth slices as a csv file
        for i in range(n_depths):
            nn = len(Tzs_cropped[i][0])

            tcols = list((runtimes[output_steps] / year).astype(int))
            cols = tcols
            xi = xzs[i]

            dfts = pd.DataFrame(index=list(xi), columns=cols)

            for j in range(n_ts):
                Ti = Tzs_cropped[i][j]
                dfts[tcols[j]] = Ti

            fnout = 'temp_depth_slice_%i_%0.0f_%i_runs_%s_%s.csv' % (i, target_depths[i],
                                                                  len(param_list),
                                                                   mp.output_fn_adj,
                                                                    today_str)
            print('saving modeled temperatures for surface slice to %s' % fnout)
            dfts.to_csv(os.path.join(output_folder, fnout), index_label='x')

    except BaseException as err:
        print('warning, failed to save temperature data')
        print(err)

    # save borehole temperature history
    if mp.analyse_borehole_temp is True:
        fn_new = os.path.split(mp.temperature_file)[-1].split('.')[:-1]
        fn_new = ''.join(fn_new)
        fn_new += '_modeled_%i_runs_%s_%s.csv' % (len(param_list),
                                                  mp.output_fn_adj,
                                                  today_str)
        print('saving modeled temperatures for boreholes to %s' % fn_new)
        dft.to_csv(os.path.join(output_folder, fn_new))

    # save AHe ages to file
    if Ahe_ages_surface_all is not None and mp.save_AHe_ages is True and AHe_ages_surface_all != []:

        nxs = np.max(np.array([AHe_ii.shape[0]
                               for AHe_i in AHe_ages_surface_all
                               for AHe_ii in AHe_i]))
        nts = len(AHe_ages_surface)

        cols = []
        for j in range(n_model_runs):
            for i in range(nts):
                cols += ['x_run_%i_ts%i' % (j, i),
                         'AHe_age_uncorr_run_%i_ts%i' % (j, i),
                         'AHe_age_corrected_run_%i_ts%i' % (j, i)]

        dfh = pd.DataFrame(columns=cols, index=np.arange(nxs))

        for j in range(n_model_runs):
            for i in range(nts):
                a = len(dfh)
                b = len(AHe_xcoords_surface_all[j][i])

                n_ahe_data = np.min((a, b))
                try:
                    dfh.loc[:n_ahe_data-1, 'x_run_%i_ts%i' % (j, i)] = \
                        AHe_xcoords_surface_all[j][i][:n_ahe_data]
                    dfh.loc[:n_ahe_data-1, 'AHe_age_uncorr_run_%i_ts%i' % (j, i)] = \
                        AHe_ages_surface_all[j][i][:n_ahe_data] / My
                    dfh.loc[:n_ahe_data-1, 'AHe_age_corrected_run_%i_ts%i' % (j, i)] = \
                        AHe_ages_surface_corr_all[j][i][:n_ahe_data] / My
                except Exception as msg:
                    print('error, something went wrong with saving AHe data ' \
                          'to a .csv file for model run %i and timestep %i' \
                          % (j, i))
                    print(msg)
                    print('continuing with next timestep')

        fnh = 'AHe_surface_modeled_%i_runs_%s_%s.csv' % (n_model_runs,
                                                         mp.output_fn_adj,
                                                         today_str)

        print('saving modeled AHe ages at the surface to %s' % fnh)
        dfh.to_csv(os.path.join(output_folder, fnh), index_label='row')

    if Ahe_ages_surface_all is not None and mp.model_AHe_samples is True:

        # save AHe ages
        for timestep in range(len(output_steps)):

            # and store in datafile
            col_name = 'modeled_AHe_age_uncorr_run_%i_timestep_%i' % (model_run, output_steps[timestep])
            col_name_corr = 'modeled_AHe_age_corrected_run_%i_timestep_%i' % (model_run, output_steps[timestep])

            if mp.profile_number in dfhs['profile']:
                profiles = [mp.profile_number]
            else:
                profiles = np.unique(dfhs['profile'])

            for profile in profiles:
                profile_loc = dfhs['profile'] == profile

                if True in profile_loc:
                    try:
                        dfhs.loc[profile_loc, col_name] = AHe_ages_samples_surface[timestep] / My
                        dfhs.loc[profile_loc, col_name_corr] = AHe_ages_samples_surface_corr[timestep] / My
                    except Exception as msg:
                        print('error, something went wrong with saving AHe ' \
                              'sample data to a .csv file for model run %i ' \
                              'and timestep %i' % (model_run, timestep))
                        print(msg)
                        print('continuing with next timestep')

        output_fn1 = os.path.split(mp.AHe_data_file)[-1]
        output_fn2 = output_fn1[:-4] + '_modeled_%i_runs_%s_%s.csv' \
                                       % (n_model_runs,
                                          mp.output_fn_adj,
                                          today_str)
        output_fn = os.path.join(output_folder, output_fn2)
        print('saving modeled AHe ages samples to %s' % output_fn)
        dfhs.to_csv(output_fn, index=False)

print('-' * 30)
print('done with all model runs')
