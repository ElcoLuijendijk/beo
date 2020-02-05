"""
file containing all parameters for Beo.py

this file contains the parameters to run the example for the Baden springs as discussed in the
Geoscientific Model Development paper

"""

__author__ = 'Elco Luijendijk'

import numpy as np


class ModelParams:

    # constants
    degree_symbol = unichr(176)
    day = 24.0 * 60.0 * 60.0
    year = 365.25 * day
    My = year * 1e6
    Kelvin = 273.15

    # directory to store model output
    output_folder = 'model_output'

    #
    output_fn_adj = 'baden_7km_slow'

    #
    save_VTK_file = True

    # solver, see escript documentation for details
    # available choices: 'PCG', 'DIRECT', 'GMRES', 'ROWSUM_LUMPING'
    solver = 'GMRES'

    # steady state or transient model
    # note that regardless of this setting, the initial condition of transient model is
    # the steady-state solution without any advection
    steady_state = False

    # iterations for steady-state runs
    # these are needed if either vapour correction or variable K air are used
    n_iterations_steady_state = 10

    # keep the temperature below the max T in the vapour pressure curve
    vapour_correction = False

    # model dimensions
    width = 3000.0
    total_depth = 10000.0
    air_height = 2.0

    # depth to fine discretization near surface:
    z_fine = -100

    # default cellsize
    cellsize = 500.0

    # cellsize in the air layer:
    cellsize_air = 100.0

    # cellsize at fault surface:
    cellsize_fault_surface = 0.5

    # cellsize at land surface:
    cellsize_surface = 100.0

    # fine cellsize near surface (up to depth = z_fine)
    cellsize_fine = 100.0

    # in fault zone:
    cellsize_fault = 2.5

    # cellsize at the lower left and right corners:
    cellsize_base = 1000.0

    # new: buffer zone around fault with the same cell size as the fault
    # this is to reduce model instability
    use_mesh_with_buffer = False
    fault_buffer_zone = 20.0

    # exhumation parameters
    # add exhumation or not
    add_exhumation = False

    # exhumation rate in m/yr
    exhumation_rate = 0.0

    # number of grid layers between initial and final surface level
    # the more layers, the more smooth and accurate the exhumation history,
    # but this also slows the model down somewhat
    exhumation_steps = 25

    # minimum layer thickness, if the exhumation steps result in surfaces that
    # are less than the min thickness apart, the number of steps is reduced
    # default value is 1.0 m, reduce this value if gmsh returns an error while
    # creating the mesh
    min_layer_thickness = 1.0

    # number of timesteps after which the surface level is recalculated
    # ideally this should be 1 (ie recalculate at each timestep)
    # higher number means faster model
    # note: this parameter is no longer used. Now surface level is recalculated
    # at each timestep....
    # exhumation_interval = 10

    # temperature bnd conditions
    air_temperature = 10.0

    # calculate bottom T using a fixed geothermal gradient.
    # use None is you want to use a specified basal heat flux instead
    thermal_gradient = None

    # bottom flux bnd condition, set to None if T bnd is used
    basal_heat_flux = 70e-3

    # elevation of layers either side of the fault
    # structured like this:
    # [[depth layer 1 left, depth layer 1 right],
    #  [depth layer 2 left, depth layer 2 right],
    #  [depth layer 3 left, depth layer 3 right],
    #  etc...
    # ]
    # layers are counted from bottom to top
    # leave depth of first layer at arbitrarily high value to make sure the entire
    # model domain is covered
    # note that currently only 1 fault is taken into account...

    layer_bottom = [[-20000, -20000]]

    # porosity for each layer
    porosities = [0.15]

    # thermal parameters
    # note that only thermal conductivity is varied between layers,
    # the other parameters are constant
    K_solids = [2.5]

    # thermal properties air, solid matrix and porewater:
    # note K_air is not used if variable_K_air is set to True
    K_air = 50.0
    K_water = 0.58

    rho_air = 1.29
    rho_f = 1000.0
    rho_s = 2650.0

    c_air = 1000.0
    c_f = 4000.
    c_s = 900.

    # variable heat transfer coefficient for air
    variable_K_air = True

    # parameters to estimate heat transfer coefficient of air:
    ra = 80

    # measurement height for aerodynamic resistance
    dz = 2.0

    # number of output steps for each timeslice
    N_outputs = [15]

    # size of a single timestep
    dt = 0.005 * year

    # size of timestep to store model results. make this higher than dt if you want to conserve memory,
    # otherwise make this the same as dt
    dt_stored = 10.0 * year

    # duration of each timestep_slice
    durations = [1e3 * year]

    # repeat timesclices x times, use this to model repeated episodic heating events
    # set this to zero or None to not use this
    repeat_timeslices = 0

    # target depth slices for calculating temperature and U-Th/He
    # in case of exhumation, this values is overridden and
    # set equal to each exhumation step layer
    # in this way one can track the AHe response to each layer that
    # comes to the surface in a longer time period
    # note, values should go from low/bottom to high/top
    # warning: the lowest value in target_zs should always be higher than z_fine
    target_zs = [0.0]

    # U-Th/He params
    calculate_he_ages = False

    # model-data comparison AHe samples
    model_AHe_samples = True
    AHe_data_file = 'model_parameters/AHe_data.csv'
    profile_number = 0

    # save the AHe ages at the surface to a separate file
    save_AHe_ages = False

    # method to calculate helium diffusivity, use Wolf1996, Farley2000 or RDAAM
    AHe_method = 'RDAAM'

    # use temperature at each x timestep for the calculation of AHe ages
    # this should be an integer. Ideally this should be 1, but higher numbers
    # significantly speed up the model code
    # !! new parameter (3 march 2017)
    AHe_timestep_reduction = 1

    # crystallization age
    t0 = 15.6 * My

    # temperature of apatites after crystallization and before hydrothermal heating
    T0 = 10.0
    T_surface = 10.0

    # apatite params
    radius = 100.0 * 1e-6
    U238 = 8.98e-6
    Th232 = 161.3e-6

    # alpha ejection parameters:
    alpha_ejection = True

    # alpha ejection stopping distance (um), see Ketcham (2011) for estimates
    stopping_distance = 21e-6

    # x location of fault (m):
    fault_xs = [0.0]

    # fault width (m)
    fault_widths = [10.0]

    # angle of the fault zone (degrees), dip of normal faults ~60-70 degrees
    fault_angles = [-65.0]

    # elevation of bottom of fault
    fault_bottoms = [-7000.0]

    # different segments of the fault, list of the top bnd of each segments starting from the bottom
    # nested list: [[segment_top1_fault1, segment_top2_fault1], [segment_top1_fault2, segment_top2_fault2], etc...]
    fault_segments = [[5000.0]]

    # fluid advection rates in faults:
    # nested list,
    # [[[fault1_segment1_t1, fault1_segment2_t1], [fault2_segment1_t1, fault2_segment2_t1], etc...]
    # note units are m2/sec, ie the integrated flux over the entire width of the fault zone
    fault_fluxes = [[[-2e-5]]]

    # aquifers, used for modeling horizontal advective flow
    # use aquifer_top = [None] to not use this:
    # note for multiple aquifers start at the lowest aquifer
    aquifer_tops = [None]
    aquifer_bottoms = [-50.0]
    # aquifer flux: nested list
    # [[aquifer1_timestep1, aquifer2_timestep1], [aquifer1_timestep2, aquifer2_timestep2]]
    aquifer_fluxes = [[-250.0 / year], [0.0]]
    # left side of aquifer. right hand bnd is assumed to be the fault zone
    aquifer_left_bnds = [-1000.0, -1000.0]

    # relative limit to consider a sample partial reset or not, ie if 0.95
    # a sample will be considered partially reset if the modeled uncorrected
    # AHe age is less than 0.95 x the maximum age in the system.
    partial_reset_limit = 0.75

    # absolute age limit below which samples are considered reset (ie. AHe age ~0 My)
    reset_limit = 0.1

    # option to calculate temperature data for one or several boreholes
    # note that there seems to be a bug in the output timesteps for the temperature calculation
    # avoid using this for now...
    analyse_borehole_temp = False

    # file that contains temperature data
    temperature_file = 'model_parameters/temperature_data.csv'
    borehole_names = ['dummy']
    report_borehole_xcoords = False

    # locations of boreholes for temperature data,
    # !! note location is now relative to the location of the first fault
    # ie, -100 m means 100 m to the left of the fault.
    # the model code automatically calculates the correct position to take
    # into account the changing position of the fault surface over time
    # due to exhumation
    borehole_xs = [-250.0]

    # temperature changes to report, report area in depth slices with temperature changes > x degrees
    T_change_report = [10.0]


class ParameterRanges:

    """
    parameter ranges for sensitivity or uncertainty analysis

    beo.py will look for any variable ending with _s below and then look for the
    corresponding variable in model_parameters.py

    each _s variable should be a list of values, beo.py will replace the variable
    in model_parameters.py with each item in the list consecutively
    """

    year = 365.25 * 24 * 60 * 60.0

    # option whether to vary one model parameter at a time
    # (ie for a sensitivtiy analysis)
    # or to run all parameter combinations, using the parameter ranges specified below
    parameter_combinations = False

    # option to add a first base run with unchanged parameters to the list of model
    # runs
    initial_base_run = True

    # parameters that will be changed in the sensitivity analysis runs:
    #fault_bottoms_s = [[-3000.0], [-4000.0], [-5000.0], [-6000.0], [-7000.0], [-8000.0], [-9000.0]]

    #dt_s = [10 * year]

    #cellsize_fault_s = [10.0, 5.0, 2.5]