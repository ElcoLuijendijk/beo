"""
file containing all parameters for Beo.py

this file contains the base parameters for a standard model run

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
    output_fn_adj = 'beowawe_AHe_cont'

    solver = 'GMRES'
    # steady state or transient model
    # note that regardless of this setting, the initial condition of transient model is
    # the steady-state solution without any advection
    steady_state = False

    # iterations for steady-state runs
    # these are needed if either vapour correction or variable K air are used
    n_iterations_steady_state = 10

    # keep the temperature below the max T in the vapour pressure curve
    vapour_correction = True

    # model dimensions
    width = 2000.0
    total_depth = 8000.0
    air_height = 2.0

    # depth to fine discretization near surface:
    z_fine = -150

    # default cellsize
    cellsize = 250.0

    # cellsize in the air layer:
    cellsize_air = 5.0

    # cellsize at surface layers:
    cellsize_surface = 50.0

    # fine cellsize near surface (up to depth = z_fine)
    cellsize_fine = 100.0

    # in fault zone:
    cellsize_fault = 2.5

    # cellsize at the lower left and right corners:
    cellsize_base = 500.0

    # new: buffer zone around fault with the same cell size as the fault
    # this is to reduce model instability
    use_mesh_with_buffer = False
    fault_buffer_zone = 20.0

    ## exhumation parameters
    # add exhumation or not
    add_exhumation = True

    # exhumation rate in m/yr
    exhumation_rate = 1e-4

    # number of grid layers between initial and final surface level
    # the more layers, the more smooth and accurate the exhumation history,
    # but this also slows the model down somewhat
    exhumation_steps = 50

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
    #exhumation_interval = 10

    # temperature bnd conditions
    air_temperature = 10.0

    # calculate bottom T using a fixed geothermal gradient.
    # use None is you want to use a specified basal heat flux instead
    #thermal_gradient = None
    thermal_gradient = 0.04

    # bottom flux bnd condition, set to None if T bnd is used
    #basal_heat_flux = 65e-3
    #basal_heat_flux = 100e-3

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

    # layers: valmy, basaltic andesite, tuff, basalt, dacite
    # note that the alluvial sediments are not included, since theyre only found in the valley

    #layer_bottom = [[-20000.0, -20250.0],
    #                [-810.0, -1060.0],
    #                [-530.0, -780.0],
    #                [-440.0, -690.0],
    #                [-210.0, -460.0]]

    # well Collins 76-17: Dacite at 370 m, basalt: 550 m, tuff: 720 m, basaltic andesite: 950 m, conglormerate: 1040 m, Valmy fm....

    layer_bottom = [[-20230, -20000],
                    [-1180., -950],
                    [-950., -720],
                    [-780., -550],
                    [-600., -370]]

    # porosity for each layer
    porosities = [0.08, 0.05, 0.20, 0.05, 0.05]

    # thermal parameters
    # note that only thermal conductivity is varied between layers,
    # the other parameters are constant
    K_solids = [4.44, 2.26, 1.50, 1.6, 2.0]

    # thermal properties air, solid matrix and porewater:
    # note K_air is not used if variable_K_air is set to True
    K_air = 100.0
    K_water = 0.58

    rho_air = 1.29
    rho_f = 1000.0
    rho_s = 2650.0

    c_air = 1000.0
    c_f = 4000.
    c_s = 900.

    ## variable heat transfer coefficient for air
    variable_K_air = True

    ## parameters to estimate heat transfer coefficient of air:
    # aerodynamic resistance, see Liu et al (2007), Hydrology and Earth System Sciences 11 (2)
    ra = 80
    # measurement height for aerodynamic resistance
    dz = 1.8

    # number of output steps for each timeslice
    N_outputs = [30]

    # size of a single timestep
    dt = 250.0 * year

    # size of timestep to store model results. make this higher than dt if you want to conserve memory,
    # otherwise make this the same as dt
    dt_stored = 1000.0 * year

    # duration of each timestep_slice
    #durations = [2e3 * year, 8e3 * year]
    durations = [3e5 * year]

    # repeat timesclices x times, use this to model repeated episodic heating events
    # set this to zero or None to not use this
    repeat_timeslices = None

    # target depth slices for calculating temperature and U-Th/He
    # in case of exhumation, this values is overridden and
    # set equal to each exhumation step layer
    # in this way one can track the AHe response to each layer that
    # comes to the surface in a longer time period
    # note, values should go from low/bottom to high/top
    # warning: the lowest value in target_zs should always be higher than z_fine
    target_zs = [0.0]

    # U-Th/He params
    calculate_he_ages = True

    # model-data comparison AHe samples
    model_AHe_samples = True
    AHe_data_file = 'model_parameters/AHe_data.csv'
    profile_number = 0

    #save the AHe ages at the surface to a separate file
    save_AHe_ages = False

    # method to calculate helium diffusivity, use Wolf1996, Farley2000 or RDAAM
    AHe_method = 'RDAAM'

    # use temperature at each x timestep for the calculation of AHe ages
    # this should be an integer. Ideally this should be 1, but higher numbers
    # significantly speed up the model code
    # !! new parameter (3 march 2017)
    AHe_timestep_reduction = 2

    # crystallization age
    t0 = 15.6 * My

    # temperature of apatites after crystallization and before hydrothermal heating
    T0 = 10.0
    T_surface = 10.0

    # apatite params
    radius = 100.0 * 1e-6
    U238 = 8.98e-6
    Th232 = 161.3e-6

    #distance_samples = [1500.0, 2000.0]
    #radius_samples = [100.0 * 1e-6, 150e-6]
    #U238_samples = [8.98e-6]
    #Th232_samples = [161.3e-6]

    #AHe_age_samples = [10e6, 15e6]

    # alpha ejection parameters:
    alpha_ejection = True

    # alpha ejection stopping distance (um), see Ketcham (2011) for estimates
    stopping_distance = 21e-6

    ## fault data for multiple faults:

    # x location of fault (m):
    fault_xs = [5000.0]

    # fault width (m)
    fault_widths = [10.0]

    # angle of the fault zone (degrees), dip of normal faults ~60-70 degrees
    fault_angles = [-65.0]

    # elevation of bottom of fault
    fault_bottoms = [-5000.0]

    # different segments of the fault, list of the top bnd of each segments starting from the bottom
    # nested list: [[segment_top1_fault1, segment_top2_fault1], [segment_top1_fault2, segment_top2_fault2], etc...]
    fault_segments = [[5000.0]]

    # fluid advection rates in faults:
    # nested list,
    # [[[fault1_segment1_t1, fault1_segment2_t1], [fault2_segment1_t1, fault2_segment2_t1], etc...]
    # note units are m2/sec, ie the integrated flux over the entire width of the
    # fault zone
    # length of fault = 1.5 km, total discharge is 18 kg/s = 0.018 m3/s = 0.018 / 1500 = 378.4 m2/year
    fault_fluxes = [[[-380 / year]], [[0.0]]]

    # aquifers, used for modeling horizontal advective flow
    # use aquifer_top = [None] to not use this:
    # note for multiple aquifers start at the lowest aquifer
    # aquifer_tops = [-350.0, -50.0]
    aquifer_tops = [40.0]
    #aquifer_bottoms = [-550.0, -200.0]
    aquifer_bottoms = [-50.0]
    #aquifer_fluxes = [[100.0 / year, -350.0 / year]]
    # sideways flux is 2/3 of total (Olmsted & Rush) = 2/3. * 380.0 = 250 m2/year
    # aquifer flux: nested list
    # [[aquifer1_timestep1, aquifer2_timestep1], [aquifer1_timestep2, aquifer2_timestep2]]
    aquifer_fluxes = [[-250.0 / year], [0.0]]
    ## left side of aquifer. right hand bnd is assumed to be the fault zone
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
    borehole_names = ['85-18']
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

    # option to add a first base run with unchanged parameters to the lsit of model
    # runs
    initial_base_run = True

    ###################################################################
    # parameters that will be changed in the sensitivity analysis runs:
    ###################################################################

    #profile_number_s = [1, 2]
    #durations_s = [[1e3 * year, 9e3 * year], [2e3 * year, 8e3 * year], [3e3 * year, 7e3 * year]]
