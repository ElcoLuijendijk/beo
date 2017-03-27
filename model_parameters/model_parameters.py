"""
file containing all parameters for Beo.py

"""

__author__ = 'Elco Luijendijk'


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
    output_fn_adj = 'profile1'

    # steady state or transient model
    # note that regardless of this setting, the initial condition of transient model is
    # the steady-state solution without any advection
    steady_state = False

    # keep the temperature below the max T in the vapour pressure curve
    vapour_correction = True

    # model dimensions
    width = 8000.0
    total_depth = 8000.0
    air_height = 40.0

    # depth to fine discretization at surface:
    # deprecated, doesnt work, maybe need to add this in again to increase numerical stability
    z_fine = -100

    # grid size
    cellsize = 200.0
    cellsize_air = 5.0
    cellsize_fault = 5.0
    cellsize_base = 1000.0

    # new: buffer zone around fault with the same cell size as the fault
    # this is to reduce model instability
    use_mesh_with_buffer = False
    fault_buffer_zone = 25.0
    cellsize_fine = 10.0

    # exhumation parameters
    # exhumation rate in m/yr
    # assuming the AHe was not reset, the max exhumation is ~1500 m in 15 My = 1e-4 m/yr
    # look up regional AFT, AHe and cosmogenic nuclide work for realistic range of exhumation rates
    exhumation_rate = 1e-4

    # number of grid layers between initial and final surface level
    # the more layers, the more smooth and accurate the exhumation history,
    # but this also slows the model down somewhat
    exhumation_steps = 20

    # number of timesteps after which the surface level is recalculated
    # ideally this should be 1 (ie recalculate at each timestep)
    # higher number means faster model
    # note: this parameter is no longer used. Now surface level is recalculated
    # at each timestep....
    #exhumation_interval = 10

    # temperature bnd conditions
    air_temperature = 10.0
    #bottom_temperature = total_depth * 0.03 + air_temperature
    #bottom_temperature = total_depth * 0.03
    # new version: calculate bottom T using a fixed geothermal gradient./r
    thermal_gradient = 0.04

    # bottom flux bnd condition, set to None if T bnd is used
    basal_heat_flux = None
    #basal_heat_flux = 65e-3

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
    #85-18
    layer_bottom = [[-20000.0, -20250.0],
                    [-810.0, -1060.0],
                    [-530.0, -780.0],
                    [-440.0, -690.0],
                    [-210.0, -460.0]]

    # porosity for each layer
    porosities = [0.08, 0.05, 0.25, 0.17, 0.05]

    # thermal parameters
    # note that only thermal conductivity is varied between layers, the rest
    # is constant
    K_solids = [4.44, 2.26, 1.58, 1.6, 2.0]

    # wikipedia: thermal properties air
    # heat transfer coefficient = 10- 100 W / (m2 K))
    # heat capacity = 1000 J kg-1 K-1
    # density = 1.29 kg m-3
    K_air = 50.0
    K_water = 0.58

    rho_air = 1.29
    rho_f = 1000.0
    rho_s = 2650.0

    c_air = 1000.0
    c_f = 4000.
    c_s = 900.

    # timesteps
    # number of output steps
    # this is not used when exhumation > 0, in this case output is generated
    # once each new surface level is reached
    # the number of surfaces is controlled by the exhumation_steps parameter
    N_outputs = [50]

    # size of timestep
    dt = 500.0 * year

    # duration of each timestep_slice
    durations = [5e5 * year]

    # target depth slices for calculating temperature and U-Th/He
    # in case of exhumation, this values is overridden and
    # set equal to each exhumation step layer
    # in this way one can track the AHe response to each layer that
    # comes to the surface in a longer time period
    target_zs = [10.0, 5.0, 0.0]

    # U-Th/He params
    calculate_he_ages = True

    # model-data comparison AHe samples
    model_AHe_samples = True
    AHe_data_file = 'model_parameters/AHe_data.csv'
    profile_number = 1

    #save the AHe ages qat the surface to a separate file
    # do not turn this on yet, there's a bug in this part of the code....
    save_AHe_ages = True

    # method to calculate helium diffusivity, use Wolf1996, Farley2000 or RDAAM
    AHe_method = 'RDAAM'

    # use temperature at each x timestep for the calculation of AHe ages
    # this should be an integer. Ideally this should be 1, but higher numbers
    # significantly speed up the model
    # !! new parameter (3 march 2017)
    AHe_timestep_reduction = 4

    # temperature after crystallization and before hydrothermal heating
    T0 = 10.0
    T_surface = 10.0

    # crystallization age
    t0 = 15.6 * My

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
    # alpha ejection stopping distance, see Ketcham (2011) for estimates
    stopping_distance = 21e-6

    ## fault data for multiple faults:

    # x location of fault:
    fault_xs = [5000]

    # fault width
    fault_widths = [20.0]

    # angle of the fault zone (degrees), dip of normal faults ~60-70 degrees
    fault_angles = [-65.0]

    # elevation of bottom of fault
    fault_bottoms = [-5000.0]

    # fluid advection rates in faults:
    # nested list,
    # [[fault1_t1, fault2_t1], [fault1_t2, fault2_t2], etc...]
    # note units are m2/sec, ie the integrated flux over the entire width of the
    # fault zone
    fault_fluxes = [[-400.0 / year]]

    # aquifers, use aquifer_top = [None] to not use this:
    #aquifer_tops = [-0.0]
    aquifer_tops = [-50.0]
    aquifer_bottoms = [-100.0]
    aquifer_fluxes = [[-250.0 / year]]

    # left side of aquifer. right hand bnd is assumed to be the fault zone
    aquifer_left_bnds = [2000.0]

    # relative limit to consider a sample partial reset or not, ie if 0.95
    # a sample will be considered partially reset if the modeled uncorrected
    # AHe age is less than 0.95 x the maximum age in the system.
    partial_reset_limit = 0.75

    # absolute limit below which samples are considered reset (ie. AHe age ~0 My)
    reset_limit = 0.1

    # option to calculate temperature data for one or several boreholes
    # note that there seems to be a bug in the output timesteps for the temperature calculation
    # avoid using this for now...
    analyse_borehole_temp = True

    # file that contains temperature data
    temperature_file = 'model_parameters/temperature_data.csv'
    borehole_names = ['85-18']

    # locations of boreholes for temperature data
    borehole_xs = [4800.0]
