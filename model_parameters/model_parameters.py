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

    # steady state or transient model
    # note that initial condition of transient model =
    # steady-state solution without any advection
    steady_state = False

    # keep the temperature below the max T in the vapour pressure curve
    vapour_correction = True

    # model dimensions
    width = 6000.0
    total_depth = 6000.0
    air_height = 100.0

    # depth to fine discretization at surface:
    z_fine = -100

    # grid size
    cellsize = 500.0
    cellsize_air = 10.0
    cellsize_fault = 5.0
    cellsize_fine = 20.0
    cellsize_base = 1000.0

    # exhumation parameters
    # exhumation rate in m/yr
    # assuming the AHe was not reset, the max exhumation is ~1500 m in 15 My = 1e-4 m/yr
    # look up regional AFT, AHe and cosmogenic nuclide work for realistic range of exhumation rates
    exhumation_rate = 1e-4

    # number of grid layers between initial and final surface level
    # the more layers, the more smooth and accurate the exhumation history,
    # but this also slows the model down somehwat
    exhumation_steps = 3

    # number of timesteps after which the surface level is recalculated
    # ideally this should be 1 (ie recalculate at each timestep)
    # higher number means faster model
    exhumation_interval = 10

    # temperature bnd conditions
    air_temperature = 10.0
    #bottom_temperature = total_depth * 0.03 + air_temperature
    #bottom_temperature = total_depth * 0.03
    # new version: calculate bottom T using a fixed geothermal gradient./r
    thermal_gradient = 0.05

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
    layer_bottom = [[-20000, -20000],
                    [-500.0, - 450.0],
                    [-200.0, -150.0]]

    # porosity for each layer
    porosities = [0.1, 0.15, 0.25]

    # thermal parameters
    # note that only thermal conductivity is varied between layers, the rest
    # is constant
    K_solids = [2.5, 2.5, 2.5]

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
    N_outputs = [10]

    # size of timestep
    dt = 1000 * year

    # duration of each timestep_slice
    durations = [1e4 * year]

    # target depth slices for calculating temperature and U-Th/He
    # in case of exhumation, this values is overridden and
    # set equal to each exhumation step layer
    # in this way one can track the AHe response to each layer that
    # comes to the surface in a longer time period
    target_zs = [10.0, 5.0, 0.0]

    # U-Th/He params
    calculate_he_ages = True

    # method to calculate helium diffusivity, use Wolf1996, Farley2000 or RDAAM
    AHe_method = 'RDAAM'

    # temperature after crystallization and before hydrothermal heating
    T0 = 10.0
    T_surface = 10.0

    # crystallization age
    t0 = 15.2 * My

    # apatite params
    radius = 100.0 * 1e-6
    U238 = 8.98e-6
    Th232 = 161.3e-6

    # alpha ejection parameters:
    alpha_ejection = True
    # alpha ejection stopping distance, see Ketcham (2011) for estimates
    stopping_distance = 21e-6

    ## fault data for multiple faults:

    # x location of fault:
    fault_xs = [4000]

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

    aquifer_bottoms = [None]
    aquifer_tops = [None]
    aquifer_fluxes = [None]

    # relative limit to consider a sample partial reset or not, ie if 0.95
    # a sample will be considered partially reset if the modeled uncorrected
    # AHe age is less than 0.95 x the maximum age in the system.
    partial_reset_limit = 0.75

    # absolute limit below which samples are considered reset (ie. AHe age ~0 My)
    reset_limit = 0.01
