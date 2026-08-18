"""
file containing all parameters for Beo.py

this file contains an example set of parameters for a radially symmetric model of a submarine hydrothermal vent system, 
with a fault that acts as a conduit for the hot water. The fault simplified as a cylindrical vertical conduit.

"""

__author__ = 'Elco Luijendijk'

import numpy as np


class ModelParams:

    # constants
    degree_symbol = chr(176)
    day = 24.0 * 60.0 * 60.0
    year = 365.25 * day
    My = year * 1e6
    Kelvin = 273.15

    # directory to store model output
    output_folder = 'model_output'

    #
    output_fn_adj = 'radial_submarine'

    #
    save_VTK_file = False

    # numerical backend used to solve the heat flow equation
    # available choices: 'escript' and 'fipy'
    # the escript backend uses the finite element method, the fipy backend
    # uses the finite volume method and does not require escript to be
    # installed. install the fipy backend with pip install fipy gmsh
    backend = 'fipy'

    # solve the heat flow equation assuming radial symmetry around a vertical
    # axis at r = 0, instead of a cartesian cross section. the model domain is
    # then a vertical section from the symmetry axis to a radius of width, and
    # the fault is a cylindrical conduit centered on the axis.
    # only supported by the fipy backend
    radial_symmetry = True

    # discretization of the advection term for the fipy backend
    # available choices: 'powerlaw', 'exponential', 'upwind' and 'central'
    fipy_convection_scheme = 'powerlaw'

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
    # note that the vapour correction assumes atmospheric pressure at the top
    # of the model domain and is therefore not suitable for a submarine model,
    # where the pressure at the seafloor is much higher. boiling is strongly
    # suppressed at these pressures anyway
    vapour_correction = False

    # model dimensions
    # for a radially symmetric model width is the outer radius of the model
    # domain, measured from the symmetry axis at r = 0
    width = 4000.0

    # deep enough to contain the thermal drawdown below the fault
    total_depth = 4000.0

    # thickness of the seawater layer above the seafloor
    air_height = 2.0

    # depth to fine discretization near surface:
    z_fine = -100

    # default cellsize
    cellsize = 300.0

    # cellsize in the air layer:
    cellsize_air = 20.0

    # cellsize at fault surface:
    # the fault conduit has a radius of only 0.15 m, so the cells at the
    # seafloor end of the conduit have to be much smaller than that
    cellsize_fault_surface = 0.05

    # cellsize at land surface:
    cellsize_surface = 20.0

    # fine cellsize near surface (up to depth = z_fine)
    cellsize_fine = 30.0

    # in fault zone:
    # keep this larger than the width of the conduit. the width of the
    # conduit then sets the cell size at depth, which avoids a mesh with
    # millions of cells in the 850 m long conduit
    cellsize_fault = 1.0

    # cellsize at the lower left and right corners:
    cellsize_base = 150.0

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
    # temperature of the seawater above the seafloor
    air_temperature = 8.0

    # calculate bottom T using a fixed geothermal gradient.
    # use None is you want to use a specified basal heat flux instead
    thermal_gradient = None

    # bottom flux bnd condition, set to None if T bnd is used
    # heat flux at the base of the model domain. 
    basal_heat_flux = 60.24e-3

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
    # thermal conductivity of the layer above the seafloor. this is the
    # thermal conductivity of seawater, so heat is only transported by
    # conduction in the seawater and not by convection
    K_air = 0.6
    K_water = 0.58

    rho_air = 1025.0
    rho_f = 1000.0
    rho_s = 2650.0

    c_air = 3990.0
    c_f = 4000.
    c_s = 900.

    # variable heat transfer coefficient for air
    # do not use the heat transfer coefficient for sensible and latent heat
    # flux in air, that equation only applies to a land surface in contact with
    # the atmosphere. for a submarine model the layer above the seafloor is
    # seawater and only conducts heat
    variable_K_air = False

    # parameters to estimate heat transfer coefficient of air:
    ra = 80

    # measurement height for aerodynamic resistance
    dz = 2.0

    # number of output steps for each timeslice
    N_outputs = [15]

    # size of a single timestep. if dt_growth_factor is set below, this is
    # instead the size of the first timestep of each timeslice
    dt = 0.1 * year

    # optionally grow the timestep dt geometrically, by a factor of
    # dt_growth_factor each timestep, until it reaches dt_max, after which
    # dt stays constant at dt_max. this is useful to resolve a fast initial
    # transient, for instance right after a fault zone starts discharging
    # water, with a small timestep, without having to run the entire
    # timeslice at that small timestep. the timestep starts growing again at
    # the start of every timeslice. only supported by the fipy backend, and
    # only in combination with stored_timesteps, since the regular dt_stored
    # output interval assumes a constant dt. leave dt_growth_factor at None
    # to use a constant timestep dt instead
    # dt_growth_factor = 1.1
    # dt_max = 250.0 * year
    dt_growth_factor = 1.2
    dt_max = 250 * year

    # size of timestep to store model results. make this higher than dt if you want to conserve memory,
    # otherwise make this the same as dt
    # every stored timestep is interpolated to the mesh vertices, while only
    # N_outputs timeslices end up in the model output file. storing the results
    # once every duration / N_outputs, so once every 1000 years here, stores 16
    # temperature fields instead of 151 and gives exactly the same output file,
    # which saves about 40 percent of the runtime of this model. note that the
    # AHe ages are calculated from the full stored temperature history, so a
    # larger value of dt_stored changes the modeled AHe ages when
    # calculate_he_ages is switched on
    dt_stored = 1000.0 * year

    # alternative to dt_stored: store results at an explicit list of times
    # instead of a regular interval, useful to resolve a part of the model
    # run in more detail than the rest, for instance the first few timesteps
    # after a fault zone starts discharging water. this is a single list of
    # times in seconds covering the entire model run, not one list per
    # timeslice, and overrides dt_stored and N_outputs when set. only
    # supported by the fipy backend. leave at None to use dt_stored instead.
    # the initial condition at t=0 is always stored regardless of this list,
    # so there is no need to include 0 here
    # stored_timesteps = [10 * year, 50 * year, 100 * year, 500 * year,
    #                     1000 * year, 5000 * year, 15e3 * year]

    stored_timesteps = [0.1 * year, 0.5 * year, 1 * year, 5 * year, 
                        10 * year, 50 * year, 100 * year, 500 * year,
                        1000 * year, 5000 * year, 10e3 * year, 15e3 * year]

    # duration of each timestep_slice
    durations = [15e3 * year]

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

    # calculate thermochron ages for surface layer, including the effects of exhumation
    model_thermochron_surface = False

    # model thermochron ages in a borehole section
    model_thermochron_borehole = True

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

    # cooling rate for calculating pre-hydrothermal activity ages (degr. C / sec)
    # this is only used when modeling borehole temperature samples
    t_cooling = [4.0 * My, 3.0 * My, 1.0 * My]
    cooling_rates = [(1.5 * 27.) / My, 5.0 / My, (1.25 * 27.) / My]

    # apatite params
    radius = 100.0 * 1e-6
    U238 = 8.98e-6
    Th232 = 161.3e-6

    # alpha ejection parameters:
    alpha_ejection = True

    # alpha ejection stopping distance (um), see Ketcham (2011) for estimates
    stopping_distance = 21e-6

    # helium diffusion parameters Farley 2000 model:
    D0 = 50.0 / 1e4
    Ea = 32.9 * 4184.0

    # RDAAM parameters (see Flowers et al. 2009, table 1):
    log_omega_p = -22.0
    log_phi_p = -13.0
    Etrap = 34.0 * 1000.0
    ln_D0_L_div_a2 = 9.733
    E_L = 122.3 * 1000.0

    # x location of fault (m):
    # for a radially symmetric model the fault has to be on the symmetry axis
    fault_xs = [0.0]

    # fault width (m)
    # for a radially symmetric model the fault is a cylindrical conduit with a
    # radius of half the fault width, so a fault width of 50 m gives a conduit
    # with a radius of 25 m
    # a conduit radius of 0.15 m, so a fault width of 0.3 m
    fault_widths = [0.3]

    # angle of the fault zone (degrees), dip of normal faults ~60-70 degrees
    # for a radially symmetric model the fault has to be vertical
    fault_angles = [90.0]

    # elevation of bottom of fault
    # depth of the fault.
    fault_bottoms = [-470.0]

    # different segments of the fault, list of the top bnd of each segments starting from the bottom
    # nested list: [[segment_top1_fault1, segment_top2_fault1], [segment_top1_fault2, segment_top2_fault2], etc...]
    fault_segments = [[5000.0]]

    # fluid advection rates in faults:
    # nested list,
    # [[[fault1_segment1_t1, fault1_segment2_t1], [fault2_segment1_t1, fault2_segment2_t1], etc...]
    # note that for a radially symmetric model the units are m3/sec, ie the
    # total volume of water flowing up through the fault conduit. this is
    # divided by the cross sectional area of the conduit (pi r^2) to get the
    # darcy flux in m/sec. a positive value means upward flow.
    # 1e-3 m3/sec is a discharge of 1 litre per second
    fault_fluxes = [[[35e-3]]]

    # aquifers, used for modeling horizontal advective flow
    # use aquifer_top = [None] to not use this:
    # note for multiple aquifers start at the lowest aquifer
    aquifer_tops = [-370.0]
    aquifer_bottoms = [-470.0]
    # aquifer flux: nested list
    # [[aquifer1_timestep1, aquifer2_timestep1], [aquifer1_timestep2, aquifer2_timestep2]]
    # the aquifer feeds the fault with the same flow of 35 litres per second
    # that the fault discharges at the vent. for a radially symmetric model the
    # units are m3/sec and a negative value means flow towards the symmetry
    # axis, so towards the fault. the flow is divided by the circumference at
    # each radius and by the thickness of the aquifer, so the darcy flux
    # decreases with the distance from the fault while the total flow stays the
    # same at every radius
    aquifer_fluxes = [[-35e-3], [0.0]]
    # left side of aquifer. right hand bnd is assumed to be the fault zone
    # for a radially symmetric model this is the inner edge of the aquifer,
    # None means that the aquifer starts at the fault
    aquifer_left_bnds = [None]
    # right side of aquifer. Use None to use the fault as the right-hand boundary
    # the aquifer extends to half of the radius of the model domain
    aquifer_right_bnds = [2000.0]

    # a radially symmetric model requires horizontal aquifers
    aquifer_angles = [0.0]

    # by default Beo assumes that an aquifer drains the fault and switches off
    # the flow in the fault above the aquifer. here the aquifer feeds the fault
    # instead, so the flow in the fault above the aquifer has to be kept
    aquifer_cuts_off_fault = False

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

    # add additional mesh points for the borehole
    discretize_borehole = False

    # depth of boreholes
    # this is only used for adding mesh node points to resolve the borehole
    borehole_depths = [500.0]

    #
    borehole_cellsize = 50.0

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
    # runs. switched off here so that the output contains only the fault depths
    # of the sweep below
    initial_base_run = True

    # parameters that will be changed in the sensitivity analysis runs:
    #fault_bottoms_s = [[-300.0], [-400.0], [-450.0], [-500.0],
    #                   [-550.0], [-600.0], [-700.0], [-850.0]]

    #fault_bottoms_s = [[-3000.0], [-4000.0], [-5000.0], [-6000.0], [-7000.0], [-8000.0], [-9000.0]]

    #dt_s = [10 * year]

    #cellsize_fault_s = [10.0, 5.0, 2.5]
