"""
parameter ranges for sensitivity or uncertainty analysis

beo.py will look for any variable ending with _s below and then look for the
corresponding variable in model_parameters.py

each _s variable should be a list of values, beo.py will replace the variable
in model_parameters.py with each item in the list consecutively

"""

import numpy as np

__author__ = 'Elco Luijendijk'


year = 365.25 * 24 * 60 * 60.0

# option whether to vary one model parameter at a time
# (ie for a sensitivtiy analysis)
# or to run all parameter combinations, using the parameter ranges specified
# in the parameter_ranges.py file
parameter_combinations = False

# option to add a first base run with unchanged parameters to the lsit of model
# runs
initial_base_run = False

###################################################################
# parameters that will be changes in the sensitivity analysis runs:
###################################################################

#fault_bottoms_s = [[-5000], [-4000]]

#thermal_gradient_s = [0.04]

#exhumation_rate_s = [1.0e-4, 1.0e-5]

#fault_fluxes_s = [[[-200.0 / year]], [[-600.0 / year]]]

#aquifer_fluxes_s = [[[-250.0 / year]]]

#fault_widths_s = [[30.0], [40.0]]

#total_depth_s = [8000.0]


# low exhumation rates (< 1e-4) result in solver errors.
# Not yet sure why. Potentially there are problems with the grid, with low exhumation rates the layers between the
# initial surface and final surface get really thin, which may result with grid elements that have angles that are
# too low. reducing exhumation_steps may help, or reducing the grid cell size. Reducing the time step size may also
# help with unstable model runs in general
#exhumation_rate_s = [1e-4]

#radius_s = [60e-6, 100e-6, 150e-6]



dt_s = [500 * year]

#cellsize_fault_s = [2.5]
#cellsize_s = [400.0, 200.0, 100.0]

#dt_s = [1000 * year, 500 * year, 250 * year, 100 * year]
