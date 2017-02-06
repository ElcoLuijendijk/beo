"""
parameter ranges for sensitivity or uncertainty analysis

beo.py will look for any variable ending with _s below and then look for the
corresponding variable in model_parameters.py

each _s variable should be a list of values, beo.py will replace the variable
in model_parameters.py with each item in the list consecutively

"""

import numpy as np

__author__ = 'elcopone'


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

#exhumation_rate_s = [1.0e-5, 1.0e-4]

fault_fluxes_s = [[[10.0 / year]], [[15.0 / year]]]