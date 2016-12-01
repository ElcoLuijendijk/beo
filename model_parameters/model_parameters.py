"""
file containing all parameters for hydrotherm escript

"""

__author__ = 'elco'

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

# model dimensions
width = 10000.0
total_depth = 6000.0
air_height = 100.0

# depth to fine discretization at surface:
z_fine = -200

# grid size
cellsize = 500.0
cellsize_air = 10.0
cellsize_fault = 5.0
cellsize_fine = 25.0
cellsize_base = 1000.0

# temperature bnd conditions
air_temperature = 10.0
#bottom_temperature = total_depth * 0.03 + air_temperature
#bottom_temperature = total_depth * 0.03
# new version: calculate bottom T using a fixed geothermal gradient./r
thermal_gradient = 0.03

# bottom flux bnd condition, set to None if T bnd is used
basal_heat_flux = None
#basal_heat_flux = 65e-3

porosity = 0.15

# thermal parameters

# wikipedia: thermal properties air
# heat transfer coefficient = 10- 100 W / (m2 K))
# heat capacity = 1000 J kg-1 K-1
# density = 1.29 kg m-3
K_solid = 3.5
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
N_outputs = [11]
# size of timestep
dt = 1000 * year

# duration of each timestep
durations_all = [[2e5 * year]]

# target depth slices for calculating temperature and U-Th/He
target_zs = [0, -200]

# U-Th/He params
calculate_he_ages = True

T0 = 10.0
T_surface = 10.0
t0 = 100.0 * My
radius = 60.0 * 1e-6
U238 = 8.98e-6
Th232 = 161.3e-6

## fault data for multiple faults:

# x location of fault:
fault_xs = [2000]

# fault width
fault_widths = [20.0]

# angle of the fault zone (degrees), dip of normal faults ~60-70 degrees
# but start with vertical fault to keep things simple
fault_angles = [65.0]

# elevation of bottom of fault
fault_bottoms = [-5000.0]

# fluid advection rates in faults:
# nested list,
# [[fault1_t1, fault2_t1], [fault1_t2, fault2_t2], etc...]
fault_fluxes_all = [[[9.0 / year]]]

aquifer_bottoms = [None]
aquifer_tops = [None]
aquifer_fluxes = [None]

############################################
# variable params, for sensitivity analysis
###########################################