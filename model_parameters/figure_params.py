"""

"""

__author__ = 'Elco Luijendijk'


# size of plot in inches
xsize = 10.0
ysize = 6.0

# resolution in dpi
figure_resolution = 150

# figure format, png, jpg, svg, eps, etc..
fig_format = 'png'

# limits of x and y axis of figure
xlim = [4750, 5250]
ylim = [-500, 100]

# grid size for interpolating model results
dx = 5.0
dy = 5.0

# timeslices to show
# timeslices = [0] make a figure of the first timestep
# timeslices = [-1] make a figure of the last timestep
timeslices = [-5, -1]

# scale for flow arrows, usually 5 is a good value
scale_multiplier = 5.0

# show only one in every x flow arrows:
skip_arrows = 40

arrow_transparency = 0.75

# option to show nodes where water vapour is present
show_vapour = True

# marker size and color for AHe data:
marker_size = 4
AHe_color = 'blue'

# add a legend or not
add_legend = True
legend_font_size = 'x-small'

# interval temperature contours
cnt_int = 10.0

# min and max of temperature contours
# set to cnt_range = None to use the min and max values of the data
#cnt_range = [10.0, 200.0]
cnt_range = None

# option to show corrected or uncorrected AHe ages
show_corrected_AHe_ages = True