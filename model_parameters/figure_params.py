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
xlim = [0, 6000]
ylim = [-6000, 100]

# grid size for interpolating model results
dx = 10.0
dy = 10.0

# timeslices to show
# timeslices = [0] make a figure of the first timestep
# timeslices = [-1] make a figure of the last timestep
timeslices = [1, -5, -1]

# scale for flow arrows, usually 5 is a good value
scale_multiplier = 5.0

# option to show nodes where water vapour is present
show_vapour = True

# marker size and color for AHe data:
marker_size = 4
AHe_color = 'blue'

# add a legend or not
add_legend = True
legend_font_size = 'x-small'