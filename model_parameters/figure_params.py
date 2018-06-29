"""

"""

__author__ = 'Elco Luijendijk'

import matplotlib.cm

# size of plot in inches
xsize = 10.0
ysize = 6.0

# resolution in dpi
figure_resolution = 150

# figure format, png, jpg, svg, eps, etc..
fig_format = 'png'

# ratios of the different rows in the figure, relative sizes (numbers dont matter only the differences
# [top panel, main panel, space between main panel and colorbar, colorbar]
# adjust this to move the legend upwards, reduce the empty space between the colorbar and main panel, etc...
height_ratios = [25, 75, 7, 3]

# limits of x and y axis of figure, set xlim = None and ylim = None to show the entire model domain
#xlim = None
#ylim = None
xlim = [-200, 200]
ylim = [-300, 100]

# add an extra panel on the right with the modeled temperatures
add_temperature_panel = False
relative_size_temp_panel = 0.5

# grid size for interpolating model results
dx = 5.0
dy = 2.0

# timeslices to show
# timeslices = [0] make a figure of the first timestep
# timeslices = [-1] make a figure of the last timestep
timeslices = [1, -1]

# scale for flow arrows, usually 5 is a good value, note higher value is smaller arrows (!)
scale_multiplier = 15.0

# show only one in every x flow arrows:
skip_arrows = 5

arrow_transparency = 0.75

# option to show nodes where water vapour is present
show_vapour = True

# marker size and color for AHe data:
marker_size = 4
AHe_color = 'blue'
#AHe_color = 'None'
show_AHe_data = True


# add a legend or not
add_legend = True
legend_font_size = 'medium'

# interval temperature contours
cnt_int = 10.0

# min and max of temperature contours
# set to cnt_range = None to use the min and max values of the data
#cnt_range = [10.0, 200.0]
cnt_range = None

# option to show corrected or uncorrected AHe ages
show_corrected_AHe_ages = True

# option to show mesh nodes
show_mesh = False

# linestyles for the modeled surface/borehole temperature
#linestyles = [':', '--', '-', '-.']
linestyles = ['-', '-', '-', '-']

# colors for the modeled surface/borehole temperature
#colors = ['green', 'orange', 'red', 'yellow']
colors = ['black', 'black', 'black', 'black']

# colors for the temperature contours:
cmap = matplotlib.cm.get_cmap('coolwarm')
