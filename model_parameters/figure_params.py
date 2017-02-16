__author__ = 'elco'


# size of plot in inches
xsize = 10.0
ysize = 6.0

# limits of x and y axis of figure
xlim = [3500, 4500]
ylim = [-1000, 100]

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