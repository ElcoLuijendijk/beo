"""
script to help prepare a geological cross-section

"""


__author__ = 'elco'

import numpy as np


def intersectLinePoint(xy_line, xy_point):

    """

    Parameters:
    -----------


    Returns:
    --------

    """

    ((x1, y1), (x2, y2)) = xy_line
    (xp, yp) = xy_point


    # slope of line
    ul = (y2-y1) / (x2-x1)
    # slope of normal intersecting line
    um = -1.0 / ul
    # coordinates of intersection of normal and line
    xint = (y1-yp + um*xp - ul*x1) / (um-ul)
    yint = yp + um*(xint-xp)

    # check if intersection on line:
    xs = np.array([x1, x2])
    ys = np.array([y1, y2])
    if xint < xs.min() or xint > xs.max() or yint < ys.min() or yint > ys.max():
        onLine = False
    else:
        onLine = True

    return (xint, yint), onLine


xy_line = [(0.0, 0.0),
           (1000.0, 1000.0)]

xy_point = [1000.0, 0.0]

xy_int, onLine = intersectLinePoint(xy_line, xy_point)

print 'the coordinates of the intersection point of line ', xy_line
print 'and point ', xy_point
print 'is: ', xy_int

if onLine is True:
    print 'this point is located on the line'
else:
    print 'note, this point is not located on the line'

print 'done'