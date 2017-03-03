"""
2D model of advective and conductive heat flow in hydrothermal systems

Elco Luijendijk, McGill university & Goettingen Unviersity, 2013-2017
"""

__author__ = 'Elco Luijendijk'

###############
# load modules
###############

import time
import os
import sys
import pickle
import pdb
import datetime

import numpy as np
import matplotlib

# escript/Finley modules:
import esys.escript as es
import esys.finley as fl
import esys.pycad as pc
import esys.pycad.gmsh as gmsh
import esys.escript.linearPDEs as linearPDEs

# helium diffusion algorithm by Meesters and Dunai (2003)
import lib.helium_diffusion_models as he


def calculate_vapour_pressure(T,
                              c1=8e-8,
                              c2=-7e-5,
                              c3=0.028,
                              c4=-3.1597):

    """
    calculate the vapour pressure curve and check whether there should be
    vapour or not in the model domain

    based on a 3rd order polyonmial fit to vapour curve data by the NIST
    publication "Thermophysical Properties of Fluid Systems"
    found at http://webbook.nist.gov/chemistry/fluid/
    """

    log_Pv = c1 * T**3 + c2 * T**2 + c3 * T + c4

    Pv = 10**log_Pv * 1e6

    return Pv


def calculate_boiling_temp(P,
                      c1=3.866,
                      c2=25.151,
                      c3=103.28,
                      c4=179.99):
    """
    find the maximum temperature for a given pressure at which there is one
    liquid phase only

    based on a 3rd order polyonmial fit to vapour curve data by the NIST
    publication "Thermophysical Properties of Fluid Systems"
    found at http://webbook.nist.gov/chemistry/fluid/

    input pressure in Pa
    returns T in degrees C

    """

    logP = es.log10(P / 1.0e6)

    Tmax = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4

    return Tmax


def convert_to_array(u, no_coords=False):

    """
    return the x,y coordinates and the value of escript variable u as
    a numpy array
    """

    coords = u.getFunctionSpace().getX()
    x, y = coords[0], coords[1]

    xy = np.array([x.toListOfTuples(), y.toListOfTuples()]).T

    u_array = np.array(u.toListOfTuples())

    assert len(u_array.shape) == 1

    if no_coords is True:
        return u_array
    else:
        return xy, u_array


def interpolate_data(xyz_array, Ti, dx, dy):

    """

    """
    xi = np.arange(xyz_array[:, 0].min(), xyz_array[:, 0].max() + dx, dx)
    yi = np.arange(xyz_array[:, 1].min(), xyz_array[:, 1].max() + dy, dy)
    xg, yg = np.meshgrid(xi, yi)
    #xgf, ygf = xg.flatten(), yg.flatten()
    #zgf = scipy.interpolate.griddata(xyz_array, Ti, np.vstack((xgf, ygf)).T,
    #                                 method='linear')
    zg = matplotlib.mlab.griddata(xyz_array[:, 0], xyz_array[:, 1], Ti, xi, yi)

    return xg, yg, zg


def calculate_fault_x(z_flt, fault_angle, x_flt_surface):

    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    return x_flt


def setup_mesh(width, x_flt_surface, fault_width, fault_angle, z_air,
               z_surface, z_fine, z_base, cellsize,
               cellsize_air, cellsize_fault, cellsize_fine, cellsize_base):

    """

    create a mesh for the model of the bewowawe hydrothermal system

    """

    ###############################
    # use gmsh to construct domain
    ##############################

    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain
    z_flt = np.array([z_surface, z_fine, z_base])
    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    print 'fault locations in mesh:'
    for x, z in zip(x_flt, z_flt):
        print x, z

    xys = [[0, z_air], [width, z_air],
           [0, z_surface], [x_flt[0], z_surface], [x_flt[0] + fault_width, z_surface], [width, z_surface],
           [0, z_fine], [x_flt[1], z_fine], [x_flt[1] + fault_width, z_fine], [width, z_fine],
           [0, z_base], [x_flt[2], z_base], [x_flt[2] + fault_width, z_base], [width, z_base]]

    #points = create_points(xs,zs)
    points = [pc.Point(x, z) for x, z in xys]

    # coarse cellsize at lower right
    points[0].setLocalScale(cellsize_air / cellsize)
    points[1].setLocalScale(cellsize_air / cellsize)
    #points[7].setLocalScale(cellsize_fine / cellsize)

    points[3].setLocalScale(cellsize_fault / cellsize)
    points[4].setLocalScale(cellsize_fault / cellsize)
    points[7].setLocalScale(cellsize_fault / cellsize)
    points[8].setLocalScale(cellsize_fault / cellsize)
    points[11].setLocalScale(cellsize_fault / cellsize)
    points[12].setLocalScale(cellsize_fault / cellsize)

    #points[5].setLocalScale(cellsize_fine / cellsize)
    #points[8].setLocalScale(cellsize_fine / cellsize)
    points[10].setLocalScale(cellsize_base / cellsize)
    points[13].setLocalScale(cellsize_base / cellsize)

    # horizontal lines:
    lineh1 = pc.Line(points[0], points[1])
    lineh2 = pc.Line(points[2], points[3])
    lineh3 = pc.Line(points[3], points[4])
    lineh4 = pc.Line(points[4], points[5])
    lineh5 = pc.Line(points[6], points[7])
    lineh6 = pc.Line(points[7], points[8])
    lineh7 = pc.Line(points[8], points[9])
    lineh8 = pc.Line(points[10], points[11])
    lineh9 = pc.Line(points[11], points[12])
    lineh10 = pc.Line(points[12], points[13])

    # vertical lines:
    linev1 = pc.Line(points[0], points[2])
    linev2 = pc.Line(points[1], points[5])
    linev3 = pc.Line(points[2], points[6])
    linev4 = pc.Line(points[3], points[7])
    linev5 = pc.Line(points[4], points[8])
    linev6 = pc.Line(points[5], points[9])
    linev7 = pc.Line(points[6], points[10])
    linev8 = pc.Line(points[7], points[11])
    linev9 = pc.Line(points[8], points[12])
    linev10 = pc.Line(points[9], points[13])

    # closed curves for different segments of the mesh:
    curve_air = pc.CurveLoop(lineh1, linev2, -lineh4, -lineh3, -lineh2, -linev1)
    curve_fine_left = pc.CurveLoop(lineh2, linev4, -lineh5, -linev3)
    curve_flt_fine = pc.CurveLoop(lineh3, linev5, -lineh6, -linev4)
    curve_fine_right = pc.CurveLoop(lineh4, linev6, -lineh7, -linev5)
    curve_base_left = pc.CurveLoop(lineh5, linev8, -lineh8, -linev7)
    curve_flt_base = pc.CurveLoop(lineh6, linev9, -lineh9, -linev8)
    curve_base_right = pc.CurveLoop(lineh7, linev10, -lineh10, -linev9)

    surface_air = pc.PlaneSurface(curve_air)
    surface_flt_fine = pc.PlaneSurface(curve_flt_fine)
    surface_fine_left = pc.PlaneSurface(curve_fine_left)
    surface_fine_right = pc.PlaneSurface(curve_fine_right)
    surface_flt_base = pc.PlaneSurface(curve_flt_base)
    surface_base_left = pc.PlaneSurface(curve_base_left)
    surface_base_right = pc.PlaneSurface(curve_base_right)

    #surface_air.setLocalScale(factor=cellsize_air / cellsize)
    #surface_flt_fine.setLocalScale(factor=cellsize_fault / cellsize)
    #surface_flt_base.setLocalScale(factor=cellsize_fault / cellsize)
    #surface_fine.setLocalScale(factor=cellsize_fine / cellsize)
    #surface_base.setLocalScale(factor=cellsize_base / cellsize)
    #surface_c.setLocalScale(factor=Parameters.grid_refinement_factor)

    #if fine_mesh is True:
    #    print 'assigning refined grid to entire landward side of model domain'
    #    surface_d.setLocalScale(factor=Parameters.grid_refinement_factor)

    d = gmsh.Design(dim=2, element_size=cellsize)

    d.setMeshFileName('beowawe_mesh')
    ps1 = pc.PropertySet("bottomleft", surface_base_left)
    ps2 = pc.PropertySet("bottommid", surface_flt_base)
    ps3 = pc.PropertySet("bottomright", surface_base_right)

    d.addItems(surface_air, surface_flt_fine,
               surface_fine_left, surface_fine_right,
               ps1, ps2, ps3)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    mesh.write('mesh.fly')

    return mesh


def setup_mesh_with_exhumation(width, x_flt_surface, fault_width, fault_angle,
                               z_air,
                               z_surface_initial, z_surface_final,
                               z_surface_steps,
                               z_fine, z_base, cellsize,
                               cellsize_air, cellsize_fault,
                               cellsize_fine, cellsize_base):

    """

    create a mesh for the model of the bewowawe hydrothermal system

    """

    ###############################
    # use gmsh to construct domain
    ##############################

    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain
    zs_surface = np.linspace(z_surface_initial, z_surface_final, z_surface_steps + 1)
    z_flt = np.concatenate((zs_surface, np.array([z_fine, z_base])))
    #z_flt = np.array([z_surface, z_fine, z_base])
    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    print 'fault locations in mesh:'
    for x, z in zip(x_flt, z_flt):
        print x, z

    #xys = [[0, z_air], [width, z_air],
    #       [0, z_surface], [x_flt[0], z_surface], [x_flt[0] + fault_width, z_surface], [width, z_surface],
    #       [0, z_fine], [x_flt[1], z_fine], [x_flt[1] + fault_width, z_fine], [width, z_fine],
    #       [0, z_base], [x_flt[2], z_base], [x_flt[2] + fault_width, z_base], [width, z_base]]

    xys = [[[0, z_air], [width, z_air]]]

    for xf, zf in zip(x_flt, z_flt):
        xys.append([[0, zf], [xf, zf], [xf + fault_width, zf], [width, zf]])

    points = []
    for xyi in xys:
        points_local = [pc.Point(x, z) for x, z in xyi]
        points.append(points_local)

    # fine cellsize in air layer:
    for point in points[0]:
        point.setLocalScale(cellsize_air / cellsize)

    # small cellsize in fault:
    for point in points[1:]:
        point[1].setLocalScale(cellsize_fault / cellsize)
        point[2].setLocalScale(cellsize_fault / cellsize)

    points[-1][0].setLocalScale(cellsize_base / cellsize)
    points[-1][-1].setLocalScale(cellsize_base / cellsize)

    # horizontal lines:
    hlines = [[pc.Line(points[0][0], points[0][1])]]
    for point in points[1:]:
        hline_local = [pc.Line(point[0], point[1]),
                       pc.Line(point[1], point[2]),
                       pc.Line(point[2], point[3])]
        hlines.append(hline_local)

    # vertical lines:
    vlines = [[pc.Line(points[0][0], points[1][0]), pc.Line(points[0][1], points[1][3])]]
    for point, point_below in zip(points[1:], points[2:]):
        vline_local = [pc.Line(point[0], point_below[0]),
                       pc.Line(point[1], point_below[1]),
                       pc.Line(point[2], point_below[2]),
                       pc.Line(point[3], point_below[3])]
        vlines.append(vline_local)

    curves = [pc.CurveLoop(hlines[0][0], vlines[0][1],
                           -hlines[1][2], -hlines[1][1], -hlines[1][0],
                           -vlines[0][0])]
    for hline, hline_below, vline in zip(hlines[1:-1], hlines[2:], vlines[1:]):
        curve_local_left = pc.CurveLoop(hline[0], vline[1], -hline_below[0], -vline[0])
        curve_local_fault = pc.CurveLoop(hline[1], vline[2], -hline_below[1], -vline[1])
        curve_local_right = pc.CurveLoop(hline[2], vline[3], -hline_below[2], -vline[2])

        curves += [curve_local_left, curve_local_fault, curve_local_right]

    surfaces = [pc.PlaneSurface(curve) for curve in curves]

    d = gmsh.Design(dim=2, element_size=cellsize)

    d.setMeshFileName('beowawe_mesh')

    d.addItems(*surfaces)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    mesh.write('mesh.fly')

    return mesh


def setup_mesh_with_exhumation_v2(width, x_flt_surface, fault_width, fault_angle,
                                  z_air,
                                  z_surface_initial, z_surface_final,
                                  z_surface_steps,
                                  z_fine, z_base, cellsize,
                                  cellsize_air, cellsize_fault,
                                  cellsize_fine, cellsize_base, fault_buffer_zone):

    """

    create a mesh for the model of the bewowawe hydrothermal system

    """

    ###############################
    # use gmsh to construct domain
    ##############################

    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain
    zs_surface = np.linspace(z_surface_initial, z_surface_final, z_surface_steps + 1)
    z_flt = np.concatenate((zs_surface, np.array([z_fine, z_base])))
    #z_flt = np.array([z_surface, z_fine, z_base])
    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    print 'fault locations in mesh:'
    for x, z in zip(x_flt, z_flt):
        print x, z

    #xys = [[0, z_air], [width, z_air],
    #       [0, z_surface], [x_flt[0], z_surface], [x_flt[0] + fault_width, z_surface], [width, z_surface],
    #       [0, z_fine], [x_flt[1], z_fine], [x_flt[1] + fault_width, z_fine], [width, z_fine],
    #       [0, z_base], [x_flt[2], z_base], [x_flt[2] + fault_width, z_base], [width, z_base]]

    xys = [[[0, z_air], [width, z_air]]]

    for xf, zf in zip(x_flt, z_flt):
        xys.append([[0, zf], [xf - fault_buffer_zone, zf], [xf, zf], [xf + fault_width, zf], [xf + fault_buffer_zone * 2, zf],  [width, zf]])

    points = []
    for xyi in xys:
        points_local = [pc.Point(x, z) for x, z in xyi]
        points.append(points_local)

    # fine cellsize in air layer:
    for point in points[0]:
        point.setLocalScale(cellsize_air / cellsize)

    # small cellsize in fault:
    for point in points[1:]:
        point[1].setLocalScale(cellsize_fine / cellsize)
        point[2].setLocalScale(cellsize_fault / cellsize)
        point[3].setLocalScale(cellsize_fault / cellsize)
        point[4].setLocalScale(cellsize_fine / cellsize)

    points[-1][0].setLocalScale(cellsize_base / cellsize)
    points[-1][-1].setLocalScale(cellsize_base / cellsize)

    # horizontal lines:
    hlines = [[pc.Line(points[0][0], points[0][1])]]
    for point in points[1:]:
        hline_local = [pc.Line(point[0], point[1]),
                       pc.Line(point[1], point[2]),
                       pc.Line(point[2], point[3]),
                       pc.Line(point[3], point[4]),
                       pc.Line(point[4], point[5])]
        hlines.append(hline_local)

    # vertical lines:
    vlines = [[pc.Line(points[0][0], points[1][0]), pc.Line(points[0][1], points[1][3])]]
    for point, point_below in zip(points[1:], points[2:]):
        vline_local = [pc.Line(point[0], point_below[0]),
                       pc.Line(point[1], point_below[1]),
                       pc.Line(point[2], point_below[2]),
                       pc.Line(point[3], point_below[3]),
                       pc.Line(point[4], point_below[4]),
                       pc.Line(point[5], point_below[5])]
        vlines.append(vline_local)

    curves = [pc.CurveLoop(hlines[0][0], vlines[0][1],
                           -hlines[1][4], -hlines[1][3], -hlines[1][2], -hlines[1][1], -hlines[1][0],
                           -vlines[0][0])]
    for hline, hline_below, vline in zip(hlines[1:-1], hlines[2:], vlines[1:]):
        curve_local_left = pc.CurveLoop(hline[0], vline[1], -hline_below[0], -vline[0])
        curve_local_fault_buffer_left = pc.CurveLoop(hline[1], vline[2], -hline_below[1], -vline[1])
        curve_local_fault = pc.CurveLoop(hline[2], vline[3], -hline_below[2], -vline[2])
        curve_local_fault_buffer_right = pc.CurveLoop(hline[3], vline[4], -hline_below[3], -vline[3])
        curve_local_right = pc.CurveLoop(hline[4], vline[5], -hline_below[4], -vline[4])

        curves += [curve_local_left, curve_local_fault_buffer_left, curve_local_fault, curve_local_fault_buffer_right, curve_local_right]

    surfaces = [pc.PlaneSurface(curve) for curve in curves]

    d = gmsh.Design(dim=2, element_size=cellsize)

    d.setMeshFileName('beowawe_mesh')

    d.addItems(*surfaces)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    mesh.write('mesh.fly')

    return mesh


def setup_mesh_2faults(width, x_flt_surface, fault_width, fault_angle, z_air,
                       z_surface, z_fine, z_base, cellsize,
                       cellsize_air, cellsize_fault, cellsize_fine, cellsize_base):

    #mp.fault_xs[0], mp.fault_widths[0], mp.fault_angles[0], mp.air_height,
    #              z_surface, mp.z_fine, z_base, mp.cellsize,
    #              mp.cellsize_air, mp.cellsize_fault,
    #              mp.cellsize_fine, mp.cellsize_base


    """
    create a mesh for a hydrothermal model containing 2 faults
    and one shallow aquifer

    """

    ###############################
    # use gmsh to construct domain
    ##############################

    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain
    z_flt = np.array([z_surface, z_fine, z_base])
    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    xys = [[0, z_air], [width, z_air],
           [x_flt[0], z_surface], [x_flt[0] + fault_width, z_surface], [width, z_surface],
           [0, z_fine], [x_flt[1], z_fine], [x_flt[1] + fault_width, z_fine], [width, z_fine],
           [0, z_base], [x_flt[2], z_base], [x_flt[2] + fault_width, z_base], [width, z_base]]

    #points = create_points(xs,zs)
    points = [pc.Point(x, z) for x, z in xys]

    # coarse cellsize at lower right
    points[0].setLocalScale(cellsize_air / cellsize)
    points[1].setLocalScale(cellsize_air / cellsize)
    #points[7].setLocalScale(cellsize_fine / cellsize)

    points[2].setLocalScale(cellsize_fault / cellsize)
    points[3].setLocalScale(cellsize_fault / cellsize)
    points[6].setLocalScale(cellsize_fault / cellsize)
    points[7].setLocalScale(cellsize_fault / cellsize)
    points[10].setLocalScale(cellsize_fault / cellsize)
    points[11].setLocalScale(cellsize_fault / cellsize)

    #points[5].setLocalScale(cellsize_fine / cellsize)
    #points[8].setLocalScale(cellsize_fine / cellsize)
    points[9].setLocalScale(cellsize_base / cellsize)
    points[12].setLocalScale(cellsize_base / cellsize)

    lineh1 = pc.Line(points[0], points[1])
    lineh2 = pc.Line(points[2], points[3])
    lineh3 = pc.Line(points[3], points[4])
    lineh4 = pc.Line(points[5], points[6])
    lineh5 = pc.Line(points[6], points[7])
    lineh6 = pc.Line(points[7], points[8])
    lineh7 = pc.Line(points[9], points[10])
    lineh8 = pc.Line(points[10], points[11])
    lineh9 = pc.Line(points[11], points[12])
    linev1 = pc.Line(points[0], points[2])
    linev2 = pc.Line(points[1], points[4])
    linev3 = pc.Line(points[2], points[5])
    linev4 = pc.Line(points[2], points[6])
    linev5 = pc.Line(points[3], points[7])
    linev6 = pc.Line(points[4], points[8])
    linev7 = pc.Line(points[5], points[9])
    linev8 = pc.Line(points[6], points[10])
    linev9 = pc.Line(points[7], points[11])
    linev10 = pc.Line(points[8], points[12])

    # finer grid cell size around fresh-salt water interface
    curve_air = pc.CurveLoop(lineh1, linev2, -lineh3, -lineh2, -linev1)
    curve_flt_fine = pc.CurveLoop(lineh2, linev5, -lineh5, -linev4)
    curve_fine_left = pc.CurveLoop(linev4, -lineh4, -linev3)
    curve_fine_right = pc.CurveLoop(lineh3, linev6, -lineh6, -linev5)
    curve_flt_base = pc.CurveLoop(lineh5, linev9, -lineh8, -linev8)
    curve_base_left = pc.CurveLoop(lineh4, linev8, -lineh7, -linev7)
    curve_base_right = pc.CurveLoop(lineh6, linev10, -lineh9, -linev9)

    surface_air = pc.PlaneSurface(curve_air)
    surface_flt_fine = pc.PlaneSurface(curve_flt_fine)
    surface_fine_left = pc.PlaneSurface(curve_fine_left)
    surface_fine_right = pc.PlaneSurface(curve_fine_right)
    surface_flt_base = pc.PlaneSurface(curve_flt_base)
    surface_base_left = pc.PlaneSurface(curve_base_left)
    surface_base_right = pc.PlaneSurface(curve_base_right)

    #surface_air.setLocalScale(factor=cellsize_air / cellsize)
    #surface_flt_fine.setLocalScale(factor=cellsize_fault / cellsize)
    #surface_flt_base.setLocalScale(factor=cellsize_fault / cellsize)
    #surface_fine.setLocalScale(factor=cellsize_fine / cellsize)
    #surface_base.setLocalScale(factor=cellsize_base / cellsize)
    #surface_c.setLocalScale(factor=Parameters.grid_refinement_factor)


    #if fine_mesh is True:
    #    print 'assigning refined grid to entire landward side of model domain'
    #    surface_d.setLocalScale(factor=Parameters.grid_refinement_factor)

    d = gmsh.Design(dim=2, element_size=cellsize)

    d.setMeshFileName('beowawe_mesh')

    ps1 = es.PropertySet("bottomleft", surface_base_left)
    ps2 = es.PropertySet("bottommid", surface_flt_base)
    ps3 = es.PropertySet("bottomright", surface_base_right)


    d.addItems(surface_air, surface_flt_fine, surface_flt_base,
               surface_fine_left, surface_fine_right,
               surface_base_left, surface_base_right, ps1, ps2, ps3)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    mesh.write('mesh.fly')

    return mesh


def interpolate_var_to_2d_grid(model_var):

    """
    convert escript variable to 2d numpy array
    only works for regular 2d grids
    """

    xyz, model_array = convert_to_array(model_var)

    # interpolate variable back to 2d grid
    xi, xi_pos = np.unique(xyz[:, 0], return_inverse=True)
    yi, yi_pos = np.unique(xyz[:, 1], return_inverse=True)

    # resize doesnt work, escript screws up row and column order of data...
    #raster_array = np.resize(model_array, [len(xi), len(yi)])
    raster_array = np.zeros((len(xi), len(yi)))
    raster_array[xi_pos, yi_pos] = model_array

    return xi, yi, raster_array


def model_hydrothermal_temperatures(mesh, hf_pde,
                                    fault_zones, fault_angles,
                                    specified_T_loc, specified_flux_loc,
                                    durations, fault_fluxes,
                                    K_var, rho_var, c_var,
                                    rho_f, c_f,
                                    specified_temperature, specified_flux,
                                    dt,
                                    top_bnd, bottom_bnd,
                                    air_temperature, bottom_temperature,
                                    solve_as_steady_state=True,
                                    surface_level_init=0, exhumation_rate=0,
                                    target_depths=None,
                                    K_b=None, c_b=None, rho_b=None,
                                    K_air=None, c_air=None, rho_air=None,
                                    vapour_correction=True):

    """

    """
    
    day = 24.0 * 60.0 * 60.0
    year = 365.25 * day

    fluid_density = 1000.0
    g = 9.81
    atmospheric_pressure = 101325

    c1 = 3.866
    c2 = 25.151
    c3 = 103.28
    c4 = 179.99

    P_buffer = 10.0

    xyz = mesh.getX()

    exceed_max_liquid_T_old = None

    ######################################
    # model steady-state temperature field
    ######################################
    # set PDE coefficients, steady-state heat flow equation
    A = K_var * es.kronecker(mesh)
    C = 0
    D = 0

    if specified_flux is not None:
        print 'solving with specified heat flux bnd'
        #specified_heat_flux = specified_flux * specified_flux_loc
        specified_heat_flux = specified_flux

        pdb.set_trace()

        hf_pde.setValue(A=A, D=D,
                        r=specified_temperature,
                        q=specified_T_loc,
                        y=specified_heat_flux)
    else:
        print 'only fixed T bnd'
        hf_pde.setValue(A=A, D=D,
                        r=specified_temperature,
                        q=specified_T_loc)

    T_steady = hf_pde.getSolution()
    T = T_steady

    print 'modeled steady state temperatures ', T_steady

    print 'starting transient heat flow calculations'
    runtimes = []
    t_total = 0
    q_vectors = []
    Ts = []
    surface_levels = []
    surface_level = surface_level_init

    depth = -(xyz[1] - surface_level)
    subsurface = es.whereNonPositive(xyz[1] - surface_level)
    air = es.wherePositive(xyz[1] - surface_level)
    surface = es.whereZero(xyz[1] - surface_level)

    if vapour_correction is True:
        P_init = fluid_density * g * depth + atmospheric_pressure
        P = P_init * es.wherePositive(depth) + atmospheric_pressure * es.whereNonPositive(depth)

        vapour_pressure = calculate_vapour_pressure(T)
        logP = es.log10(P / 1.0e6)
        boiling_temp = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4
        vapour = es.whereNegative(P - vapour_pressure)
        xmin_vapour = es.inf(vapour * xyz[0])
        xmax_vapour = es.sup(vapour * xyz[0])
        ymin_vapour = -es.sup(-(vapour * xyz[1]))
        ymax_vapour = es.sup(vapour * xyz[1])

        boiling_temps = []
        exceed_boiling_temps = []

        if es.sup(vapour) >= 1:
            print 'warning, vapour present at initial steady-state P-T conditions'

    for fault_flux, duration in zip(fault_fluxes, durations):

        print 'duration of %i' % duration
        print 'fluxes in faults: ', fault_flux

        # set up advective flow field in faults
        q_vector = es.Vector((0, 0), es.Function(mesh))

        for fault_zone, fault_angle, q_fault_zone in \
                zip(fault_zones, fault_angles, fault_flux):

            # add heat advection in fault zone
            qh_fault_zone = - q_fault_zone * np.cos(np.deg2rad(fault_angle))
            qv_fault_zone = q_fault_zone * np.sin(np.deg2rad(fault_angle))
            q_vector[0] += fault_zone * qh_fault_zone
            q_vector[1] += fault_zone * qv_fault_zone

        ###############################################
        # model transient response to fluid convection
        ###############################################
        if solve_as_steady_state is False:
            # set PDE coefficients, transient heat flow equation
            A = dt * K_var * es.kronecker(mesh)
            C = dt * rho_f * c_f * q_vector
            D = rho_var * c_var
            Y = rho_var * c_var * T

        else:
            print 'solving steady-state, with advective flux'

            # set PDE coefficients, steady-state flow equation
            A = K_var * es.kronecker(mesh)
            C = rho_f * c_f * q_vector
            D = 0
            Y = 0

        # update bnd cond if spec flux bnd
        if specified_flux is not None:

            print 'adding specified flux bnd'
            specified_heat_flux = specified_flux * specified_flux_loc * dt

            #
            hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                            r=specified_temperature,
                            q=specified_T_loc,
                            y=specified_heat_flux)
        else:
            hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                            r=specified_temperature,
                            q=specified_T_loc)

        # calculate courant number
        cfl_cond = q_vector * dt / mesh.getSize()
        cfl_cond_max = es.sup(cfl_cond)
        print 'CFL number: ', cfl_cond

        if cfl_cond_max > 1:

            print 'warning, cfl condition not met, cfl number = %0.2f' % cfl_cond_max

            max_cfl_number = 0.5

            #print 'dt was %0.2e' % dt

            #dt = max_cfl_number * dy / es.Lsup(q_vector)

            #print 'reducing dt to %0.2e' % dt

        nt = int(duration / dt)
        #nt_recovery = 1 * nt_heating

        # caclulate grid peclet number
        Pe = rho_f * c_f * es.Lsup(q_vector) * mesh.getSize() / K_var
        print 'grid peclet number = ', Pe

        #############################################
        # iterate heat flow eq.s
        #############################################
        print 'x' * 10
        print 'starting iterations'

        start = time.time()

        if solve_as_steady_state is True:
            nt = 1

        # solve transient heat flux
        for t in range(nt):

            if t / 10 == t / 10.0:

                end = time.time()
                comptime = end - start

                start = end

                # find closest nodes to the actual surface level

                try:
                    surface_level_mesh_id = np.where(target_depths <= surface_level)[0][-1]
                    surface_level_mesh = target_depths[surface_level_mesh_id]
                except:
                    surface_level_mesh = surface_level
                    print '\twarning could not find land surface nodes'

                land_surface = es.whereZero(xyz[1] - surface_level_mesh)
                print 'step %i of %i' % (t, nt)
                print '\tcomputational time for one timestep = %0.3f sec' \
                      % (comptime / 10.0)
                print '\tactual surface level ', surface_level
                print '\tclosest surface in mesh ', surface_level_mesh
                print '\ttemperature: ', T
                if es.sup(land_surface) > 0:
                    print '\tmax. temperature at land surface: ', \
                        es.sup(T * land_surface)
                else:
                    print '\tcould not find land surface nodes'

                if vapour_correction is True:
                    vapour_pressure = calculate_vapour_pressure(T)
                    vapour = subsurface * es.whereNegative(P - vapour_pressure + P_buffer)
                    #vapour = es.whereNegative(P - vapour_pressure)

                    xmin_vapour = es.inf(vapour * xyz[0] + es.whereZero(vapour) * 999999.9)
                    xmax_vapour = es.sup(vapour * xyz[0])
                    ymin_vapour = -es.sup(-(vapour * xyz[1]))
                    ymax_vapour = es.sup(vapour * xyz[1])

                    if es.sup(vapour) > 0:
                        print '\tvapour present in: ', es.integrate(vapour), ' m^2'
                        print '\t\tfrom x = ', xmin_vapour, ' to x = ', \
                            xmax_vapour
                        print '\t\tand from y = ', ymin_vapour, ' to y = ', \
                            ymax_vapour
                        print '\tmax. liquid T at the surface = ', es.sup(boiling_temp * land_surface)
                    else:
                        print '\tno vapour present'

            surface_level = surface_level_init - t_total / year * exhumation_rate

            if exhumation_rate != 0 and surface_level in target_depths:

                print 'exhumation, new surface level at %0.2f' % surface_level
                subsurface = es.whereNonPositive(xyz[1] - surface_level)
                air = es.wherePositive(xyz[1] - surface_level)
                surface = es.whereZero(xyz[1] - surface_level)

                q_vector = q_vector * subsurface

                # populate K, c and rho scalar fields
                K_var = subsurface * K_b + air * K_air
                c_var = subsurface * c_b + air * c_air
                rho_var = subsurface * rho_b + air * rho_air

                # reset heatflow PDE coefficients
                A = dt * K_var * es.kronecker(mesh)
                C = dt * rho_f * c_f * q_vector
                D = rho_var * c_var
                Y = rho_var * c_var * T

                # update bnd cond if spec flux bnd
                if specified_flux is not None:

                    print 'adding specified flux bnd'
                    specified_heat_flux = specified_flux * specified_flux_loc * dt

                    #
                    hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                    r=specified_temperature,
                                    q=specified_T_loc,
                                    y=specified_heat_flux)
                else:
                    hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                    r=specified_temperature,
                                    q=specified_T_loc)

                # recalculate fluid pressure
                xyz = mesh.getX()
                depth = -(xyz[1] - surface_level)

                if vapour_correction is True:
                    P_init = fluid_density * g * depth + atmospheric_pressure
                    P = P_init * es.wherePositive(depth) + atmospheric_pressure * es.whereNonPositive(depth)

                    logP = es.log10(P / 1.0e6)
                    boiling_temp = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4


            # recalculate vapour pressure and max liquid temperature
            if vapour_correction is True:
                #vapour_pressure = calculate_vapour_pressure(T)
                #boiling_temp = calculate_boiling_temp(P)
                #exceed_boiling_temp = subsurface * es.wherePositive(T - boiling_temp)

                if exceed_max_liquid_T_old is None:
                    exceed_boiling_temp = subsurface * es.wherePositive(T - boiling_temp)
                else:
                    exceed_boiling_temp = \
                        subsurface * es.whereZero(exceed_max_liquid_T_old) \
                        * es.wherePositive(T - boiling_temp) + subsurface * exceed_max_liquid_T_old

                exceed_max_liquid_T_old =  exceed_boiling_temp

            if vapour_correction is True:
                specified_T_loc = es.wherePositive(top_bnd) + es.wherePositive(bottom_bnd) + exceed_boiling_temp
                specified_temperature = es.wherePositive(top_bnd) * air_temperature \
                      + es.wherePositive(bottom_bnd) * bottom_temperature + exceed_boiling_temp * boiling_temp

            # solve PDE for temperature
            T = hf_pde.getSolution()

            # update PDE coefficients
            if solve_as_steady_state is False:
                Y = rho_var * c_var * T
                hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                r=specified_temperature,
                                q=specified_T_loc)

            t_total += dt

            surface_levels.append(surface_level)

            # store output
            Ts.append(T)
            q_vectors.append(q_vector)

            if vapour_correction is True:
                boiling_temps.append(boiling_temp)
                exceed_boiling_temps.append(exceed_boiling_temp)

            #ti = output_steps.index(t)
            #print 'surface T: ', T * surface

            runtimes.append(t_total)

        print 'T after advective heating ', T

    return (np.array(runtimes), T_steady, Ts, q_vectors,
            np.array(surface_levels), boiling_temps, exceed_boiling_temps)


def model_run(mp):

    """
    setup and run a single model

    :param mp:
    :return:
    """

    #
    year = 365.25 * 24 * 60 * 60
    Myr = year * 1e6

    ############################
    # construct rectangular mesh
    ############################
    print 'constructing mesh'
    #mesh = fl.Rectangle(l0=width, l1=height, n0=nx, n1=ny)

    z_surface = 0
    z_base = -mp.total_depth
    #mesh = setup_mesh(mp.width, mp.fault_xs[0], mp.fault_widths[0],
    #                  mp.fault_angles[0], mp.air_height,
    #                  z_surface, mp.z_fine, z_base, mp.cellsize,
   #                   mp.cellsize_air, mp.cellsize_fault,
    #                  mp.cellsize_fine, mp.cellsize_base)

    exhumed_thickness = mp.exhumation_rate * (np.sum(np.array(mp.durations)) / mp.year)
    exhumation_steps = mp.exhumation_steps

    min_layer_thickness = 1.0
    if exhumed_thickness / exhumation_steps < min_layer_thickness:
        print 'warning, exhumation levels would be smaller than %0.2f m' % min_layer_thickness
        exhumation_steps = int(np.ceil(exhumed_thickness) / min_layer_thickness)
        if exhumation_steps < 1:
            exhumation_steps = 1

        print 'reducing exhumation steps to %i' % exhumation_steps

        mp.exhumation_steps = exhumation_steps

    if exhumed_thickness != 0:
        # track AHe and temperature in each exhumed layer in the model domain:
        mp.target_zs = np.linspace(0, exhumed_thickness, exhumation_steps + 1)

    elevation_top = z_surface + exhumed_thickness + mp.air_height

    if mp.use_mesh_with_buffer is False:
        mesh = setup_mesh_with_exhumation(mp.width, mp.fault_xs[0],
                                             mp.fault_widths[0],
                                             mp.fault_angles[0], elevation_top,
                                             z_surface + exhumed_thickness, z_surface,
                                             exhumation_steps,
                                             mp.z_fine, z_base, mp.cellsize,
                                             mp.cellsize_air, mp.cellsize_fault,
                                             mp.cellsize_fine, mp.cellsize_base)
                                             #,mp.fault_widths)

    else:
        print 'use mesh with a buffer zone with small cell sizes around faults'
        mesh = setup_mesh_with_exhumation_v2(mp.width, mp.fault_xs[0],
                                             mp.fault_widths[0],
                                             mp.fault_angles[0], elevation_top,
                                             z_surface + exhumed_thickness, z_surface,
                                             exhumation_steps,
                                             mp.z_fine, z_base, mp.cellsize,
                                             mp.cellsize_air, mp.cellsize_fault,
                                             mp.cellsize_fine, mp.cellsize_base,
                                             mp.fault_buffer_zone)


    ###############################################################
    # convert input params to escript variables
    # see here for more info:
    # https://answers.launchpad.net/escript-finley/+question/189076
    ###############################################################
    print 'converting input data to Escript variables'
    xyz = mesh.getX()

    #######################################
    # set up PDE solver
    #######################################
    print 'setting up PDE and boundary conditions'
    hf_pde = linearPDEs.LinearPDE(mesh)

    solver = ''
    if solver == 'GMRES':
        print 'using GMRES solver for heat transport PDE'
        hf_pde.getSolverOptions().setSolverMethod(es.SolverOptions.GMRES)
    elif solver is 'DIRECT':
        print 'using direct solver for heat transport PDE'
        hf_pde.getSolverOptions().setSolverMethod(
            es.SolverOptions.DIRECT)

    # find which nodes are on top & bottom boundaries
    surface = es.whereZero(xyz[1])
    top_bnd = es.whereZero(xyz[1] - es.sup(xyz[1]))
    bottom_bnd = es.whereZero(xyz[1] - es.inf(xyz[1]))

    # find which nodes are in the subsurface
    surface_level = exhumed_thickness
    subsurface = es.whereNonPositive(xyz[1] - surface_level)
    air = es.wherePositive(xyz[1] - surface_level)

    # set boundary conditions
    if mp.thermal_gradient is not None:
        bottom_temperature = (mp.air_temperature +
                              mp.thermal_gradient * mp.total_depth)
        specified_T_loc = es.wherePositive(top_bnd) + es.wherePositive(bottom_bnd)
        specified_T = es.wherePositive(top_bnd) * mp.air_temperature \
                      + es.wherePositive(bottom_bnd) * bottom_temperature
        specified_flux_loc = None
        specified_flux = None
    else:
        specified_T_loc = es.wherePositive(top_bnd)
        specified_T = es.wherePositive(top_bnd) * mp.air_temperature
        #specified_flux_loc = es.wherePositive(bottom_bnd)
        specified_flux_loc = es.Scalar(0, es.FunctionOnBoundary(mesh))
        specified_flux_loc.setTaggedValue("bottomleft", 1)
        specified_flux_loc.setTaggedValue("bottommid", 1)
        specified_flux_loc.setTaggedValue("bottomright", 1)
        #specified_flux = specified_flux_loc * mp.basal_heat_flux
        #specified_flux = mp.basal_heat_flux

        specified_flux = es.Scalar(0, es.FunctionOnBoundary(mesh))
        specified_flux.setTaggedValue("bottomleft", mp.basal_heat_flux)
        specified_flux.setTaggedValue("bottommid", mp.basal_heat_flux)
        specified_flux.setTaggedValue("bottomright", mp.basal_heat_flux)

    # populate porosity and K_solid values
    fault_x = calculate_fault_x(xyz[1], mp.fault_angles[0], mp.fault_xs[0])
    #K_solid = es.Scalar(0, es.FunctionOnBoundary(mesh))

    K_solid = xyz[0] * 0.0
    porosity = xyz[0] * 0.0

    for i, layer_bottom_i in enumerate(mp.layer_bottom):
        indl = es.whereNonPositive(xyz[0] - fault_x)
        indr = es.wherePositive(xyz[0] - fault_x)
        ind_layer_l = es.wherePositive(xyz[1] - layer_bottom_i[0])
        ind_layer_r = es.wherePositive(xyz[1] - layer_bottom_i[1])
        K_solid = K_solid * es.whereNonPositive(ind_layer_l * indl) \
                  + es.wherePositive(ind_layer_l * indl) * mp.K_solids[i]
        K_solid = K_solid * es.whereNonPositive(ind_layer_r * indr) \
                  + es.wherePositive(ind_layer_r * indr) * mp.K_solids[i]
        porosity = porosity * es.whereNonPositive(ind_layer_l * indl) \
                  + es.wherePositive(ind_layer_l * indl) * mp.porosities[i]
        porosity = porosity * es.whereNonPositive(ind_layer_r * indr) \
                  + es.wherePositive(ind_layer_r * indr) * mp.porosities[i]


    print 'matrix thermal conductivity: ', K_solid
    print 'porosity: ', porosity

    # calculate bulk material parameters
    K_b = K_solid ** (1 - porosity) * mp.K_water ** porosity
    rho_b = porosity * mp.rho_f + (1 - porosity) * mp.rho_s
    c_b = porosity * mp.c_f + (1 - porosity) * mp.c_s
    C1 = (mp.rho_f * mp.c_f) / (rho_b * c_b)

    # populate K, c and rho scalar fields
    K_var = subsurface * K_b + air * mp.K_air
    c_var = subsurface * c_b + air * mp.c_air
    rho_var = subsurface * rho_b + air * mp.rho_air

    # specify the location of faults
    #fault_zone = (subsurface * es.whereNegative(xyz[0] - fault_right)
    #              * es.wherePositive(xyz[1] - fault_bottom))

    fault_zones = []

    for fault_x, fault_angle, fault_width, fault_bottom \
            in zip(mp.fault_xs, mp.fault_angles,
                   mp.fault_widths, mp.fault_bottoms):

        depth = xyz[1]
        #x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

        fault_left = -depth * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + fault_x
        fault_right = fault_left + fault_width + 0.02
        fault_zone = ((subsurface * es.wherePositive(xyz[0] - fault_left))
                      * (subsurface * es.whereNegative(xyz[0] - fault_right))
                      * (subsurface * es.wherePositive(xyz[1] - fault_bottom)))

        print 'fault zone locs x ', fault_left, fault_right

        fault_zones.append(fault_zone)

        # TODO: find way to stop faults from crossing over,
        # perhaps have an order of faults, ie first main fault,
        # second fault cannot cross the first fault, etc....

    # add horizontal aquifers:
    aquifer_locs = []
    for aquifer_top, aquifer_bottom in zip(mp.aquifer_tops, mp.aquifer_bottoms):
        if aquifer_top is not None:
            aquifer_loc = (subsurface * es.whereNegative(xyz[1] - mp.total_depth - aquifer_top)
                           * es.wherePositive(xyz[1] - mp.total_depth - aquifer_bottom))

            aquifer_locs.append(aquifer_loc)

    # convert fault fluxes from m2/sec to m/sec
    # by dividing by fault zone width
    fault_fluxes_m_per_sec = []
    for fault_width, fault_flux in zip(mp.fault_widths, mp.fault_fluxes):

        fault_flux_i = [f / fault_width for f in fault_flux]
        fault_fluxes_m_per_sec.append(fault_flux_i)

    # model hydrothermal heating
    runtimes, T_steady, Ts, q_vectors, surface_levels, boiling_temps, exceed_boiling_temps = \
        model_hydrothermal_temperatures(
            mesh, hf_pde,
            fault_zones, mp.fault_angles, specified_T_loc, specified_flux_loc,
            mp.durations, fault_fluxes_m_per_sec,
            K_var, rho_var, c_var,
            mp.rho_f, mp.c_f,
            specified_T, specified_flux,
            mp.dt,
            top_bnd, bottom_bnd,
            mp.air_temperature, bottom_temperature,
            solve_as_steady_state=mp.steady_state,
            surface_level_init=surface_level, exhumation_rate=mp.exhumation_rate,
            target_depths=mp.target_zs,
            K_b=K_b, c_b=c_b, rho_b=rho_b,
            K_air=mp.K_air, c_air=mp.c_air, rho_air=mp.rho_air,
            vapour_correction=mp.vapour_correction)

    print 'T after thermal recovery ', Ts[-1]
    print 'done modeling'

    # convert modeled T field and vectors to arrays
    xyz_array, T0 = convert_to_array(Ts[0])
    T_list = [convert_to_array(T, no_coords=True) for T in Ts]
    T_array = np.array(T_list)
    T_init_array = convert_to_array(T_steady, no_coords=True)

    xyz_element_array, qh0 = convert_to_array(q_vectors[0][0])
    qh_list = [convert_to_array(q_vector[0], no_coords=True)
               for q_vector in q_vectors]
    qv_list = [convert_to_array(q_vector[1], no_coords=True)
               for q_vector in q_vectors]
    qh_array = np.array(qh_list)
    qv_array = np.array(qv_list)

    xyz_array_bt, b0 = convert_to_array(boiling_temps[-1])
    boiling_temp_list = [convert_to_array(maxT, no_coords=True) for maxT in boiling_temps]
    boiling_temp_array = np.array(boiling_temp_list)

    if np.max(xyz_array_bt - xyz_array_bt) > 0:
        print 'warning, node coords for boiling and T parameter are not the same'

    xyz_array_exc, bte_last = convert_to_array(exceed_boiling_temps[-1])
    bt_list = [convert_to_array(bt, no_coords=True) for bt in exceed_boiling_temps]
    exceed_boiling_temp_array = np.array(bt_list)

    ##############################################################
    # calculate temperature at depth slices (surface or otherwise)
    ##############################################################
    xzs = []
    Tzs = []
    nt, a = T_array.shape

    z_tolerance = 0.01

    for target_z in mp.target_zs:

        ind = np.abs(xyz_array[:, 1] - target_z) < z_tolerance
        xz = xyz_array[:, 0][ind]

        # new array for temperatures at depth
        Tzs_array = np.zeros((nt, len(xz)))

        # sort nodes in order of increasing x
        a = np.argsort(xz)
        xz = xz[a]

        # find temperature at target depth for each timestep:
        for i in range(nt):
            Tz = T_array[i, :][ind]
            Tzs_array[i, :] = Tz[a]

        xzs.append(xz)
        Tzs.append(Tzs_array)

    ##########################################
    # calculate helium ages at surface outcrop
    ##########################################
    #calculate_he_ages = True

    if mp.calculate_he_ages is False:
        Ahe_ages_all = None
        xs_Ahe_all = None

    else:

        # convert escript vars to arrays

        # calculate U-Th/He age for entire model domain.
        # assume no T change after K-Ar determined crystallization age
        # of 16.5 Ma
        ind_surface = np.where(xyz_array[:, 1] == 0)[0]
        nx = len(ind_surface)

        start_age = mp.t0 - runtimes[-1]

        t_prov = np.linspace(0, mp.t0, 31)
        T_prov = np.linspace(mp.T0, mp.T_surface, 31)

        print 'calculating helium ages'
        xs_Ahe_all = []
        Ahe_ages_all = []
        T_surf_mod_all = []

        for target_depth in mp.target_zs:
            #target_depth = 0
            nt = len(Ts)

            #ind_surface1 = xyz_array[:, 1] == target_depth
            z_tolerance = 0.01
            ind_surface = np.where(np.abs(xyz_array[:, 1] - target_depth) < z_tolerance)[0]
            nx = len(ind_surface)

            he_ages_surface = np.zeros((nt, nx))
            T_surf_mod = np.zeros((nt, nx))

            for xii in range(nx):

                t_he = np.concatenate((t_prov[:], t_prov[-1] + runtimes))
                T_he = np.concatenate((T_prov[:], T_array[:, ind_surface[xii]]))

                nt_prov = len(t_prov)
                #T_he *= 2

                T_he += mp.Kelvin

                he_age_i = he.calculate_he_age_meesters_dunai_2002(
                    t_he, T_he,
                    mp.radius, mp.U238, mp.Th232,
                    alpha_ejection=mp.alpha_ejection,
                    stopping_distance=mp.stopping_distance,
                    method=mp.AHe_method,
                    n_eigenmodes=50)

                # copy only He ages after provenance:
                for i in xrange(nt):
                    he_ages_surface[:, xii] = he_age_i[nt_prov:]

                #My = 1e6 * 365.25 * 24 * 60 * 60
                #print zip(t_he/My, T_he - mp.Kelvin, he_age_i / My)

            # get surface locs and T
            xs = xyz_array[:, 0][ind_surface]
            sort_order = np.argsort(xs)
            xs = xs[sort_order]

            # sort the data in order of increasing x
            ind_surface2 = ind_surface[sort_order]

            for i in xrange(nt):
                T_surf_mod[i, :] = T_array[i][ind_surface2]
                he_ages_surface[i] = he_ages_surface[i][sort_order]

            Ahe_ages_all.append(he_ages_surface)
            xs_Ahe_all.append(xs)
            T_surf_mod_all.append(T_surf_mod)

        AHe_data = [Ahe_ages_all, xs_Ahe_all]

        print 'done calculating helium ages'

        My = 1e6 * 365 * 24 * 60 * 60.



        print 'AHe ages: '
        for i, Ahe_ages in enumerate(Ahe_ages_all):
            try:
                print 'layer %i, min = %0.2f, mean = %0.2f, max = %0.2f My' \
                  % (i, Ahe_ages.min() / My,
                     Ahe_ages.mean() / My,
                     Ahe_ages.max() / My)
            except:
                print 'layer %i, no modeled AHe data availabel'

    print 'surface T: ', T * surface

    output = [runtimes, xyz_array, surface_levels,
              T_init_array, T_array, boiling_temp_array,
              xyz_array_exc, exceed_boiling_temp_array,
              xyz_element_array,
              qh_array, qv_array,
              fault_fluxes_m_per_sec, mp.durations, xzs, Tzs,
              Ahe_ages_all, xs_Ahe_all, mp.target_zs]

    return output


if __name__ == "__main__":

    print '-' * 20
    print 'running a single model scenario:'

    # import model parameters file
    from model_parameters.model_parameters import ModelParams

    mp = ModelParams

    scriptdir = os.path.realpath(sys.path[0])

    # run a single model scenario
    output = model_run(mp)

    (runtimes, xyz_array, surface_levels,
     T_init_array, T_array, boiling_temp_array,
     xyz_array_exc, exceed_boiling_temp_array,
     xyz_element_array,
     qh_array, qv_array,
     fault_fluxes, durations, xzs, Tzs,
     Ahe_ages_all, xs_Ahe_all, target_depths) = output

    # crop output to only the output timesteps, to limit filesize
    output_steps = []
    for duration, N_output in zip(mp.durations, mp.N_outputs):
        nt = int(duration / mp.dt)

        output_steps_i = list(np.linspace(0, nt-1, N_output).astype(int))
        output_steps += output_steps_i

    # select data for output steps only
    output_steps = np.array(output_steps)

    Tzs_cropped = [Tzi[output_steps] for Tzi in Tzs]

    if mp.calculate_he_ages is True:
        AHe_ages_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_all]
    else:
        AHe_ages_cropped = None

    N_output_steps = len(output_steps)

    # find surface temperatures
    T_surface = []
    x_surface = []

    for j in range(N_output_steps):

        # add output T at surface
        surface_elev = surface_levels[output_steps[j]]

        if surface_elev in target_depths:
            surface_ind = np.where(target_depths == surface_elev)[0]
            T_surface_i = Tzs[surface_ind][j]
            x_coords_i = xzs[surface_ind]

        else:
            # interpolate AHe age from nearest surfaces
            diff = target_depths - surface_elev
            ind_low = np.where(diff < 0)[0][-1]
            ind_high = np.where(diff > 0)[0][0]

            fraction = np.abs(diff[ind_low]) / (surface_levels[ind_high] - surface_levels[ind_low])

            T_surface_i = ((1.0-fraction) * Tzs[ind_low][j] + fraction * Tzs[ind_high][j])

            x_coords_i = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]

        T_surface.append(T_surface_i)
        x_surface.append(x_coords_i)

    # find AHe ages at surface
    if mp.calculate_he_ages is True:

        # add surface AHe data to output
        AHe_ages_surface = []
        AHe_xcoords_surface = []

        for i in range(N_output_steps):
            surface_elev = surface_levels[i]

            if surface_elev in target_depths:
                surface_ind = np.where(target_depths == surface_elev)[0]
                ages_raw = AHe_ages_cropped[surface_ind][i]
                x_coords = xzs[surface_ind]

            else:
                # interpolate AHe age from nearest surfaces
                diff = target_depths - surface_elev
                ind_low = np.where(diff < 0)[0][-1]
                ind_high = np.where(diff > 0)[0][0]

                fraction = np.abs(diff[ind_low]) / (target_depths[ind_high] - target_depths[ind_low])

                ages_raw = ((1.0-fraction) * AHe_ages_cropped[ind_low][i] + fraction * AHe_ages_cropped[ind_high][i])

                x_coords = (1.0-fraction) * xzs[ind_low] + fraction * xzs[ind_high]
            # add surface AHe data to output
            AHe_ages_surface.append(ages_raw)
            AHe_xcoords_surface.append(x_coords)
    else:
        AHe_ages_surface = None
        AHe_xcoords_surface = None

    output_selected = \
        [runtimes, runtimes[output_steps], xyz_array, surface_levels,
         T_init_array,
         T_array[output_steps], boiling_temp_array[output_steps],
         xyz_array_exc, exceed_boiling_temp_array[output_steps],
         xyz_element_array,
         qh_array[output_steps], qv_array[output_steps],
         fault_fluxes, durations,
         xzs, Tzs_cropped, x_surface, T_surface,
         AHe_ages_cropped, xs_Ahe_all, target_depths,
         AHe_ages_surface, AHe_xcoords_surface]

    #output_folder = os.path.join(folder, 'model_output')
    output_folder = os.path.join(scriptdir, mp.output_folder)

    if os.path.exists(output_folder) is False:
        print 'creating directory %s' % output_folder
        os.mkdir(output_folder)

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month,
                              today.year)

    fn = 'results_q_%s_%i_yr_grad_%0.0f_%s.pck' \
         % (str(np.array(mp.fault_fluxes) * mp.year),
            int(np.sum(np.array(durations) / mp.year)),
            mp.thermal_gradient * 1000.0,
            today_str)

    fn = fn.replace(' ', '')

    fn_path = os.path.join(output_folder, fn)

    print 'saving model results as %s' % fn_path

    fout = open(fn_path, 'w')
    pickle.dump(output_selected, fout)
    fout.close()

    print 'done with model scenario'

print 'done'



