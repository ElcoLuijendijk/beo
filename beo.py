"""
simple 2D model of advective heat flow

Elco Luijendijk, McGill university, 2013
"""

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

    #mp.fault_xs[0], mp.fault_widths[0], mp.fault_angles[0], mp.air_height,
    #              z_surface, mp.z_fine, z_base, mp.cellsize,
    #              mp.cellsize_air, mp.cellsize_fault,
    #              mp.cellsize_fine, mp.cellsize_base


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
                                    solve_as_steady_state=True):

    """

    """
    
    day = 24.0 * 60.0 * 60.0
    year = 365.25 * day

    #C1 = (rho_f * c_f) / (rho_var * c_var)

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
                print 'step %i of %i' % (t, nt)
                print 'temperature: ', T
            # solve PDE for temperature
            T = hf_pde.getSolution()

            # update PDE coefficients
            if solve_as_steady_state is False:
                Y = rho_var * c_var * T
                hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                r=specified_temperature,
                                q=specified_T_loc)

            t_total += dt

            #if t in output_steps:

            # store output
            Ts.append(T)
            q_vectors.append(q_vector)

            #ti = output_steps.index(t)
            #print 'surface T: ', T * surface

            runtimes.append(t_total)

        print 'T after convective heating ', T

    return np.array(runtimes), T_steady, Ts, q_vectors


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
    mesh = setup_mesh(mp.width, mp.fault_xs[0], mp.fault_widths[0],
                      mp.fault_angles[0], mp.air_height,
                      z_surface, mp.z_fine, z_base, mp.cellsize,
                      mp.cellsize_air, mp.cellsize_fault,
                      mp.cellsize_fine, mp.cellsize_base)

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

    # find which nodes are on top & bottom boundaries
    surface = es.whereZero(xyz[1])
    top_bnd = es.whereZero(xyz[1] - mp.air_height)
    bottom_bnd = es.whereZero(xyz[1] - es.inf(xyz[1]))

    # find which nodes are in the subsurface
    subsurface = es.whereNonPositive(xyz[1])
    air = es.wherePositive(xyz[1])

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

        #pdb.set_trace()

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

    # model hydrothermal heating
    runtimes, T_steady, Ts, q_vectors = \
        model_hydrothermal_temperatures(
            mesh, hf_pde,
            fault_zones, mp.fault_angles, specified_T_loc, specified_flux_loc,
            mp.durations, mp.fault_fluxes,
            K_var, rho_var, c_var,
            mp.rho_f, mp.c_f,
            specified_T, specified_flux,
            mp.dt,
            solve_as_steady_state=mp.steady_state)

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

    ##############################################################
    # calculate temperature at depth slices (surface or otherwise)
    ##############################################################
    xzs = []
    Tzs = []
    nt, a = T_array.shape

    for target_z in mp.target_zs:

        ind = xyz_array[:, 1] == target_z
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
        AHe_data = None

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

            ind_surface1 = xyz_array[:, 1] == target_depth
            ind_surface = np.where(xyz_array[:, 1] == target_depth)[0]
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
                #pdb.set_trace()

            # get surface locs and T
            xs = xyz_array[:, 0][ind_surface1]
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
            print 'layer %i, min = %0.2f, mean = %0.2f, max = %0.2f My' \
                  % (i, Ahe_ages.min() / My,
                     Ahe_ages.mean() / My,
                     Ahe_ages.max() / My)


    print 'surface T: ', T * surface

    output = [runtimes, xyz_array, T_init_array, T_array, xyz_element_array,
              qh_array, qv_array,
              mp.fault_fluxes, mp.durations, xzs, Tzs, Ahe_ages_all, xs_Ahe_all]

    return output


if __name__ == "__main__":

    print '-' * 20
    print 'running a single model scenario:'

    # import model parameters file
    import model_parameters.model_parameters as mp

    scriptdir = os.path.realpath(sys.path[0])

    # run a single model scenario
    output = model_run(mp)

    runtimes, xyz_array, T_init_array, T_array, xyz_element_array, \
        qh_array, qv_array, \
        fault_fluxes, durations, xzs, Tzs, Ahe_ages_all, xs_Ahe_all = output

    # crop output to only the output timesteps, to limit filesize
    output_steps = []
    for duration, N_output in zip(mp.durations, mp.N_outputs):
        nt = int(duration / mp.dt)

        output_steps_i = list(np.linspace(0, nt-1, N_output).astype(int))
        output_steps += output_steps_i

    # select data for output steps only
    output_steps = np.array(output_steps)

    Tzs_cropped = [Tzi[output_steps] for Tzi in Tzs]
    AHe_ages_cropped = [AHe_i[output_steps] for AHe_i in Ahe_ages_all]
    output_selected = \
        [runtimes, runtimes[output_steps], xyz_array, T_init_array,
         T_array[output_steps], xyz_element_array,
         qh_array[output_steps], qv_array[output_steps],
         fault_fluxes, durations, xzs, Tzs_cropped,
         AHe_ages_cropped, xs_Ahe_all]

    #output_folder = os.path.join(folder, 'model_output')
    output_folder = os.path.join(scriptdir, mp.output_folder)

    if os.path.exists(output_folder) is False:
        print 'creating directory %s' % output_folder
        os.mkdir(output_folder)

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month,
                              today.year)

    fn = 'results_q_%s_%i_yr_grad_%0.0f_%s.pck' \
         % (str(np.array(fault_fluxes) * mp.year),
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



