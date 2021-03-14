"""
2D model of advective and conductive heat flow and thermochronology in hydrothermal systems

Elco Luijendijk, McGill university & Goettingen University, 2013-2018
"""

__author__ = 'Elco Luijendijk'

# load modules

import time
import itertools
import pdb

import numpy as np
import pandas as pd
import scipy.interpolate

# escript/Finley modules:
import esys.escript as es
import esys.finley as fl
import esys.pycad as pc
import esys.pycad.gmsh as gmsh
import esys.escript.linearPDEs as linearPDEs

# import esys.weipa

# helium diffusion algorithm by Meesters and Dunai (2003)
import lib.helium_diffusion_models as he


def interpolate_data(xyz, data, target_x, y_steps=300):

    """

    """

    xi = np.array([xyz[:, 0].min(), target_x, xyz[:, 0].max()])
    yi = np.linspace(xyz[:, 1].min(), xyz[:, 1].max(), y_steps)
    nx = len(xi)
    ny = len(yi)

    xg, yg = np.meshgrid(xi, yi)
    xgf, ygf = xg.flatten(), yg.flatten()

    # interpolate u to grid
    zgf = scipy.interpolate.griddata(xyz, data, np.vstack((xgf, ygf)).T,
                                     method='linear')

    # make a 2d grid again
    zg = np.resize(zgf, (ny, nx))

    return yg[:, 1], zg[:, 1]


def calculate_vapour_pressure(T, c1=8e-8, c2=-7e-5, c3=0.028, c4=-3.1597):

    """
    Calculate the vapour pressure curve and check whether there should be
    vapour or not in the model domain

    based on a 3rd order polyonmial fit to vapour curve data by the NIST
    publication "Thermophysical Properties of Fluid Systems"
    found at http://webbook.nist.gov/chemistry/fluid/
    """

    log_Pv = c1 * T**3 + c2 * T**2 + c3 * T + c4

    Pv = 10**log_Pv * 1e6

    return Pv


def calculate_boiling_temp(P, c1=3.866, c2=25.151, c3=103.28, c4=179.99):
    """
    Find the maximum temperature for a given pressure at which there is one
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


def convert_to_array(u):

    """
    Return escript variable u as a numpy array
    """


    return np.array(u.toListOfTuples())



def convert_to_array_with_coords(u):

    """
    Return the x,y coordinates and the value of escript variable u as
    a numpy array
    """

    coords = u.getFunctionSpace().getX()
    x, y = coords[0], coords[1]

    xy = np.array([x.toListOfTuples(), y.toListOfTuples()]).T

    u_array = np.array(u.toListOfTuples())

    assert len(u_array.shape) == 1

    return xy, u_array


def Magnus_eq(T):
    """
    Calculate the saturation vapour pressure in air using the improved Magnus equation

    See Alduchov and Eskridge (1996) Journal of applied Meteorology 35 (4)

    """

    try:
        P = 0.61094 * np.exp(17.625 * T / (T + 243.04))
    except:
        P = 0.61094 * es.exp(17.625 * T / (T + 243.04))

    return P * 1000.0


def calculate_heat_transfer_coeff_sensible_hf(rho, c, ra, dz):
    """
    Calculate the heat transfer coefficient for sensible heat flux at the surface
    """

    Ksi = rho * c / ra * dz

    return Ksi


def calculate_heat_transfer_coeff_latent_hf(Ts, Ta, rho, dz, ra,
                                            L=2.264e6, Pa=1.0e5, RH_air=1.0):
    """
    Calculate the heat transfer coefficient for latent heat flux at the land surface
    """

    esa = Magnus_eq(Ta)
    qa_sat = 0.622 * esa / Pa
    qa = qa_sat * RH_air

    ess = Magnus_eq(Ts)
    qs = 0.622 * ess / Pa

    Kl = rho * L * dz / ra * (qs - qa) / (Ts - Ta)

    return Kl


def calculate_surface_heat_transfer_coeff(rho, c, ra, dz, Ts, Ta):
    """
    Calculate the heat transfer coefficient of air at the land surface

    """

    Ks = calculate_heat_transfer_coeff_sensible_hf(rho, c, ra, dz)
    Kl = calculate_heat_transfer_coeff_latent_hf(Ts, Ta, rho, dz, ra)
    Kt = Ks + Kl

    return Kt


def calculate_fault_x(z_flt, fault_angle, x_flt_surface):
    """
    Calculate the x coordinate of a fault with a given angle and intersection with the surface
    """

    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    return x_flt


def setup_mesh_with_exhumation_new(width, x_flt_surface, fault_width, fault_angle, fault_bottom,
                                   z_air,
                                   z_fine, z_base, cellsize,
                                   cellsize_air, cellsize_surface, cellsize_fault_surface, cellsize_fault,
                                   cellsize_fine, cellsize_base, target_zs,
                                   discretize_borehole, borehole_xs, borehole_depths, borehole_cellsize,
                                   x_left=0):

    """
    Create a mesh for the model of the beowawe hydrothermal system, including exhumation

    new version, fault at surface is at x=0 by definition, width denotes the width of the model domain on
    either side of the fault
    """

    ###############################
    # use gmsh to construct domain
    ##############################
    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain

    zs_surface = target_zs[::-1]

    z_flt = np.concatenate((zs_surface, np.array([z_fine, fault_bottom])))

    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    print 'x, z coords of fault locations in mesh:'
    for x, z in zip(x_flt, z_flt):
        print x, z

    x_left_bnd = np.min(x_flt) - width
    x_right_bnd = np.max(x_flt) + width

    check_x_bnds = False
    if check_x_bnds is True:
        if np.min(x_flt) <= x_left:
            print 'warning the left hand side of the fault is at %0.2f, ' \
                  'which is outside the model domain boundary at x=%0.2f' % (np.min(x_flt), x_left)
            x_left = np.min(x_flt) - fault_width * 2
            print 'relocating left hand boundary to %0.2f ' % x_left

        if np.max(x_flt) >= width:
            print 'warning the right hand side of the fault is at %0.2f, ' \
                  'which is outside the model domain boundary at x=%0.2f' % (np.max(x_flt), width)
            width = np.max(x_flt) + fault_width * 2
            print 'relocating right hand boundary to %0.2f ' % width

    xys = [[[x_left_bnd, z_air], [x_flt[0] -fault_width / 2.0, z_air],
            [x_flt[0] + fault_width / 2.0, z_air], [x_right_bnd, z_air]]]

    for xf, zf in zip(x_flt, z_flt):
        xys.append([[x_left_bnd, zf],
                    [xf - fault_width / 2.0, zf],
                    [xf + fault_width / 2.0 + 0.02, zf],
                    [x_right_bnd, zf]])
    xys.append([[x_left_bnd, z_base], [x_right_bnd, z_base]])

    print 'corner points in mesh: '
    for xysi in xys:
        print xysi

    points = []
    for xyi in xys:
        points_local = [pc.Point(x, z) for x, z in xyi]
        points.append(points_local)

    # cellsize in air layer:
    for p in points[0]:
        p.setLocalScale(cellsize_air / cellsize)

    # cellsize at land surface
    for point in points[1:-3]:
        for p in point:
            p.setLocalScale(cellsize_surface / cellsize)

    # cellsize in layer close to surface
    for point in points[-3]:
        point.setLocalScale(cellsize_fine / cellsize)

    # cellsize at fault surface
    for point_i in points[:-3]:
        point_i[1].setLocalScale(cellsize_fault_surface / cellsize)
        point_i[2].setLocalScale(cellsize_fault_surface / cellsize)

    # cellsize in fault:
    for point in points[-3:-1]:
        point[1].setLocalScale(cellsize_fault / cellsize)
        point[2].setLocalScale(cellsize_fault / cellsize)

    # cellsize at lower left and right bnds
    points[-2][0].setLocalScale(cellsize_base / cellsize)
    points[-2][-1].setLocalScale(cellsize_base / cellsize)

    points[-1][0].setLocalScale(cellsize_base / cellsize)
    points[-1][-1].setLocalScale(cellsize_base / cellsize)

    # horizontal lines:
    #hlines = [[pc.Line(points[0][0], points[0][1])]]
    hlines = []
    for point in points[0:-1]:
        hline_local = [pc.Line(point[0], point[1]),
                       pc.Line(point[1], point[2]),
                       pc.Line(point[2], point[3])]
        hlines.append(hline_local)

    hlines.append([pc.Line(points[-1][0], points[-1][1])])

    # vertical lines:
    #vlines = [[pc.Line(points[0][0], points[1][0]), pc.Line(points[0][1], points[1][3])]]
    vlines = []
    for point, point_below in zip(points[0:-2], points[1:]):
        vline_local = [pc.Line(point[0], point_below[0]),
                       pc.Line(point[1], point_below[1]),
                       pc.Line(point[2], point_below[2]),
                       pc.Line(point[3], point_below[3])]
        vlines.append(vline_local)

    # bottom lines
    vlines.append([pc.Line(points[-2][0], points[-1][0]), pc.Line(points[-2][3], points[-1][1])])

    #curves = [pc.CurveLoop(hlines[0][0], vlines[0][1],
    #                       -hlines[1][2], -hlines[1][1], -hlines[1][0],
    #                       -vlines[0][0])]
    curves = []
    for hline, hline_below, vline in zip(hlines[0:-2], hlines[1:-1], vlines[:]):
        curve_local_left = pc.CurveLoop(hline[0], vline[1], -hline_below[0], -vline[0])
        curve_local_fault = pc.CurveLoop(hline[1], vline[2], -hline_below[1], -vline[1])
        curve_local_right = pc.CurveLoop(hline[2], vline[3], -hline_below[2], -vline[2])

        curves += [curve_local_left, curve_local_fault, curve_local_right]

    curve_bottom = pc.CurveLoop(hlines[-2][0], hlines[-2][1], hlines[-2][2], vlines[-1][1], -hlines[-1][0], -vlines[-1][0])

    curves.append(curve_bottom)

    surfaces = [pc.PlaneSurface(curve) for curve in curves]

    d = gmsh.Design(dim=2, element_size=cellsize)

    d.setMeshFileName('beowawe_mesh')

    # add additional lines for borehole locations
    if discretize_borehole is True:
        z_top = max(target_zs)
        for x, depth in zip(borehole_xs, borehole_depths):
            p_top = pc.Point(x, z_top)
            p_bottom = pc.Point(x, -depth)
            bh_line = pc.Line(p_top, p_bottom)

            bh_line.setLocalScale(borehole_cellsize)
            d.addItems(bh_line)

    d.addItems(*surfaces)

    #labels = ["bottomleft", "bottommid", "bottomright"]
    #for hlinei, label in zip(hlines[-1], labels):
    #    ps = pc.PropertySet(label, hlinei)
    #    d.addItems(ps)

    #labels = ["bottomleft", "bottommid", "bottomright"]
    #for hlinei, label in zip(hlines[-1], labels):
    #    ps = pc.PropertySet(label, hlinei)
    #    d.addItems(ps)

    labels = ["bottom"]
    for hlinei, label in zip(hlines[-1], labels):
        ps = pc.PropertySet(label, hlinei)
        d.addItems(ps)

    mesh = fl.MakeDomain(d, optimizeLabeling=True)

    mesh.write('mesh.stl')

    return mesh, x_flt[:-2], z_flt[:-2], labels


def setup_mesh_2faults(width, x_flt_surface, fault_width, fault_angle, z_air,
                       z_surface, z_fine, z_base, cellsize,
                       cellsize_air, cellsize_fault, cellsize_fine, cellsize_base):

    """
    Create a mesh for a hydrothermal model containing 2 faults
    and one shallow aquifer
    """

    ###############################
    # use GMSH to construct domain
    ##############################
    # calculate fault positions
    # TODO: enable multiple faults, right now only one fault in model domain
    z_flt = np.array([z_surface, z_fine, z_base])
    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    xys = [[0, z_air], [width, z_air],
           [x_flt[0], z_surface], [x_flt[0] + fault_width, z_surface], [width, z_surface],
           [0, z_fine], [x_flt[1], z_fine], [x_flt[1] + fault_width, z_fine], [width, z_fine],
           [0, z_base], [x_flt[2], z_base], [x_flt[2] + fault_width, z_base], [width, z_base]]

    points = [pc.Point(x, z) for x, z in xys]

    # coarse cellsize at lower right
    points[0].setLocalScale(cellsize_air / cellsize)
    points[1].setLocalScale(cellsize_air / cellsize)

    points[2].setLocalScale(cellsize_fault / cellsize)
    points[3].setLocalScale(cellsize_fault / cellsize)
    points[6].setLocalScale(cellsize_fault / cellsize)
    points[7].setLocalScale(cellsize_fault / cellsize)
    points[10].setLocalScale(cellsize_fault / cellsize)
    points[11].setLocalScale(cellsize_fault / cellsize)

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


def model_hydrothermal_temperatures(mesh, hf_pde,
                                    fault_zones, fault_segments_all, fault_angles,
                                    specified_T_loc, specified_flux_loc,
                                    durations, fault_fluxes,
                                    aquifer_locs, aquifer_fluxes_m_per_sec, aquifer_angles,
                                    K_var, rho_var, c_var,
                                    rho_f, c_f,
                                    specified_temperature, specified_flux,
                                    dt,
                                    top_bnd, bottom_bnd, air_height,
                                    air_temperature, bottom_temperature,
                                    solve_as_steady_state=True,
                                    surface_level_init=0, exhumation_rate=0,
                                    target_depths=None,
                                    K_b=None, c_b=None, rho_b=None,
                                    K_air=None, c_air=None, rho_air=None,
                                    vapour_correction=True,
                                    variable_K_air=False, ra=80.0, reference_z_ra=1.8,
                                    steady_state_iterations=10,
                                    store_results_interval=1,
                                    track_exact_surface_elev=False,
                                    adv_pde=None,
                                    max_screen_output=500,
                                    surface_buffer=0.1,
                                    discretization='implicit',
                                    debug=False):

    """
    Full single model run of the heat conduction and advection model
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

    kronecker_delta = es.kronecker(mesh)

    xyz = mesh.getX()

    exceed_max_liquid_T_old = None

    ######################################
    # model steady-state temperature field
    ######################################
    # set PDE coefficients, steady-state heat flow equation
    A = K_var * kronecker_delta
    C = 0
    D = 0

    store_results_count = 0

    if specified_flux is not None:
        print 'solving with specified heat flux bnd'
        specified_heat_flux = specified_flux * specified_flux_loc
        #specified_heat_flux = specified_flux

        hf_pde.setValue(A=A, D=D,
                        r=specified_temperature,
                        q=specified_T_loc,
                        y=specified_heat_flux)
    else:
        print 'only fixed T bnd'
        hf_pde.setValue(A=A, D=D,
                        r=specified_temperature,
                        q=specified_T_loc)

    print 'finding solution for steady state HF'
    T_steady = hf_pde.getSolution()
    T = T_steady

    print 'modeled steady state temperatures ', T_steady

    print 'starting transient heat flow calculations'
    runtimes = [0.0]
    t_total = 0
    q_vectors = [es.Vector((0, 0), es.Function(mesh))]
    Ts = [T_steady]
    surface_level = surface_level_init
    surface_levels = [surface_level]

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

        exceed_boiling_temp = \
            subsurface * es.wherePositive(T - boiling_temp)

        boiling_temps = [boiling_temp]
        exceed_boiling_temps = [exceed_boiling_temp]

        if es.sup(vapour) >= 1:
            print 'warning, vapour present at initial steady-state P-T conditions'
    else:
        boiling_temps = None
        exceed_boiling_temps = None

    try:
        assert len(fault_fluxes) == len(durations)
    except AssertionError:
        msg = 'error, the length of the fault_fluxes and durations arrays or list do not match'
        msg += '\ncheck you model parameters file'
        raise ValueError(msg)

    n_time_periods = len(durations)

    xyz = mesh.getX()

    for time_period, fault_flux, duration in zip(itertools.count(), fault_fluxes, durations):

        print 'time period %i of %i' % (time_period + 1, n_time_periods)
        print 'duration of %i' % duration
        print 'fluxes in faults: ', fault_flux

        # set up advective flow field in faults and aquifers
        q_vector = es.Vector((0, 0), es.Function(mesh))
        xyze = q_vector.getFunctionSpace().getX()

        # add flux in faults
        for h, fault_zone, fault_angle, q_fault_zone in \
                zip(itertools.count(), fault_zones, fault_angles, fault_flux):

            if fault_segments_all is not None:
                for n_segment, fault_zone_segment, q_fault_zone_segment in zip(itertools.count(),
                                                                               fault_segments_all[h], q_fault_zone):
                    qh_fault_zone = - q_fault_zone_segment * np.cos(np.deg2rad(fault_angle))
                    qv_fault_zone = q_fault_zone_segment * np.sin(np.deg2rad(fault_angle))
                    q_vector[0] += fault_zone_segment * qh_fault_zone
                    q_vector[1] += fault_zone_segment * qv_fault_zone

                    print 'adding fault flux of %0.2e to fault segment %i of fault %i' \
                          % (q_fault_zone_segment, n_segment, h)
            else:
                # add heat advection in fault zone
                qh_fault_zone = - q_fault_zone * np.cos(np.deg2rad(fault_angle))
                qv_fault_zone = q_fault_zone * np.sin(np.deg2rad(fault_angle))
                q_vector[0] += fault_zone * qh_fault_zone
                q_vector[1] += fault_zone * qv_fault_zone

        # add fluxes in aquifers
        if aquifer_locs != []:
            for k, aquifer_loc, aquifer_flux, aquifer_angle in zip(itertools.count(), aquifer_locs,
                                                    aquifer_fluxes_m_per_sec[time_period], aquifer_angles):
                print 'adding aquifer flux %0.2e to aquifer %i' % (aquifer_flux, k)
                aquifer_flux_x = aquifer_flux * np.cos(np.deg2rad(aquifer_angle))
                aquifer_flux_y = aquifer_flux * np.sin(np.deg2rad(aquifer_angle))

                q_vector[0] += aquifer_loc * aquifer_flux_x
                q_vector[1] += aquifer_loc * aquifer_flux_y

        print 'modeled advective flux:'
        print '\tqh = ', q_vector[0]
        print '\tqz = ', q_vector[1]

        # make sure only flow in subsurface
        depth_ele = surface_level - xyze[1]
        subsurface_ele = es.wherePositive(depth_ele - surface_buffer)

        q_vector[0] = q_vector[0] * subsurface_ele
        q_vector[1] = q_vector[1] * subsurface_ele

        ###############################################
        # model transient response to fluid convection
        ###############################################
        if solve_as_steady_state is False:
            # set PDE coefficients, transient heat flow equation
            A = dt * K_var * kronecker_delta
            C = dt * rho_f * c_f * q_vector
            D = rho_var * c_var
            Y = rho_var * c_var * T

        else:
            print 'solving steady-state, with advective flux'

            # set PDE coefficients, steady-state flow equation
            A = K_var * kronecker_delta
            C = rho_f * c_f * q_vector
            D = 0
            Y = 0

        # update bnd cond if spec flux bnd
        if specified_flux is not None:

            print '\tadding specified flux bnd'
            if solve_as_steady_state is False:
                specified_heat_flux = specified_flux * specified_flux_loc * dt
            else:
                specified_heat_flux = specified_flux * specified_flux_loc

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

        nt = int(duration / dt)

        # set screen output interval
        screen_output_interval = int(np.ceil(nt / max_screen_output))
        if screen_output_interval < 1:
            screen_output_interval = 1
        print 'screen output once every %i timesteps' % screen_output_interval

        # calculate grid Peclet number
        Pe = rho_f * c_f * q_vector * mesh.getSize() / K_var
        print 'max. grid peclet number = ', es.Lsup(Pe)

        #############################################
        # iterate heat flow eq.s
        #############################################
        print 'x' * 10
        print 'starting iterations, total timesteps = %i' % nt
        start = time.time()

        if solve_as_steady_state is True:

            if variable_K_air is True or vapour_correction is True:
                nt = steady_state_iterations
                print 'running steady-state model with %i iterations for variable K air or vapour correction' % nt
            else:
                nt = 1
                print 'running steady-state model without iterations'

        update_PDE = False

        # solve transient heat flux
        for t in range(nt):

            if vapour_correction is True:
                vapour_pressure = calculate_vapour_pressure(T)
                vapour = subsurface * es.whereNegative(P - vapour_pressure + P_buffer)

                xmin_vapour = es.inf(vapour * xyz[0] + es.whereZero(vapour) * 999999.9)
                xmax_vapour = es.sup(vapour * xyz[0])
                ymin_vapour = -es.sup(-(vapour * xyz[1]))
                ymax_vapour = es.sup(vapour * xyz[1])

            if exhumation_rate > 0:
                surface_level = surface_level_init - t_total / year * exhumation_rate

                try:
                    surface_level_mesh_id = np.where(target_depths <= surface_level)[0][-1]
                    surface_level_mesh = target_depths[surface_level_mesh_id]
                except:
                    surface_level_mesh = surface_level
                    # print '\twarning could not find land surface nodes'

            # calculate effective thermal conductivity air layer based on latent and sensible heat flux
            # eqs.
            if variable_K_air is True:

                xysa, Tsa = convert_to_array_with_coords(surface * T)
                sc = convert_to_array(surface)

                # find x coords surface nodes
                ind_s = sc == 1
                xa1 = xysa[:, 0][ind_s]

                # get temperature surface nodes
                Tsa1 = Tsa[ind_s]

                # sort in ascending order
                a = np.argsort(xa1)
                Tsa2 = Tsa1[a]

                int_table = Tsa2

                numslices = len(int_table) - 1
                minval = es.inf(xyz[0])
                maxval = es.sup(xyz[0])
                step = es.sup(maxval - minval) / numslices
                toobig = int_table.max() * 1.5

                # interpolate surface temperature
                surface_T_int = es.interpolateTable(int_table, xyz[0], minval, step, toobig)

                # calculate heat transfer coefficient air:
                K_air = calculate_surface_heat_transfer_coeff(rho_air, c_air, ra,
                                                              reference_z_ra, surface_T_int,
                                                              air_temperature)

                if debug is True:
                    print 'save as csv file?'

                    if raw_input() == 'y':

                        es.saveDataCSV('debug/interpolated_surface_T.csv', x=xyz[0], y=xyz[1],
                                       surface_T=surface_T_int)

            # Exhumation: calculate new surface level
            if exhumation_rate != 0 and solve_as_steady_state is False:

                if track_exact_surface_elev is True:

                    depth = surface_level - xyz[1]
                    depth_ele = surface_level - xyze[1]

                else:
                    depth = surface_level_mesh - xyz[1]
                    depth_ele = surface_level_mesh - xyze[1]

                subsurface = es.wherePositive(depth)
                subsurface_ele = es.wherePositive(depth_ele)
                subsurface_ele_buffer = es.wherePositive(depth_ele - surface_buffer)

                air = es.whereNegative(depth)
                #air_ele = es.whereNegative(depth_ele)
                surface = es.whereZero(depth)
                #surface_ele = es.whereZero(depth_ele)

                q_vector_old = q_vector.copy()
                q_vector[0] = q_vector[0] * subsurface_ele_buffer
                q_vector[1] = q_vector[1] * subsurface_ele_buffer

            # populate K, c and rho scalar fields in case of exhumation or variable K air
            if variable_K_air is True or exhumation_rate != 0:

                K_var = subsurface * K_b + air * K_air
                c_var = subsurface * c_b + air * c_air
                rho_var = subsurface * rho_b + air * rho_air

                # update bnd cond if spec flux bnd
                update_PDE = True

                if specified_flux is not None:
                    if solve_as_steady_state is False:
                        specified_heat_flux = specified_flux * specified_flux_loc * dt
                    else:
                        specified_heat_flux = specified_flux * specified_flux_loc

            # recalculate fluid pressure
            if vapour_correction is True:
                depth = -(xyz[1] - surface_level)
                P_init = fluid_density * g * depth + atmospheric_pressure
                P = P_init * es.wherePositive(depth) \
                    + atmospheric_pressure * es.whereNonPositive(depth)

                logP = es.log10(P / 1.0e6)
                boiling_temp = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4
            else:
                boiling_temp = 1e6

            # recalculate vapour pressure and max liquid temperature
            if vapour_correction is True:
                #vapour_pressure = calculate_vapour_pressure(T)
                #boiling_temp = calculate_boiling_temp(P)
                #exceed_boiling_temp = subsurface * es.wherePositive(T - boiling_temp)

                if exceed_max_liquid_T_old is None:
                    exceed_boiling_temp = \
                        subsurface * es.wherePositive(T - boiling_temp)
                else:
                    exceed_boiling_temp = \
                        subsurface * es.whereZero(exceed_max_liquid_T_old) \
                        * es.wherePositive(T - boiling_temp) #\
                        #+ subsurface * exceed_max_liquid_T_old

                exceed_max_liquid_T_old = exceed_boiling_temp
            else:
                exceed_boiling_temp = 0

            if vapour_correction is True:

                height = xyz[1] - surface_level
                top_bnd = es.whereNonNegative(height - air_height)

                if bottom_temperature is not None:
                    specified_T_loc = es.wherePositive(top_bnd) \
                                      + es.wherePositive(bottom_bnd) \
                                      + exceed_boiling_temp
                    specified_temperature = \
                        es.wherePositive(top_bnd) * air_temperature \
                        + es.wherePositive(bottom_bnd) * bottom_temperature \
                        + exceed_boiling_temp * boiling_temp
                else:
                    specified_T_loc = es.wherePositive(top_bnd) \
                                      + exceed_boiling_temp
                    specified_temperature = \
                        es.wherePositive(top_bnd) * air_temperature \
                        + exceed_boiling_temp * boiling_temp

                update_PDE = True

            elif exhumation_rate > 0.0:

                height = xyz[1] - surface_level
                top_bnd = es.whereNonNegative(height - air_height)

                if bottom_temperature is not None:
                    specified_T_loc = es.wherePositive(top_bnd) \
                                      + es.wherePositive(bottom_bnd)
                    specified_temperature = \
                        es.wherePositive(top_bnd) * air_temperature
                else:
                    specified_T_loc = es.wherePositive(top_bnd)
                    specified_temperature = \
                        es.wherePositive(top_bnd) * air_temperature

                update_PDE = True

            if update_PDE is True:

                if solve_as_steady_state is True:
                    # set PDE coefficients, steady-state flow equation
                    A = K_var * kronecker_delta
                    C = rho_f * c_f * q_vector
                    D = 0
                    Y = 0

                elif discretization == 'implicit':

                    # reset heatflow PDE coefficients
                    A = dt * K_var * kronecker_delta
                    C = dt * rho_f * c_f * q_vector
                    D = rho_var * c_var
                    Y = rho_var * c_var * T

                    if specified_flux is not None:
                        hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                        r=specified_temperature,
                                        q=specified_T_loc,
                                        y=specified_heat_flux)
                    else:
                        hf_pde.setValue(A=A, C=C, D=D, Y=Y,
                                        r=specified_temperature,
                                        q=specified_T_loc)

                elif discretization == 'explicit':

                    # set PDE coefficients, explicit form of heat flow eq..
                    # not implemented yet...
                    pass

            # solve PDE for temperature
            #print '\tsolving for T'
            solve_time_start = time.time()
            T = hf_pde.getSolution()
            solver_time = time.time() - solve_time_start

            update_PDE = True
            #print '\tnew T ', T

            # update PDE coefficients
            #if solve_as_steady_state is False:
            #    Y = rho_var * c_var * T
            #    update_PDE = True

            t_total += dt
            if solve_as_steady_state is True:
                t_total = 0.0

            # output to screen:
            if int(t / float(screen_output_interval)) == float(t / float(screen_output_interval)) \
                    or t == nt - 1:

                end = time.time()
                comptime = end - start

                start = end

                print 'step %i of %i' % (t + 1, nt)

                print '\truntime = %0.2e yrs' % (t_total / year)
                print '\tcomputational time for one timestep = %0.3f sec' \
                      % (comptime / screen_output_interval)
                print '\tcomputational time needed by solver = %0.3f sec' \
                      % solver_time
                print '\tactual surface level ', surface_level
                if exhumation_rate > 0:
                    print '\tclosest surface in mesh ', surface_level_mesh
                print '\ttemperature: ', T
                if es.sup(surface) > 0:
                    print '\tmax. temperature at land surface: ', \
                        es.sup(T * surface)
                else:
                    print '\tcould not find land surface nodes'

                if vapour_correction is True:
                    if es.sup(vapour) > 0:
                        print '\tvapour present in: ', es.integrate(vapour), ' m^2'
                        print '\t\tfrom x = ', xmin_vapour, ' to x = ', \
                            xmax_vapour
                        print '\t\tand from y = ', ymin_vapour, ' to y = ', \
                            ymax_vapour
                        # print '\tmax. liquid T at the surface = ', \
                        #    es.sup(boiling_temp * land_surface)
                    else:
                        print '\tno vapour present'

                if variable_K_air is True:
                    print '\tinterpolated surface T ', surface_T_int
                    print '\tcalculated K air ', K_air

            # store output
            store_results_count += 1
            if store_results_count >= store_results_interval:
                Ts.append(T)
                q_vectors.append(q_vector.copy())

                surface_levels.append(surface_level)

                if vapour_correction is True:
                    boiling_temps.append(boiling_temp)
                    exceed_boiling_temps.append(exceed_boiling_temp)

                runtimes.append(t_total)

                store_results_count = 0

        print 'final T after advective heating ', T

    return (np.array(runtimes), T_steady, Ts, q_vectors,
            np.array(surface_levels), boiling_temps, exceed_boiling_temps)


def model_run(mp):

    """
    setup mesh and parameters, and run a single model experiment

    :param mp:
    :return:
    """

    #
    year = 365.25 * 24 * 60 * 60
    Myr = year * 1e6

    ############################
    # construct rectangular mesh
    ############################
    print 'constructing mesh (note, this may take a while....)'

    z_surface = 0
    z_base = -mp.total_depth

    try:
        assert np.min(np.array(mp.target_zs)) > mp.z_fine
    except AssertionError:
        msg = 'error, the lowest value of the model parameter target_z should always be higher than z_fine'
        msg += ', otherwise the grid algorithm does not work.\n'
        msg += 'target_zs = %s\n' % str(mp.target_zs)
        msg += 'z_fine = %0.1f\n' % mp.z_fine
        #raise AssertionError(msg)
        mp.z_fine = np.min(np.array(mp.target_zs)) - mp.cellsize_fine

        msg += 'lowering z_fine to %0.1f' % mp.z_fine

        print msg

    if mp.add_exhumation is False:
        print 'construct static mesh without exhumation'

        # make sure to set exhumation rate to 0, will get errors in code otherwise...
        mp.exhumation_rate = 0.0

        elevation_top = z_surface + mp.air_height
        exhumation_steps = 0
        exhumed_thickness = 0

        mesh, x_flt, z_flt, bottom_labels = setup_mesh_with_exhumation_new(
                                                        mp.width, mp.fault_xs[0],
                                                        mp.fault_widths[0],
                                                        mp.fault_angles[0], mp.fault_bottoms[0] - 500.0,
                                                        elevation_top,
                                                        mp.z_fine, z_base, mp.cellsize,
                                                        mp.cellsize_air, mp.cellsize_surface,
                                                        mp.cellsize_fault_surface, mp.cellsize_fault,
                                                        mp.cellsize_fine, mp.cellsize_base, mp.target_zs,
                                                        mp.discretize_borehole, mp.borehole_xs, mp.borehole_depths,
                                                        mp.borehole_cellsize)

        #def setup_mesh_with_exhumation_new(width, x_flt_surface, fault_width, fault_angle, fault_bottom,
        #                                   z_air,
        #                                   z_fine, z_base, cellsize,
        #                                   cellsize_air, cellsize_surface, cellsize_fault,
        #                                   cellsize_fine, cellsize_base, target_zs,
        #                                   x_left=0):

        x_flt = np.ones(len(mp.target_zs)) * mp.fault_xs[0]
        z_flt = np.zeros(len(mp.target_zs))

        #exhumed_thickness = 0
        #elevation_top = z_surface + exhumed_thickness + mp.air_height

    else:
        print 'construct dynamic mesh with exhumation'
        exhumed_thickness = mp.exhumation_rate * (np.sum(np.array(mp.durations)) / mp.year)
        exhumation_steps = mp.exhumation_steps

        min_layer_thickness = mp.min_layer_thickness
        if exhumed_thickness / exhumation_steps < min_layer_thickness:
            print 'warning, exhumation levels would be smaller than %0.2f m' % min_layer_thickness
            exhumation_steps = int(np.ceil(exhumed_thickness) / min_layer_thickness)
            if exhumation_steps < 1:
                exhumation_steps = 1

            print 'reducing exhumation steps to %i' % exhumation_steps

        if exhumed_thickness != 0:
            # track AHe and temperature in each exhumed layer in the model domain:

            #test1 = np.concatenate((np.array(mp.target_zs), np.linspace(0, exhumed_thickness, exhumation_steps + 1)))
            #test2 = np.unique(test1)
            mp.target_zs = np.linspace(0, exhumed_thickness, exhumation_steps + 1)
            #mp.target_zs = test2

        elevation_top = z_surface + exhumed_thickness + mp.air_height

        mesh, x_flt, z_flt, bottom_labels = setup_mesh_with_exhumation_new(
                                                        mp.width, mp.fault_xs[0],
                                                        mp.fault_widths[0],
                                                        mp.fault_angles[0], mp.fault_bottoms[0] - 500.0,
                                                        elevation_top,
                                                        mp.z_fine, z_base, mp.cellsize,
                                                        mp.cellsize_air, mp.cellsize_surface,
                                                        mp.cellsize_fault_surface, mp.cellsize_fault,
                                                        mp.cellsize_fine, mp.cellsize_base, mp.target_zs,
                                                        mp.discretize_borehole,
                                                        mp.borehole_xs, mp.borehole_depths, mp.borehole_cellsize)


        #def setup_mesh_with_exhumation_new(width, x_flt_surface, fault_width, fault_angle, fault_bottom,
        #                                   z_air,
        #                                   z_fine, z_base, cellsize,
        #                                   cellsize_air, cellsize_surface, cellsize_fault,
        #                                   cellsize_fine, cellsize_base, target_zs,
        #                                   x_left=0):

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
    adv_pde = None

    solver = mp.solver
    if solver == 'GMRES':
        print 'using GMRES solver for heat transport PDE'
        hf_pde.getSolverOptions().setSolverMethod(es.SolverOptions.GMRES)
    elif solver is 'DIRECT':
        print 'using direct solver for heat transport PDE'
        hf_pde.getSolverOptions().setSolverMethod(
            es.SolverOptions.DIRECT)
    elif solver is 'ROWSUM_LUMPING':
        adv_pde = linearPDEs.LinearPDE(mesh)
        adv_pde.getSolverOptions().setSolverMethod(es.SolverOptions.ROWSUM_LUMPING)
        hf_pde.getSolverOptions().setSolverMethod(es.SolverOptions.GMRES)
    elif solver is 'PCG':
        #hf_pde.getSolverOptions().setSolverMethod(es.SolverOptions.ROWSUM_LUMPING)
        hf_pde.getSolverOptions().setSolverMethod(es.SolverOptions.PCG)
        hf_pde.getSolverOptions().setPreconditioner(es.SolverOptions.AMG)

    # find which nodes are on top & bottom boundaries
    top_bnd = es.whereZero(xyz[1] - es.sup(xyz[1]))
    bottom_bnd = es.whereZero(xyz[1] - es.inf(xyz[1]))

    # find which nodes are in the subsurface
    surface_level = exhumed_thickness
    surface = es.whereZero(xyz[1] - surface_level)
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
        bottom_temperature = None
        specified_T_loc = es.wherePositive(top_bnd)
        specified_T = es.wherePositive(top_bnd) * mp.air_temperature

        specified_flux_loc = es.Scalar(0, es.FunctionOnBoundary(mesh))
        specified_flux = es.Scalar(0, es.FunctionOnBoundary(mesh))

        for bottom_label in bottom_labels:
            specified_flux_loc.setTaggedValue(bottom_label, 1)
            specified_flux.setTaggedValue(bottom_label, mp.basal_heat_flux)

    # populate porosity and K_solid values
    fault_x = calculate_fault_x(xyz[1], mp.fault_angles[0], mp.fault_xs[0])
    # K_solid = es.Scalar(0, es.FunctionOnBoundary(mesh))

    K_solid = xyz[0] * 0.0
    porosity = xyz[0] * 0.0

    try:
        assert len(mp.layer_bottom) == len(mp.K_solids) == len(mp.porosities)
    except AssertionError:
        msg = 'error, number of layers, K_solid and porosities supplied in input file are not equal'
        msg += '\nlayers = %i, K_solid = %i, porosities = %i' % (len(mp.layer_bottom), len(mp.K_solids), len(mp.porosities))
        msg += '\n check your input file'
        raise ValueError(msg)

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
    if mp.variable_K_air is True:
        # base K_air on surface temperature of underlying nodes
        # assume a small difference in initial air and land temperature (otherwise the eq. doesnt work)
        surface_T = mp.air_temperature + 0.5
        K_air = calculate_surface_heat_transfer_coeff(mp.rho_air, mp.c_air, mp.ra,
                                                      mp.dz, surface_T, mp.air_temperature)
    else:
        K_air = mp.K_air
    K_var = subsurface * K_b + air * K_air
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

        fault_left = -depth * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 - fault_width / 2.0 + fault_x
        fault_right = fault_left + fault_width + 0.02
        fault_zone = ((subsurface * es.wherePositive(xyz[0] - fault_left))
                      * (subsurface * es.whereNegative(xyz[0] - fault_right))
                      * (subsurface * es.wherePositive(xyz[1] - fault_bottom)))

        print 'fault zone locs x ', fault_left, fault_right

        fault_zones.append(fault_zone)

        # TODO: find way to stop faults from crossing over,
        # perhaps have an order of faults, ie first main fault,
        # second fault cannot cross the first fault, etc....

    fault_segments_all = None
    if mp.fault_segments is not None:
        fault_segments_all = []
        for h, fault_x, fault_angle, fault_width, fault_bottom, fault_segments \
                in zip(itertools.count(), mp.fault_xs, mp.fault_angles,
                       mp.fault_widths, mp.fault_bottoms, mp.fault_segments):

            fault_segment_tops = fault_segments

            depth = xyz[1]

            fault_left = -depth * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 \
                         - fault_width / 2.0 + fault_x
            fault_right = fault_left + fault_width + 0.02

            fault_segments_i = []
            segment_bottoms = [fault_bottom] + list(fault_segment_tops[:-1])
            for segment_top, segment_bottom in zip(fault_segment_tops, segment_bottoms):
                fault_zone_segment = ((subsurface * es.wherePositive(xyz[0] - fault_left))
                          * (subsurface * es.whereNegative(xyz[0] - fault_right))
                          * (subsurface * es.wherePositive(xyz[1] - segment_bottom))
                          * (subsurface * es.whereNonPositive(xyz[1] - segment_top)))

                fault_segments_i.append(fault_zone_segment)

                print 'adding segment from z=%0.2f to %0.2f for fault %i ' \
                      % (segment_bottom, segment_top, h)
            fault_segments_all.append(fault_segments_i)

    # add horizontal aquifers:
    aquifer_locs = []
    for aquifer_top, aquifer_bottom, aquifer_left_bnd, aquifer_right_bnd, aquifer_angle \
            in zip(mp.aquifer_tops, mp.aquifer_bottoms,
                   mp.aquifer_left_bnds, mp.aquifer_right_bnds, mp.aquifer_angles):

        if aquifer_top is not None:

            try:
                assert aquifer_top > aquifer_bottom
            except AssertionError, msg:
                new_msg = 'error, the aquifer bottom is higher than the top, check input file'
                raise AssertionError(new_msg)

            aquifer_left_bnd_i = aquifer_left_bnd
            aquifer_right_bnd_i = aquifer_right_bnd

            if aquifer_left_bnd is None:
                aquifer_left_bnd_i = fault_right
            if aquifer_right_bnd is None:
                aquifer_right_bnd_i = fault_left

            aquifer_thickness_i = aquifer_top - aquifer_bottom
            aquifer_x = xyz[0] - aquifer_left_bnd_i
            aquifer_top_i = aquifer_top + aquifer_x * np.tan(np.deg2rad(aquifer_angle))
            aquifer_bottom_i = aquifer_top_i - aquifer_thickness_i

            aquifer_loc = (subsurface
                           * es.whereNegative(xyz[1] - aquifer_top_i)
                           * es.wherePositive(xyz[1] - aquifer_bottom_i)
                           * es.whereNegative(xyz[0] - aquifer_right_bnd_i)
                           * es.wherePositive(xyz[0] - aquifer_left_bnd_i))

            # remove overlap aquifer and fault zone
            for fault_segments_i in fault_segments_all:
                for fault_segment_ii in fault_segments_i:
                    aquifer_loc = aquifer_loc * es.whereZero(fault_segment_ii)

            aquifer_locs.append(aquifer_loc)

            aquifer_center_x = es.integrate(aquifer_loc * xyz[0]) / es.integrate(xyz[0])
            aquifer_center_y = es.integrate(aquifer_loc * xyz[1]) / es.integrate(xyz[1])

            print 'added aquifer centered on x= ', aquifer_center_x,' and z = ', aquifer_center_y

            try:
                assert(es.Lsup(aquifer_loc) >= 1.0)
            except AssertionError, msg:
                print 'error, something wrong with assigning the aquifer nodes. ' \
                      'check the aquifer params in the input file and the grid cell size'

    # find intersections fault and aquifers
    fault_int_locs = []
    fault_int_angles = []
    for i, fault_zone, fault_angle in zip(itertools.count(), fault_zones, mp.fault_angles):

        for aquifer_loc, aquifer_bottom, aquifer_top in zip(aquifer_locs, mp.aquifer_bottoms, mp.aquifer_tops):

            flow_cutoff_level = (aquifer_top + aquifer_bottom) / 2.0
            #flow_cutoff_level = aquifer_top

            fault_int_i = es.wherePositive(fault_zone) * es.wherePositive(xyz[1] - flow_cutoff_level)
            fault_int_locs.append(fault_int_i)
            fault_int_angles.append(fault_angle)

            # remove everything that crosses aquifer zone and above from fault zone
            if es.Lsup(fault_int_i) > 0:

                fault_zones[i] = fault_zones[i] * es.whereNonPositive(xyz[1] - flow_cutoff_level)

    #print 'fault fluxes: ', [fi * year for f in mp.fault_fluxes for fi in f]

    # convert fault fluxes from m2/sec to m/sec
    # by dividing by fault zone width
    fault_fluxes_m_per_sec = []
    for duration, fault_flux_timeslice in zip(mp.durations, mp.fault_fluxes):
        fts_timeslice = []
        for fault_width, fault_flux in zip(mp.fault_widths,  fault_flux_timeslice):
            #fts = []
            #for fault_flux_ts in fault_flux:
            fseg = []
            for fault_flux_segment in fault_flux:
                fault_flux_i = fault_flux_segment / fault_width
                fseg.append(fault_flux_i)
                #fts.append(fseg)
            fts_timeslice.append(fseg)
        fault_fluxes_m_per_sec.append(fts_timeslice)

    # convert aquifer fluxes from m2/sec to m/sec
    aquifer_fluxes_m_per_sec = []

    if aquifer_top is not None:

        for duration, aquifer_flux_timeslice in zip(mp.durations, mp.aquifer_fluxes):

            aquifer_flux_ii = []

            for aquifer_top, aquifer_bottom, aquifer_flux in zip(mp.aquifer_tops,
                                                                 mp.aquifer_bottoms,
                                                                 aquifer_flux_timeslice):
                if aquifer_top is not None:

                    aquifer_thickness = aquifer_top - aquifer_bottom

                    aquifer_flux_i = aquifer_flux / aquifer_thickness

                aquifer_flux_ii.append(aquifer_flux_i)
            aquifer_fluxes_m_per_sec.append(aquifer_flux_ii)

    store_results_interval = mp.dt_stored / mp.dt
    if float(store_results_interval) != int(store_results_interval) and store_results_interval != 1.0:
        msg = 'warning, dt_stored divided by dt is %s but should be an integer.' % str(store_results_interval)
        msg += '\nstoring results at each timestep instead.' % str(store_results_interval)
        print(msg)
        #raise ValueError(msg)

    # model hydrothermal heating
    runtimes, T_steady, Ts, q_vectors, surface_levels, boiling_temps, exceed_boiling_temps = \
        model_hydrothermal_temperatures(
            mesh, hf_pde,
            fault_zones, fault_segments_all, mp.fault_angles,
            specified_T_loc, specified_flux_loc,
            mp.durations, fault_fluxes_m_per_sec,
            aquifer_locs, aquifer_fluxes_m_per_sec, mp.aquifer_angles,
            K_var, rho_var, c_var,
            mp.rho_f, mp.c_f,
            specified_T, specified_flux,
            mp.dt,
            top_bnd, bottom_bnd, mp.air_height,
            mp.air_temperature, bottom_temperature,
            solve_as_steady_state=mp.steady_state,
            adv_pde=adv_pde,
            surface_level_init=surface_level,
            exhumation_rate=mp.exhumation_rate,
            target_depths=mp.target_zs,
            K_b=K_b, c_b=c_b, rho_b=rho_b,
            K_air=K_air, c_air=mp.c_air, rho_air=mp.rho_air,
            vapour_correction=mp.vapour_correction,
            variable_K_air=mp.variable_K_air, ra=mp.ra, reference_z_ra=mp.dz,
            steady_state_iterations=mp.n_iterations_steady_state,
            store_results_interval=store_results_interval)

    print 'T after model runs: ', Ts[-1]
    print 'number of saved temperature fields = %i' % (len(Ts))
    print 'done modeling'

    # convert modeled T field and vectors to arrays
    xyz_array, T0 = convert_to_array_with_coords(Ts[0])
    T_list = [convert_to_array(T) for T in Ts]
    T_array = np.array(T_list)

    T_list = None

    T_init_array = convert_to_array(T_steady)

    xyz_element_array, qh0 = convert_to_array_with_coords(q_vectors[0][0])
    qh_list = [convert_to_array(q_vector[0])
               for q_vector in q_vectors]
    qv_list = [convert_to_array(q_vector[1])
               for q_vector in q_vectors]
    qh_array = np.array(qh_list)
    qv_array = np.array(qv_list)

    qh_list = None
    qv_list = None

    if boiling_temps is not None:
        xyz_array_bt, b0 = convert_to_array_with_coords(boiling_temps[-1])
        boiling_temp_list = [convert_to_array(maxT) for maxT in boiling_temps]
        boiling_temp_array = np.array(boiling_temp_list)

        boiling_temp_list = None

        if np.max(xyz_array_bt - xyz_array_bt) > 0:
            print 'warning, node coords for boiling and T parameter are not the same'

        xyz_array_exc, bte_last = convert_to_array_with_coords(exceed_boiling_temps[-1])
        bt_list = [convert_to_array(bt) for bt in exceed_boiling_temps]
        exceed_boiling_temp_array = np.array(bt_list)
        bt_list = None

    else:
        boiling_temp_array = None
        xyz_array_exc = None
        exceed_boiling_temp_array = None

    ##############################################################
    # calculate temperature at depth slices (surface or otherwise)
    ##############################################################
    xzs = []
    Tzs = []
    Tzs_diff = []
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

        Tzs_diff_array = np.zeros_like(Tzs_array)
        # find temperature difference with initial steady-state temperatures
        for i in range(nt):
            Tzs_diff_array[i, :] = Tzs_array[i, :] - T_init_array[ind]

        Tzs_diff.append(Tzs_diff_array)

    #######################
    # calculate helium ages
    #######################

    if mp.calculate_he_ages is False:
        Ahe_ages_surface_all = None
        Ahe_ages_surface_corr_all = None
        AHe_ages_surface_samples_all = None
        AHe_ages_surface_samples_all_corr = None
        xs_AHe_surface_all = None
        z_borehole, he_ages_borehole, he_ages_borehole_corrected = None, None, None

    else:

        # find locations of AHe samples:
        if mp.model_AHe_samples is True:

            # read sample data
            dfhs = pd.read_csv(mp.AHe_data_file)

            if mp.profile_number in dfhs['profile'].values:
                print 'selecting profile %i'
                locs = dfhs['profile'] == mp.profile_number
                dfhs = dfhs.loc[locs]
            else:
                print 'not selecting profile, keeping all AHe data'

            AHe_sample_names = dfhs['sample']
            AHe_sample_distances = dfhs['distance'].values
            AHe_relative_sample_distances = dfhs['distance_to_fault'].values
            U_conc = dfhs['U'].values
            Th_conc = dfhs['Th'].values
            sphere_radius = dfhs['sphere_radius'].values
            sample_names = dfhs['sample'].tolist()
            ngrains = len(AHe_sample_distances)
            unique_AHe_samples = np.unique(AHe_sample_names)
            n_ahe_samples = len(unique_AHe_samples)

        start_age = mp.t0 - runtimes[-1]

        t_prov = np.linspace(0, mp.t0, 31)
        T_prov = np.linspace(mp.T0, mp.T_surface, 31)

        print 'calculating helium ages'

        if mp.model_thermochron_surface is True:

            z_borehole, he_ages_borehole, he_ages_borehole_corrected = None, None, None

            xs_AHe_surface_all = []
            Ahe_ages_surface_all = []
            T_surf_mod_all = []

            if mp.model_AHe_samples is True:
                AHe_ages_surface_samples_all = []
                AHe_ages_surface_samples_all_corr = []
            else:
                AHe_ages_surface_samples_all = None
                AHe_ages_surface_samples_all_corr = None

            ind_surface = np.where(xyz_array[:, 1] == 0)[0]
            nx = len(ind_surface)

            for target_depth in mp.target_zs:

                print 'modeling AHe for samples at surface level = %0.2f m' % target_depth

                z_tolerance = 0.01
                ind_surface = np.where(np.abs(xyz_array[:, 1] - target_depth) < z_tolerance)[0]
                nx = len(ind_surface)

                he_ages_surface = np.zeros((nt, nx))
                T_surf_mod = np.zeros((nt, nx))

                if mp.model_AHe_samples is True:
                    he_ages_grains = np.zeros((nt, ngrains))
                    he_ages_grains_corr = np.zeros((nt, ngrains))

                for xii in range(nx):

                    # reduce the number of timesteps for the AHe algorithm
                    runtimes_filtered = runtimes[::mp.AHe_timestep_reduction]
                    T_filtered = T_array[:, ind_surface[xii]][::mp.AHe_timestep_reduction]

                    # make sure the last timestep is also included in the new
                    # temperature and time arrays:
                    if runtimes_filtered[-1] != runtimes[-1]:
                        runtimes_filtered = np.append(runtimes_filtered, runtimes[-1:])
                        T_filtered = np.append(T_filtered, T_array[:, ind_surface[xii]][-1])

                    # t_he = np.concatenate((t_prov[:], t_prov[-1] + runtimes))
                    # T_he = np.concatenate((T_prov[:], T_array[:, ind_surface[xii]]))
                    t_he = np.concatenate((t_prov[:], t_prov[-1] + runtimes_filtered[1:]))
                    T_he = np.concatenate((T_prov[:], T_filtered[1:]))

                    nt_prov = len(t_prov)
                    # T_he *= 2

                    T_he += mp.Kelvin

                    # calculate the AHe age:
                    try:
                        assert len(t_he) == len(T_he)
                    except AssertionError:
                        msg = 'warning, length temperature and time arrays for AHe model are not equal'
                        raise ValueError(msg)

                    he_age_i = he.calculate_he_age_meesters_dunai_2002(
                        t_he, T_he,
                        mp.radius, mp.U238, mp.Th232,
                        D0=mp.D0,
                        Ea=mp.Ea,
                        alpha_ejection=mp.alpha_ejection,
                        stopping_distance=mp.stopping_distance,
                        method=mp.AHe_method,
                        n_eigenmodes=50,
                        log_omega_p=mp.log_omega_p,
                        log_phi_p=mp.log_phi_p,
                        Etrap=mp.Etrap,
                        ln_D0_L_div_a2=mp.ln_D0_L_div_a2,
                        E_L=mp.E_L
                        )

                    # copy AHe ages back into array with same length as the
                    # runtime and temperature arrays
                    he_ages_run_filtered = he_age_i[nt_prov:]
                    he_ages_unfiltered = np.interp(runtimes, runtimes_filtered[1:], he_ages_run_filtered)
                    he_ages_surface[:, xii] = he_ages_unfiltered

                # get surface locs and T
                xs = xyz_array[:, 0][ind_surface]
                sort_order = np.argsort(xs)
                xs = xs[sort_order]

                # sort the data in order of increasing x
                ind_surface2 = ind_surface[sort_order]

                for i in xrange(nt):
                    T_surf_mod[i, :] = T_array[i][ind_surface2]
                    he_ages_surface[i] = he_ages_surface[i][sort_order]

                Ahe_ages_surface_all.append(he_ages_surface)
                xs_AHe_surface_all.append(xs)
                T_surf_mod_all.append(T_surf_mod)

                ###############################
                # calculate ages of AHe samples
                ###############################
                if mp.model_AHe_samples is True:
                    print 'modeling the ages for %i AHe samples' % (n_ahe_samples)
                    unique_dist = np.unique(AHe_relative_sample_distances)

                    #z_tolerance = 0.01
                    #ind_surface = np.where(np.abs(xyz_array[:, 1] - target_depth) < z_tolerance)[0]

                    # find two nearest nodes for samples
                    x_nodes_unsorted = xyz_array[:, 0][ind_surface]
                    x_nodes_sort_order = np.argsort(x_nodes_unsorted)
                    x_nodes = x_nodes_unsorted[x_nodes_sort_order]

                    fault_ind = np.where(target_depth >= z_flt)[0][0]
                    x_flt_timestep = x_flt[fault_ind]

                    #for n_sample, rel_distance in enumerate(unique_dist):

                    for n_sample, sample in enumerate(unique_AHe_samples):

                        inds = dfhs['sample'] == sample

                        rel_distance = dfhs.loc[inds, 'distance_to_fault'][0]
                        distance = x_flt_timestep + rel_distance

                        ind_node_left = np.where(x_nodes < distance)[0][-1]
                        ind_node_right = np.where(x_nodes > distance)[0][0]
                        x_factor = (distance - x_nodes[ind_node_left]) / (x_nodes[ind_node_right] - x_nodes[ind_node_left])

                        # extract temperature history for both nodes:
                        t_hes = []
                        T_hes = []

                        for xii in (ind_node_left, ind_node_right):

                            x_ind = x_nodes_sort_order[xii]

                            # reduce the number of timesteps for the AHe algorithm
                            runtimes_filtered = runtimes[::mp.AHe_timestep_reduction]
                            T_filtered = T_array[:, ind_surface[x_ind]][::mp.AHe_timestep_reduction]

                            # make sure the last timestep is also included in the new
                            # temperature and time arrays:
                            if runtimes_filtered[-1] != runtimes[-1]:
                                runtimes_filtered = \
                                    np.append(runtimes_filtered, runtimes[-1:])
                                T_filtered = \
                                    np.append(T_filtered,
                                              T_array[:, ind_surface[x_ind]][-1])

                            t_he = np.concatenate((t_prov[:],
                                                   t_prov[-1] + runtimes_filtered[1:]))
                            T_he = np.concatenate((T_prov[:], T_filtered[1:]))

                            nt_prov = len(t_prov)

                            T_he += mp.Kelvin

                            t_hes.append(t_he)
                            T_hes.append(T_he)

                        if np.max(np.abs(t_hes[1] - t_hes[0])) > 0.0:
                            print 'warning, something wrong with taking time ' \
                                  'history of left and right nodes'

                        T_he = x_factor * T_hes[1] + (1 - x_factor) * T_hes[0]
                        t_he = t_hes[0]

                        # find AHe grain data
                        grain_inds = np.where(AHe_relative_sample_distances == rel_distance)[0]
                        print 'AHe grains: ', grain_inds
                        print 'distance to fault (m): ', rel_distance
                        print 'x coord of fault (m): ', x_flt_timestep
                        print 'absolute distance for layer at z= %0.2f , x = %0.2f m' % (target_depth, distance)
                        print 'total time = %0.2e yr' % (t_he[-1] / year)

                        for grain_ind in grain_inds:

                            try:
                                assert len(t_he) == len(T_he)
                            except AssertionError:
                                msg = 'warning, length temperature and time arrays for ' \
                                      'AHe model are not equal'
                                raise ValueError(msg)

                            # calculate AHe ages of two nearest nodes
                            he_age_i = \
                                he.calculate_he_age_meesters_dunai_2002(
                                    t_he, T_he,
                                    sphere_radius[grain_ind] * 1e-06,
                                    U_conc[grain_ind] * 1e-6,
                                    Th_conc[grain_ind] * 1e-6,
                                    D0=mp.D0,
                                    Ea=mp.Ea,
                                    alpha_ejection=mp.alpha_ejection,
                                    stopping_distance=mp.stopping_distance,
                                    method=mp.AHe_method,
                                    n_eigenmodes=50,
                                    log_omega_p=mp.log_omega_p,
                                    log_phi_p=mp.log_phi_p,
                                    Etrap=mp.Etrap,
                                    ln_D0_L_div_a2=mp.ln_D0_L_div_a2,
                                    E_L=mp.E_L
                                    )

                            R = sphere_radius[grain_ind] * 1e-06
                            S = mp.stopping_distance
                            Ft_i = 1 - 3 * S / (4 * R) + S**3 / (16 * R**3)

                            he_age_i_corr = he_age_i / Ft_i

                            My = 1e6 * 365 * 24 * 60 * 60.0

                            print 'min, mean, max T = %0.2f, %0.2f, %0.2f C' % (T_he.min() - 273.15,
                                                                                T_he.mean() - 273.15,
                                                                                T_he.max() - 273.15)
                            print 'modeled AHe age uncorr = %0.2f Ma' % (he_age_i[-1] / My)
                            print 'modeled AHe age corr = %0.2f Ma' % (he_age_i_corr[-1] / My)

                            # copy AHe ages back into array with same length as the
                            # runtime and temperature arrays
                            he_ages_run_filtered = he_age_i[nt_prov:]
                            he_ages_unfiltered = np.interp(runtimes, runtimes_filtered[1:],
                                                           he_ages_run_filtered)
                            he_ages_grains[:, grain_ind] = he_ages_unfiltered

                            he_ages_run_filtered = he_age_i_corr[nt_prov:]
                            he_ages_unfiltered = np.interp(runtimes, runtimes_filtered[1:],
                                                           he_ages_run_filtered)
                            he_ages_grains_corr[:, grain_ind] = he_ages_unfiltered

                    AHe_ages_surface_samples_all.append(he_ages_grains)
                    AHe_ages_surface_samples_all_corr.append(he_ages_grains_corr)

            # calculate corrected ages, for surface data
            R = mp.radius
            S = mp.stopping_distance
            Ft_surf = 1 - 3 * S / (4 * R) + S ** 3 / (16 * R ** 3)

            #
            Ahe_ages_surface_corr_all = [age / Ft_surf for age in Ahe_ages_surface_all]

            Ahe_ages_all_output = Ahe_ages_surface_all

        elif mp.model_thermochron_borehole is True:

            Ahe_ages_surface_all = None
            Ahe_ages_surface_corr_all = None
            AHe_ages_surface_samples_all = None
            AHe_ages_surface_samples_all_corr = None
            xs_AHe_surface_all = None

            #Ahe_ages_surface_corr_all, xs_AHe_surface_all

            # find borehole nodes
            #x_tolerance = 0.01
            borehole_x = mp.borehole_xs[0]
            lowest_ahe_sample = -np.max(dfhs['depth'])

            print 'modeling AHe for samples at borehole loc = %0.2f m up to a depth of %0.2f m' \
                  % (borehole_x, lowest_ahe_sample)
            #ind_borehole_x = np.abs(xyz_array[:, 0] - borehole_x) < x_tolerance

            #ind_borehole_y = (xyz_array[:, 1] < 0) & (xyz_array[:, 1] > lowest_ahe_sample)
            #ind_borehole = np.where(ind_borehole_x & ind_borehole_y)[0]
            #ny = len(ind_borehole)

            x_buffer = 200.0
            pre_select = np.abs(xyz_array[:, 0] - borehole_x) < x_buffer

            # generate array with discretized borehole
            #z_obs1 = np.arange(lowest_ahe_sample, 0, mp.borehole_cellsize)
            z_borehole = np.append(np.arange(lowest_ahe_sample, 0, mp.borehole_cellsize), 0)
            ny = len(z_borehole)
            borehole_temp_modeled = np.zeros((nt, ny))

            for j in range(nt):

                zi, Ti = interpolate_data(xyz_array[pre_select], T_array[j][pre_select], borehole_x, y_steps=300)

                T_mod = np.interp(z_borehole, zi, Ti)

                borehole_temp_modeled[j, :] = T_mod

            he_ages_borehole = np.zeros((nt, ny))
            T_borehole_mod = np.zeros((nt, ny))

            #if mp.model_AHe_samples is True:
            #    he_ages_grains = np.zeros((nt, ngrains))
            #    he_ages_grains_corr = np.zeros((nt, ngrains))

            t_he_bh = []
            T_he_bh = []

            for yii in range(ny):

                # set up a pre-hydrothermal activity thermal history

                T_model_start = borehole_temp_modeled[0, yii]

                t_cooling_mod = mp.t_cooling + [0]

                T_cooling = np.zeros(len(t_cooling_mod))

                for ci in range(len(t_cooling_mod) - 1, 0, -1):
                    dt = t_cooling_mod[ci - 1] - t_cooling_mod[ci]
                    T_cooling[ci - 1] = T_cooling[ci] + mp.cooling_rates[ci - 1] * dt

                #print 'temperature difference before start of hydrothermal model'
                My = 1e6 * 365.24 * 24 * 3600.
                #for ti, Ti in zip(t_cooling_mod, T_cooling):
                #    print (ti / My), Ti

                T_cooling_mod1 = T_cooling + T_model_start
                T_cooling_mod = np.append(T_cooling_mod1, T_model_start)

                t_prov_bh = []
                T_prov_bh = []

                for tp1, tp2, Tp1, Tp2 in zip(t_cooling_mod[:-1], t_cooling_mod[1:],
                                              T_cooling_mod[:-1], T_cooling_mod[1:]):

                    t_prov_bh_i = np.linspace(tp1, tp2, 11)[:-1]
                    T_prov_bh_i = np.linspace(Tp1, Tp2, 11)[:-1]

                    t_prov_bh += list(t_prov_bh_i)
                    T_prov_bh += list(T_prov_bh_i)

                #T_start = T_model_start + mp.t0 * mp.cooling_rate
                #T_prov_bh = np.linspace(T_start, T_model_start, 31)

                # correct times to increasing order
                t_prov_bh = t_prov_bh[0] - np.array(t_prov_bh)

                # substract model run time
                t_prov_bh -= runtimes[-1]

                # reduce the number of timesteps for the AHe algorithm
                runtimes_filtered = runtimes[::mp.AHe_timestep_reduction]
                bhT_filtered = borehole_temp_modeled[::mp.AHe_timestep_reduction]

                # make sure the last timestep is also included in the new
                # temperature and time arrays:
                if runtimes_filtered[-1] != runtimes[-1]:
                    runtimes_filtered = np.append(runtimes_filtered, runtimes[-1:])
                    bhT_filtered = np.append(bhT_filtered, T_borehole_mod[-1])

                # t_he = np.concatenate((t_prov[:], t_prov[-1] + runtimes))
                # T_he = np.concatenate((T_prov[:], T_array[:, ind_surface[xii]]))
                t_he = np.concatenate((t_prov_bh[:], t_prov_bh[-1] + runtimes_filtered[1:]))
                T_he = np.concatenate((T_prov_bh[:], bhT_filtered[1:, yii]))

                T_he += mp.Kelvin

                t_he_bh.append(t_he)
                T_he_bh.append(T_he)

                nt_prov_bh = len(t_prov_bh)

                # calculate the AHe age:
                try:
                    assert len(t_he) == len(T_he)
                except AssertionError:
                    msg = 'warning, length temperature and time arrays for AHe model are not equal'
                    raise ValueError(msg)

                he_age_i = he.calculate_he_age_meesters_dunai_2002(
                    t_he, T_he,
                    mp.radius, mp.U238, mp.Th232,
                    D0=mp.D0,
                    Ea=mp.Ea,
                    alpha_ejection=mp.alpha_ejection,
                    stopping_distance=mp.stopping_distance,
                    method=mp.AHe_method,
                    n_eigenmodes=50,
                    log_omega_p=mp.log_omega_p,
                    log_phi_p=mp.log_phi_p,
                    Etrap=mp.Etrap,
                    ln_D0_L_div_a2=mp.ln_D0_L_div_a2,
                    E_L=mp.E_L
                    )

                # copy AHe ages back into array with same length as the
                # runtime and temperature arrays
                he_ages_run_filtered = he_age_i[nt_prov_bh:]
                he_ages_unfiltered = np.interp(runtimes, runtimes_filtered[1:], he_ages_run_filtered)
                he_ages_borehole[:, yii] = he_ages_unfiltered

            ####################
            # go through samples
            ####################
            n_grains_total = len(dfhs)
            he_ages_samples_borehole = np.zeros((nt, n_grains_total))

            grain_count = 0

            for n_sample, sample in enumerate(unique_AHe_samples):

                # find sample depth
                ind_sample = dfhs['sample'] == sample
                sample_depth = (dfhs.loc[ind_sample, 'depth'].values)[0]

                # select overlying and underlying thermal history
                if -sample_depth in z_borehole:
                    ind_node = np.where(z_borehole == -sample_depth)[0][0]
                    T_he_int = T_he_bh[ind_node]
                    t_he_int = t_he_bh[ind_node]
                else:
                    ind_node_top = np.where(z_borehole < -sample_depth)[0][-1]
                    ind_node_bottom = np.where(z_borehole > -sample_depth)[0][0]
                    x_factor = (sample_depth - z_borehole[ind_node_top]) / (z_borehole[ind_node_bottom] - z_borehole[ind_node_top])

                    # get average of two thermal histories
                    T_he_int = x_factor * T_he_bh[ind_node_bottom] + (1 - x_factor) * T_he_bh[ind_node_top]
                    t_he_int = t_he_bh[ind_node_top]

                # find n grains
                n_grains_sample = len(dfhs.loc[ind_sample])
                grains_in_sample = dfhs.loc[ind_sample, 'grain_number']

                # go through all the grains
                for grain_number, grain_name in enumerate(grains_in_sample):

                    ind_grain = dfhs['grain_number'] == grain_name

                    try:
                        assert len(dfhs.loc[ind_grain]) == 1
                    except AssertionError:
                        print 'error, either more than one or no grain with the name %s' % grain_name
                        print dfhs.loc[ind_grain]
                        pdb.set_trace()

                    he_age_i = he.calculate_he_age_meesters_dunai_2002(
                        t_he_int, T_he_int,
                        dfhs.loc[ind_grain, 'sphere_radius'].values[0] * 1e-6,
                        dfhs.loc[ind_grain, 'U'].values[0] * 1e-6,
                        dfhs.loc[ind_grain, 'Th'].values[0] * 1e-6,
                        D0=mp.D0,
                        Ea=mp.Ea,
                        alpha_ejection=mp.alpha_ejection,
                        stopping_distance=mp.stopping_distance,
                        method=mp.AHe_method,
                        n_eigenmodes=50,
                        log_omega_p=mp.log_omega_p,
                        log_phi_p=mp.log_phi_p,
                        Etrap=mp.Etrap,
                        ln_D0_L_div_a2=mp.ln_D0_L_div_a2,
                        E_L=mp.E_L
                        )

                    # copy AHe ages back into array with same length as the
                    # runtime and temperature arrays
                    he_ages_run_filtered = he_age_i[nt_prov_bh:]
                    he_ages_unfiltered = np.interp(runtimes, runtimes_filtered[1:], he_ages_run_filtered)
                    he_ages_samples_borehole[:, grain_count] = he_ages_unfiltered

                    grain_count += 1

            # calculate corrected ages, for borehole data
            R = mp.radius
            S = mp.stopping_distance
            Ft_bh = 1 - 3 * S / (4 * R) + S ** 3 / (16 * R ** 3)

            #
            he_ages_borehole_corrected = he_ages_borehole / Ft_bh


        #if mp.report_corrected_AHe_ages is True:
        #    Ahe_ages_all_output = Ahe_ages_corr_all

        print 'done calculating helium ages'

        My = 1e6 * 365 * 24 * 60 * 60.

        if mp.model_thermochron_surface is True:

            print '\nAHe ages: '
            for i, Ahe_ages in enumerate(Ahe_ages_surface_all):
                try:
                    print '\tlayer %i, min = %0.2f, mean = %0.2f, max = %0.2f My' \
                      % (i, Ahe_ages.min() / My,
                         Ahe_ages.mean() / My,
                         Ahe_ages.max() / My)
                except:
                    print '\tlayer %i, no modeled AHe data available'

            print '\nmodeled AHe ages samples'
            if mp.model_AHe_samples is True:
                print '\tname, distance, layer, modeled AHe age uncorr, corrected: '
                for i, age_i, age_i_corr in zip(itertools.count(), AHe_ages_surface_samples_all, AHe_ages_surface_samples_all_corr):
                    for sample_name, rel_distance, age, age_corr in \
                            zip(sample_names, AHe_relative_sample_distances,
                                age_i[-1], age_i_corr[-1]):
                        print '\t%s, %0.1f m, %i, %0.2f My, %0.2f My' \
                              % (sample_name, rel_distance, i, age / My, age_corr / My)

    print 'surface T: ', T * surface

    output = [runtimes, xyz_array, surface_levels, x_flt, z_flt,
              Ts, q_vectors, T_init_array, T_array, boiling_temp_array,
              xyz_array_exc, exceed_boiling_temp_array,
              xyz_element_array,
              qh_array, qv_array,
              fault_fluxes_m_per_sec, mp.durations, xzs, Tzs, Tzs_diff,
              Ahe_ages_surface_all, Ahe_ages_surface_corr_all, xs_AHe_surface_all, mp.target_zs,
              AHe_ages_surface_samples_all, AHe_ages_surface_samples_all_corr,
              z_borehole, he_ages_borehole, he_ages_borehole_corrected]

    return output



