"""
Mesh generation for the fipy backend of Beo.

This module recreates the model geometry of Beo using the gmsh Python API
instead of esys.pycad, so that a mesh can be generated without escript being
installed. The resulting gmsh file is read by fipy.

The geometry mirrors setup_mesh_with_exhumation_new in beo_core.py: a series
of stacked horizontal layers, each split into a left block, a fault block and
a right block, with an additional single block at the base of the model
domain. The horizontal lines separating the layers are located at the
elevations in target_zs, which allows the land surface to be lowered
step by step to simulate exhumation.
"""

import numpy as np

FIPY_AVAILABLE = False
try:
    import fipy
    FIPY_AVAILABLE = True
except ImportError:
    fipy = None

GMSH_AVAILABLE = False
try:
    import gmsh
    GMSH_AVAILABLE = True
except ImportError:
    gmsh = None


def check_mesh_dependencies():

    """
    Check that fipy and the gmsh python module are both installed
    """

    if FIPY_AVAILABLE is False:
        raise ImportError('error, the fipy backend of Beo requires fipy. '
                          'install this with pip install fipy')

    if GMSH_AVAILABLE is False:
        raise ImportError('error, the fipy backend of Beo requires the gmsh '
                          'python module. install this with pip install gmsh')


def setup_mesh_with_exhumation_fipy(width, x_flt_surface, fault_width, fault_angle, fault_bottom,
                                    z_air,
                                    z_fine, z_base, cellsize,
                                    cellsize_air, cellsize_surface,
                                    cellsize_fault_surface, cellsize_fault,
                                    cellsize_fine, cellsize_base, target_zs,
                                    discretize_borehole, borehole_xs,
                                    borehole_depths, borehole_cellsize,
                                    mesh_filename='beo_mesh_fipy.msh'):

    """
    Create a fipy mesh of the model domain, including horizontal lines at each
    exhumation level.

    This is the fipy equivalent of setup_mesh_with_exhumation_new in beo_core.py
    and takes the same arguments, with the addition of the filename that the
    gmsh mesh file is written to.

    Returns the fipy mesh, the x and z coordinates of the fault and the
    elevation of the base of the model domain.
    """

    check_mesh_dependencies()

    zs_surface = target_zs[::-1]

    z_flt = np.concatenate((zs_surface, np.array([z_fine, fault_bottom])))

    x_flt = (-z_flt) * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface

    print('x, z coords of fault locations in mesh:')
    for x, z in zip(x_flt, z_flt):
        print(x, z)

    x_left_bnd = np.min(x_flt) - width
    x_right_bnd = np.max(x_flt) + width

    # assemble the corner points of each horizontal layer in the model domain,
    # the first row is the top of the air layer, the last row is the base of
    # the model domain
    xys = [[[x_left_bnd, z_air], [x_flt[0] - fault_width / 2.0, z_air],
            [x_flt[0] + fault_width / 2.0, z_air], [x_right_bnd, z_air]]]

    for xf, zf in zip(x_flt, z_flt):
        xys.append([[x_left_bnd, zf],
                    [xf - fault_width / 2.0, zf],
                    [xf + fault_width / 2.0 + 0.02, zf],
                    [x_right_bnd, zf]])
    xys.append([[x_left_bnd, z_base], [x_right_bnd, z_base]])

    print('corner points in mesh: ')
    for xysi in xys:
        print(xysi)

    # assign a target cell size to each corner point, following the same order
    # of assignment as the escript version of this function so that later
    # assignments override earlier ones
    cellsizes = [[cellsize] * len(xysi) for xysi in xys]

    # cellsize in air layer
    for i in range(len(cellsizes[0])):
        cellsizes[0][i] = cellsize_air

    # cellsize at land surface
    for row in cellsizes[1:-3]:
        for i in range(len(row)):
            row[i] = cellsize_surface

    # cellsize in layer close to the surface
    for i in range(len(cellsizes[-3])):
        cellsizes[-3][i] = cellsize_fine

    # cellsize at fault surface
    for row in cellsizes[:-3]:
        row[1] = cellsize_fault_surface
        row[2] = cellsize_fault_surface

    # cellsize in fault
    for row in cellsizes[-3:-1]:
        row[1] = cellsize_fault
        row[2] = cellsize_fault

    # cellsize at lower left and right bnds
    cellsizes[-2][0] = cellsize_base
    cellsizes[-2][-1] = cellsize_base
    cellsizes[-1][0] = cellsize_base
    cellsizes[-1][-1] = cellsize_base

    gmsh.initialize()

    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add('beo_mesh')

        geo = gmsh.model.geo

        # add the corner points
        points = []
        for xy_row, cellsize_row in zip(xys, cellsizes):
            points_local = [geo.addPoint(x, z, 0.0, cs)
                            for (x, z), cs in zip(xy_row, cellsize_row)]
            points.append(points_local)

        # horizontal lines, three per layer boundary except for the base of the
        # model domain, which consists of a single line
        hlines = []
        for point in points[0:-1]:
            hline_local = [geo.addLine(point[0], point[1]),
                           geo.addLine(point[1], point[2]),
                           geo.addLine(point[2], point[3])]
            hlines.append(hline_local)

        hlines.append([geo.addLine(points[-1][0], points[-1][1])])

        # vertical lines connecting each pair of layer boundaries
        vlines = []
        for point, point_below in zip(points[0:-2], points[1:]):
            vline_local = [geo.addLine(point[0], point_below[0]),
                           geo.addLine(point[1], point_below[1]),
                           geo.addLine(point[2], point_below[2]),
                           geo.addLine(point[3], point_below[3])]
            vlines.append(vline_local)

        # the two sides of the block at the base of the model domain
        vlines.append([geo.addLine(points[-2][0], points[-1][0]),
                       geo.addLine(points[-2][3], points[-1][1])])

        # close each layer into a left, fault and right surface
        surfaces = []
        for hline, hline_below, vline in zip(hlines[0:-2], hlines[1:-1], vlines[:]):
            loop_left = geo.addCurveLoop([hline[0], vline[1], -hline_below[0], -vline[0]])
            loop_fault = geo.addCurveLoop([hline[1], vline[2], -hline_below[1], -vline[1]])
            loop_right = geo.addCurveLoop([hline[2], vline[3], -hline_below[2], -vline[2]])

            surfaces.append(geo.addPlaneSurface([loop_left]))
            surfaces.append(geo.addPlaneSurface([loop_fault]))
            surfaces.append(geo.addPlaneSurface([loop_right]))

        loop_bottom = geo.addCurveLoop([hlines[-2][0], hlines[-2][1], hlines[-2][2],
                                        vlines[-1][1], -hlines[-1][0], -vlines[-1][0]])
        surface_bottom = geo.addPlaneSurface([loop_bottom])
        surfaces.append(surface_bottom)

        geo.synchronize()

        # add extra vertical lines at the borehole locations to refine the mesh
        # there. each line is split at the layer boundaries it crosses and the
        # resulting segments are embedded in the surface that contains them
        if discretize_borehole is True:
            _embed_borehole_lines(borehole_xs, borehole_depths, borehole_cellsize,
                                  xys, x_flt, fault_width, max(target_zs), surfaces)

        # write a gmsh file in a format that fipy can read
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(2)
        gmsh.write(mesh_filename)

    finally:
        gmsh.finalize()

    print('reading mesh %s with fipy' % mesh_filename)
    mesh = fipy.Gmsh2D(mesh_filename, communicator=fipy.solvers.serialComm)

    print('fipy mesh with %i cells' % mesh.numberOfCells)

    return mesh, x_flt[:-2], z_flt[:-2], z_base


def _embed_borehole_lines(borehole_xs, borehole_depths, borehole_cellsize,
                          xys, x_flt, fault_width, z_top, surfaces):

    """
    Add vertical lines at the borehole locations and embed these in the
    surfaces of the model domain, so that the mesh is refined around the
    boreholes.

    The lines are cut at each horizontal layer boundary that they cross,
    because a line can only be embedded in a single surface.
    """

    geo = gmsh.model.geo

    # elevations of the horizontal layer boundaries, from top to bottom
    z_layers = np.array([xy_row[0][1] for xy_row in xys])

    embedded_lines = {}

    for x_bh, depth in zip(borehole_xs, borehole_depths):

        z_bottom = -depth

        for i in range(len(z_layers) - 1):

            z_upper = z_layers[i]
            z_lower = z_layers[i + 1]

            # find the part of the borehole that falls inside this layer
            z_start = min(z_upper, z_top)
            z_end = max(z_lower, z_bottom)

            if z_end >= z_start:
                continue

            # skip layers that the borehole does not reach
            if z_start <= z_bottom or z_end >= z_top:
                continue

            surface_tag = _find_surface(i, x_bh, x_flt, fault_width, xys, surfaces)

            if surface_tag is None:
                continue

            # keep a small distance from the layer boundaries, embedded lines
            # cannot touch the edges of a surface
            buffer = min(0.01 * (z_start - z_end), 1.0)
            p_top = geo.addPoint(x_bh, z_start - buffer, 0.0, borehole_cellsize)
            p_bottom = geo.addPoint(x_bh, z_end + buffer, 0.0, borehole_cellsize)
            bh_line = geo.addLine(p_top, p_bottom)

            embedded_lines.setdefault(surface_tag, []).append(bh_line)

    if len(embedded_lines) == 0:
        return

    geo.synchronize()

    for surface_tag, line_tags in embedded_lines.items():
        try:
            gmsh.model.mesh.embed(1, line_tags, 2, surface_tag)
        except Exception as error:
            print('warning, could not refine the mesh at the borehole '
                  'locations in surface %i: %s' % (surface_tag, str(error)))


def _find_surface(layer_index, x_bh, x_flt, fault_width, xys, surfaces):

    """
    Find the gmsh surface tag of the block that contains x coordinate x_bh in
    layer layer_index. Returns None for the block at the base of the model
    domain, which is not subdivided.
    """

    # the last surface covers the base of the model domain and is not split
    # into a left, fault and right block
    if layer_index >= len(surfaces) // 3:
        return None

    # the fault is bounded by the second and third corner point of the layer
    # boundary below the layer, use the widest part of the fault in this layer
    fault_left = min(xys[layer_index][1][0], xys[layer_index + 1][1][0])
    fault_right = max(xys[layer_index][2][0], xys[layer_index + 1][2][0])

    if x_bh < fault_left:
        return surfaces[layer_index * 3]
    elif x_bh > fault_right:
        return surfaces[layer_index * 3 + 2]
    else:
        return surfaces[layer_index * 3 + 1]


def setup_radial_mesh_fipy(r_outer, fault_width, fault_bottom,
                           z_air, z_fine, z_base, cellsize,
                           cellsize_air, cellsize_surface,
                           cellsize_fault_surface, cellsize_fault,
                           cellsize_fine, cellsize_base, target_zs,
                           mesh_filename='beo_mesh_radial.msh'):

    """
    Create a fipy mesh for a radially symmetric model of a vertical fault.

    The model domain is a vertical cross section from the symmetry axis at
    r = 0 to an outer radius r_outer. The fault is a cylindrical conduit
    centered on the symmetry axis with a radius of half the fault width, so
    that rotating the conduit around the axis gives a fault of the same width
    as in the cartesian version of the model.

    Unlike the cartesian mesh, each horizontal layer is divided into two
    blocks instead of three, the fault conduit and the rock surrounding it.
    There is no block to the left of the fault, because there is nothing on
    the other side of the symmetry axis.

    Returns the fipy mesh, the r and z coordinates of the fault and the
    elevation of the base of the model domain.
    """

    check_mesh_dependencies()

    fault_radius = fault_width / 2.0

    if fault_radius <= 0:
        raise ValueError('error, the fault width should be larger than zero '
                         'for a radially symmetric model')

    if r_outer <= fault_radius:
        raise ValueError('error, the outer radius of the model domain (%0.1f m) '
                         'should be larger than the radius of the fault '
                         '(%0.1f m)' % (r_outer, fault_radius))

    z_flt = np.concatenate((target_zs[::-1], np.array([z_fine, fault_bottom])))

    # each row of corner points runs from the symmetry axis, over the outer
    # edge of the fault conduit, to the outer boundary of the model domain
    xys = [[[0.0, z_air], [fault_radius, z_air], [r_outer, z_air]]]

    for zf in z_flt:
        xys.append([[0.0, zf], [fault_radius, zf], [r_outer, zf]])

    xys.append([[0.0, z_base], [r_outer, z_base]])

    print('radially symmetric mesh, fault conduit radius %0.2f m, '
          'outer radius %0.1f m' % (fault_radius, r_outer))
    print('corner points in mesh: ')
    for xysi in xys:
        print(xysi)

    # assign a target cell size to each corner point, following the same order
    # of assignment as the cartesian version of this function
    cellsizes = [[cellsize] * len(xysi) for xysi in xys]

    # cellsize in air layer
    for i in range(len(cellsizes[0])):
        cellsizes[0][i] = cellsize_air

    # cellsize at land surface
    for row in cellsizes[1:-3]:
        for i in range(len(row)):
            row[i] = cellsize_surface

    # cellsize in layer close to the surface
    for i in range(len(cellsizes[-3])):
        cellsizes[-3][i] = cellsize_fine

    # cellsize at the fault, both at the symmetry axis and at the outer edge
    # of the conduit
    for row in cellsizes[:-3]:
        row[0] = cellsize_fault_surface
        row[1] = cellsize_fault_surface

    for row in cellsizes[-3:-1]:
        row[0] = cellsize_fault
        row[1] = cellsize_fault

    # cellsize at the outer boundary of the base of the model domain
    cellsizes[-2][-1] = cellsize_base
    cellsizes[-1][-1] = cellsize_base

    gmsh.initialize()

    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add('beo_radial_mesh')

        geo = gmsh.model.geo

        points = []
        for xy_row, cellsize_row in zip(xys, cellsizes):
            points_local = [geo.addPoint(r, z, 0.0, cs)
                            for (r, z), cs in zip(xy_row, cellsize_row)]
            points.append(points_local)

        # horizontal lines, two per layer boundary except for the base of the
        # model domain, which consists of a single line
        hlines = []
        for point in points[0:-1]:
            hlines.append([geo.addLine(point[0], point[1]),
                           geo.addLine(point[1], point[2])])

        hlines.append([geo.addLine(points[-1][0], points[-1][1])])

        # vertical lines connecting each pair of layer boundaries
        vlines = []
        for point, point_below in zip(points[0:-2], points[1:]):
            vlines.append([geo.addLine(point[0], point_below[0]),
                           geo.addLine(point[1], point_below[1]),
                           geo.addLine(point[2], point_below[2])])

        vlines.append([geo.addLine(points[-2][0], points[-1][0]),
                       geo.addLine(points[-2][2], points[-1][1])])

        # close each layer into a fault conduit and a surrounding rock surface
        surfaces = []
        for hline, hline_below, vline in zip(hlines[0:-2], hlines[1:-1], vlines[:]):
            loop_fault = geo.addCurveLoop([hline[0], vline[1], -hline_below[0], -vline[0]])
            loop_rock = geo.addCurveLoop([hline[1], vline[2], -hline_below[1], -vline[1]])

            surfaces.append(geo.addPlaneSurface([loop_fault]))
            surfaces.append(geo.addPlaneSurface([loop_rock]))

        loop_bottom = geo.addCurveLoop([hlines[-2][0], hlines[-2][1],
                                        vlines[-1][1], -hlines[-1][0], -vlines[-1][0]])
        surfaces.append(geo.addPlaneSurface([loop_bottom]))

        geo.synchronize()

        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(2)
        gmsh.write(mesh_filename)

    finally:
        gmsh.finalize()

    print('reading mesh %s with fipy' % mesh_filename)
    mesh = fipy.Gmsh2D(mesh_filename, communicator=fipy.solvers.serialComm)

    print('fipy mesh with %i cells' % mesh.numberOfCells)

    r_flt = np.zeros(len(z_flt))

    return mesh, r_flt[:-2], z_flt[:-2], z_base
