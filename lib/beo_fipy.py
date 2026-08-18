"""
fipy backend for Beo.

This module solves the same heat advection and conduction equation as the
escript backend of Beo, but uses the finite volume code fipy instead of the
finite element code escript. The governing equation is:

rho c dT/dt = div(K grad T) - rho_f c_f q . grad T

with T temperature, rho c the bulk volumetric heat capacity, K the bulk
thermal conductivity, rho_f c_f the volumetric heat capacity of water and q
the advective fluid flux in the fault zones and aquifers.

Note that the advection term is written in non-conservative form, ie as
q . grad T and not as div(q T). The advective flux in Beo is prescribed and
is not divergence free, since the flux is switched on abruptly at the edges
of the fault zones and aquifers. Using the conservative form would therefore
introduce spurious heat sources and sinks at the fault walls. fipy provides
the conservative form only, so the equation is assembled as

q . grad T = div(q T) - T div(q)

with the second term added as an implicit source term.

Temperatures are calculated at the cell centers and are interpolated to the
mesh vertices for the model output, because the post-processing code of Beo
extracts temperatures at exact elevations, which are present in the mesh as
horizontal lines but not as cell centers.
"""

import time
import itertools
import hashlib

import numpy as np
from scipy.sparse.linalg import splu

import lib.mesh_functions_fipy as mesh_functions_fipy

FIPY_AVAILABLE = False
try:
    import fipy
    from fipy import CellVariable, FaceVariable
    from fipy import TransientTerm, DiffusionTerm, ImplicitSourceTerm
    from fipy import PowerLawConvectionTerm, ExponentialConvectionTerm
    from fipy import UpwindConvectionTerm, CentralDifferenceConvectionTerm
    FIPY_AVAILABLE = True
except ImportError:
    fipy = None


day = 24.0 * 60.0 * 60.0
year = 365.25 * day

convection_schemes = {'powerlaw': 'PowerLawConvectionTerm',
                      'exponential': 'ExponentialConvectionTerm',
                      'upwind': 'UpwindConvectionTerm',
                      'central': 'CentralDifferenceConvectionTerm'}


def check_fipy_available():

    """
    Check whether both fipy and the gmsh python module are installed
    """

    return FIPY_AVAILABLE and mesh_functions_fipy.GMSH_AVAILABLE


def build_growing_dt_schedule(dt_initial, growth_factor, dt_max, duration):

    """
    Build a list of timestep sizes covering a single time period of length
    duration.

    The first timestep is dt_initial, and each following timestep is
    growth_factor times larger than the previous one, until the timestep
    reaches dt_max. From that point on the timestep stays constant at
    dt_max. The last timestep is shortened so that all the timesteps in the
    returned list add up to exactly duration.

    This is a plain numerical helper, with no dependency on fipy, so that it
    can also be used to precompute the exact number of stored timesteps in
    beo.py without needing the fipy backend to be installed.
    """

    if growth_factor <= 1.0:
        raise ValueError('error, dt_growth_factor should be larger than 1.0')
    if dt_max <= dt_initial:
        raise ValueError('error, dt_max should be larger than dt')

    dts = []
    t_total = 0.0
    dt = dt_initial

    while t_total < duration:
        dt = min(dt, dt_max)
        remaining = duration - t_total
        if dt > remaining:
            dt = remaining
        dts.append(dt)
        t_total += dt
        dt = dt * growth_factor

    return dts


def get_convection_term(scheme):

    """
    Return the fipy convection term class for a given discretization scheme
    """

    if scheme not in convection_schemes:
        msg = 'error, unknown convection scheme %s, choose one of %s' \
              % (str(scheme), str(list(convection_schemes.keys())))
        raise ValueError(msg)

    return {'powerlaw': PowerLawConvectionTerm,
            'exponential': ExponentialConvectionTerm,
            'upwind': UpwindConvectionTerm,
            'central': CentralDifferenceConvectionTerm}[scheme]


def get_solver():

    """
    Return a fipy solver that decides convergence by how far the residual has
    been reduced compared to the residual at the start of the solve.

    Newer versions of fipy compare the residual with the norm of the right hand
    side of the equation instead. That criterion fails for radially symmetric
    models, where the coefficients of the equation are multiplied by the
    radius. The right hand side is then dominated by the cells far away from
    the symmetry axis, while the interesting part of the solution is close to
    the axis. The solver decides that the equation is already solved and
    returns the temperatures of the previous timestep unchanged.
    """

    if FIPY_AVAILABLE is False:
        return None

    try:
        return fipy.DefaultSolver(criterion='initial')
    except (TypeError, ValueError):
        # older versions of fipy do not have a criterion argument
        return fipy.DefaultSolver()


class TransientLUCache:

    """
    Reuse the sparse LU factorization of the heat flow equation.

    For a transient model without exhumation, a variable thermal conductivity
    of the air layer or a vapour correction, the coefficients of the heat flow
    equation do not change during a time period. The matrix of the linear
    system is then the same at every timestep, while the right hand side is a
    linear function of the temperature of the previous timestep,

        b = b_const + a * T_old

    The matrix is assembled and factorized once and the factorization is reused
    for the following timesteps. The coefficient a is found by assembling the
    right hand side for two different temperatures, see the docstring of
    _measure_rhs_coefficient. The matrix is assembled and factorized again
    whenever any of the coefficients of the equation changes, so that models in
    which the equation changes over time give the same result as before.
    """

    def __init__(self, validate=False):

        self.validate = validate

        self.factor = None
        self.b_const = None
        self.diagonal = None
        self.key = None

        self.n_factorizations = 0
        self.n_reused = 0
        self.max_rhs_error = 0.0

    def _key(self, coefficient_vars, dt):

        parts = []

        for var in coefficient_vars:
            values = np.ascontiguousarray(np.array(var), dtype=float)
            parts.append(hashlib.blake2b(values.view(np.uint8),
                                         digest_size=16).digest())

        parts.append(repr(dt).encode())

        return tuple(parts)

    def _assemble(self, equation, T, dt):

        equation.cacheMatrix()
        equation.cacheRHSvector()

        solver = equation._prepareLinearSystem(var=T, solver=get_solver(),
                                               boundaryConditions=(), dt=dt)

        matrix = solver.matrix.matrix.asformat('csc')
        rhs = np.array(solver.RHSvector, dtype=float)

        return matrix, rhs

    def _assemble_with_temperature(self, equation, T, dt, values):

        """
        Assemble the right hand side for a given temperature field and restore
        the temperature afterwards
        """

        T_current = np.array(T, dtype=float)
        T_old = np.array(T.old, dtype=float)

        T.setValue(values)
        T.updateOld()

        rhs = self._assemble(equation, T, dt)[1]

        T.setValue(T_old)
        T.updateOld()
        T.setValue(T_current)

        return rhs

    def _measure_rhs_coefficient(self, equation, T, dt, T_old):

        """
        Find the coefficient of the temperature of the previous timestep in the
        right hand side of the equation and return the right hand side for the
        current temperature.

        The right hand side is a linear function of the temperature of the
        previous timestep. The transient term contributes rho c V / dt, and
        fipy adds a second contribution for the cells where the coefficient of
        the implicit source term has the sign that would destabilize the
        matrix, see ImplicitSourceTerm._getWeight in fipy. Instead of
        reproducing that split here, the coefficient is measured by assembling
        the right hand side twice with temperatures that differ by one degree.
        """

        rhs = self._assemble(equation, T, dt)[1]
        rhs_probe = self._assemble_with_temperature(equation, T, dt, T_old + 1.0)

        self.diagonal = rhs_probe - rhs
        self.b_const = rhs - self.diagonal * T_old

        return rhs

    def solve(self, equation, T, dt, matrix_vars, rhs_vars):

        key = self._key(list(matrix_vars) + list(rhs_vars), dt)
        T_old = np.array(T.old, dtype=float)

        if self.factor is not None and key == self.key:

            if self.diagonal is None:
                # the coefficients of the equation have not changed since the
                # last timestep, so from here on the right hand side can be
                # reconstructed. the coefficient of the temperature of the
                # previous timestep is measured once, see the docstring of
                # _measure_rhs_coefficient
                rhs = self._measure_rhs_coefficient(equation, T, dt, T_old)
            else:
                rhs = self.b_const + self.diagonal * T_old

                if self.validate is True:
                    rhs_check = self._assemble(equation, T, dt)[1]
                    scale = np.max(np.abs(rhs_check))
                    error = np.max(np.abs(rhs - rhs_check)) / scale
                    self.max_rhs_error = max(self.max_rhs_error, error)

            T.setValue(self.factor.solve(rhs))
            self.n_reused += 1

            return

        matrix, rhs = self._assemble(equation, T, dt)

        # the coefficient of the right hand side is only measured once the
        # coefficients of the equation turn out to stay the same, so that a
        # model in which the equation changes at every timestep, for instance
        # a model with exhumation, does not do any extra work
        self.diagonal = None
        self.b_const = None

        self.factor = splu(matrix)
        self.key = key
        self.n_factorizations += 1

        T.setValue(self.factor.solve(rhs))


class OutputTimeSelector:

    """
    Decide during the timestep loop whether to store the model results of the
    current timestep.

    By default results are stored at a regular interval, once every dt_stored
    seconds of model time. Instead, results can be stored at a list of
    specific times, given by the stored_timesteps parameter. This is useful
    to resolve a part of the model run in more detail than the rest, for
    instance the first few timesteps after a fault zone starts discharging
    water, without having to run the entire model at a correspondingly small
    value of dt_stored.

    stored_timesteps is a single list of times in seconds, covering the
    entire model run, not one list per timeslice. Results are stored the
    first time the cumulative model time reaches or passes each value in the
    list, so a requested time that falls between two timesteps is rounded up
    to the next timestep, in the same way that dt_stored is already rounded
    to the nearest multiple of dt. Requested times that are closer together
    than dt cannot be resolved and are skipped, so that at most one result is
    stored per timestep.

    If the timestep dt grows during the run, the regular dt_stored interval
    is recomputed whenever dt changes, so that results keep being stored at
    roughly every dt_stored seconds of simulated time once dt has reached
    its final, constant value. No results are stored at the regular interval
    while dt is still growing, since the interval in timesteps would be
    different for every value of dt. Use stored_timesteps instead of
    dt_stored to resolve part of a run with a growing timestep in detail.
    """

    def __init__(self, dt, dt_stored, stored_timesteps):

        self.stored_timesteps = stored_timesteps
        self.dt_stored = dt_stored

        if stored_timesteps is not None:
            self.targets = sorted(set(stored_timesteps))
            self.target_index = 0
        else:
            self._interval_dt = None
            self._set_interval(dt)

    def _set_interval(self, dt):

        self.interval = int(round(self.dt_stored / dt))
        if self.interval < 1:
            self.interval = 1
        self.count = 0
        self._interval_dt = dt

    def should_store(self, t_total, dt):

        """
        Return True if the results of the current timestep, with cumulative
        model time t_total and timestep size dt, should be stored
        """

        if self.stored_timesteps is not None:

            store = False

            while (self.target_index < len(self.targets)
                   and t_total >= self.targets[self.target_index]):
                store = True
                self.target_index += 1

            return store

        else:
            if dt != self._interval_dt:
                self._set_interval(dt)

            self.count += 1
            if self.count >= self.interval:
                self.count = 0
                return True

            return False

    def report_actual_times(self, dt):

        """
        Print the times that were actually requested in stored_timesteps next
        to the time they will actually be stored at, so that a difference
        caused by snapping to a reachable timestep is not mistaken for an
        error.

        should_store only stores results once t_total has reached or passed a
        requested time, so the time actually stored at is always the closest
        multiple of dt that is greater than or equal to the requested time,
        not the closest multiple of dt overall.
        """

        if self.stored_timesteps is None:
            return

        print('storing results at the following times, snapped to the next '
              'reachable multiple of dt:')

        for target in self.targets:
            actual = dt * np.ceil(target / dt)
            print('\trequested %0.3e s, stored at %0.3e s' % (target, actual))


def calculate_vapour_pressure(T, c1=8e-8, c2=-7e-5, c3=0.028, c4=-3.1597):

    """
    Calculate the vapour pressure of water

    based on a 3rd order polynomial fit to vapour curve data by the NIST
    publication "Thermophysical Properties of Fluid Systems"
    """

    log_Pv = c1 * T**3 + c2 * T**2 + c3 * T + c4

    return 10**log_Pv * 1e6


def build_cell_to_vertex_weights(mesh):

    """
    Set up the connections between the cells and the vertices of a fipy mesh,
    which are used to interpolate values at the cell centers to the mesh
    vertices.

    Returns the vertex index, cell index, weight and coordinate offset of each
    cell to vertex connection, and the number of vertices.
    """

    cell_vertex_ids = mesh._cellVertexIDs
    cell_centers = np.array(mesh.cellCenters)
    vertex_coords = np.array(mesh.vertexCoords)

    n_vertices = vertex_coords.shape[1]
    n_cells = cell_centers.shape[1]

    if np.ma.isMaskedArray(cell_vertex_ids):
        valid = ~np.ma.getmaskarray(cell_vertex_ids)
        vertex_ids_all = np.ma.getdata(cell_vertex_ids)
    else:
        vertex_ids_all = np.asarray(cell_vertex_ids)
        valid = np.ones(vertex_ids_all.shape, dtype=bool)

    cell_ids_all = np.tile(np.arange(n_cells), (vertex_ids_all.shape[0], 1))

    vertex_ids = vertex_ids_all[valid].astype(int)
    cell_ids = cell_ids_all[valid].astype(int)

    # a second set of connections that also includes the neighbours of the
    # cells around each vertex. this wider stencil is only used for vertices
    # where the cells directly around the vertex do not define a plane, which
    # happens at the corners and along parts of the edges of the model domain
    vertex_ids_wide, cell_ids_wide = _add_neighbouring_cells(
        mesh, vertex_ids, cell_ids, n_cells)

    stencils = []
    for v_ids, c_ids in [(vertex_ids, cell_ids), (vertex_ids_wide, cell_ids_wide)]:
        # offset of each cell center relative to the vertex
        dx = cell_centers[0][c_ids] - vertex_coords[0][v_ids]
        dy = cell_centers[1][c_ids] - vertex_coords[1][v_ids]
        distance = np.sqrt(dx**2 + dy**2)
        weights = 1.0 / np.maximum(distance, 1e-10)

        stencils.append((v_ids, c_ids, weights, dx, dy))

    return {'n_vertices': n_vertices,
            'stencil': stencils[0],
            'wide_stencil': stencils[1]}


def _add_neighbouring_cells(mesh, vertex_ids, cell_ids, n_cells):

    """
    Extend the set of cells around each vertex with the cells that share a
    face with those cells.
    """

    face_cells = np.array(mesh.faceCellIDs)
    interior = ~np.array(mesh.exteriorFaces, dtype=bool)

    cell_a = np.concatenate((face_cells[0][interior], face_cells[1][interior]))
    cell_b = np.concatenate((face_cells[1][interior], face_cells[0][interior]))

    order = np.argsort(cell_a, kind='stable')
    cell_a = cell_a[order].astype(int)
    cell_b = cell_b[order].astype(int)

    n_neighbours = np.bincount(cell_a, minlength=n_cells)
    start = np.concatenate(([0], np.cumsum(n_neighbours)))

    repeats = n_neighbours[cell_ids]
    n_total = int(np.sum(repeats))

    group_start = np.concatenate(([0], np.cumsum(repeats)[:-1]))
    within_group = np.arange(n_total) - np.repeat(group_start, repeats)
    index = np.repeat(start[cell_ids], repeats) + within_group

    vertex_ids_wide = np.concatenate((vertex_ids, np.repeat(vertex_ids, repeats)))
    cell_ids_wide = np.concatenate((cell_ids, cell_b[index]))

    # remove the duplicate cell to vertex connections
    key = vertex_ids_wide.astype(np.int64) * n_cells + cell_ids_wide
    unique = np.unique(key)

    return (unique // n_cells).astype(int), (unique % n_cells).astype(int)


def interpolate_to_vertices(values, vertex_weights, cell_mask=None):

    """
    Interpolate values at the cell centers to the mesh vertices.

    A weighted least squares fit of a plane through the values of the cells
    surrounding each vertex is used. Unlike inverse distance weighting this is
    exact for a linear field, also at the boundary of the model domain, where
    all surrounding cell centers lie on one side of the vertex.

    Vertices where the surrounding cells do not define a plane are calculated
    again using a wider stencil that also includes the neighbours of those
    cells.

    If cell_mask is given only the cells in the mask contribute. Vertices
    without any contributing cell are filled in using all cells instead.
    """

    n_vertices = vertex_weights['n_vertices']

    values = np.asarray(values, dtype=float)

    result, found, solvable = _vertex_values_from_stencil(
        values, vertex_weights['stencil'], n_vertices, cell_mask)

    # calculate the remaining vertices again with a wider stencil
    if np.any(~solvable):
        result_wide, found_wide, solvable_wide = _vertex_values_from_stencil(
            values, vertex_weights['wide_stencil'], n_vertices, cell_mask)
        replace = (~solvable) & (solvable_wide | (found_wide & ~found))
        result[replace] = result_wide[replace]
        found = found | found_wide

    # fill in vertices that are not surrounded by any of the selected cells
    if cell_mask is not None and np.any(~found):
        result_all = interpolate_to_vertices(values, vertex_weights)
        result[~found] = result_all[~found]

    return result


def _vertex_values_from_stencil(values, stencil, n_vertices, cell_mask):

    """
    Apply the least squares fit to the cell to vertex connections of a single
    stencil, using only the cells in cell_mask
    """

    vertex_ids, cell_ids, weights, dx, dy = stencil

    if cell_mask is None:
        return _least_squares_vertex_values(values, vertex_ids, cell_ids,
                                            weights, dx, dy, n_vertices)

    selection = cell_mask[cell_ids]

    return _least_squares_vertex_values(
        values, vertex_ids[selection], cell_ids[selection], weights[selection],
        dx[selection], dy[selection], n_vertices)


def _least_squares_vertex_values(values, vertex_ids, cell_ids, weights,
                                 dx, dy, n_vertices):

    """
    Fit a plane through the values of the cells surrounding each vertex and
    return the value of that plane at the vertex.

    Also returns a mask of the vertices that have at least one surrounding
    cell and a mask of the vertices where the least squares fit could be
    solved.
    """

    value_at_cells = values[cell_ids]

    def total(w):
        return np.bincount(vertex_ids, weights=w, minlength=n_vertices)

    # normal equations of the weighted least squares fit
    m00 = total(weights)
    m01 = total(weights * dx)
    m02 = total(weights * dy)
    m11 = total(weights * dx * dx)
    m12 = total(weights * dx * dy)
    m22 = total(weights * dy * dy)

    b0 = total(weights * value_at_cells)
    b1 = total(weights * dx * value_at_cells)
    b2 = total(weights * dy * value_at_cells)

    found = m00 > 0

    result = np.zeros(n_vertices)

    # fall back to inverse distance weighting where there are not enough cells
    result[found] = b0[found] / m00[found]

    determinant = (m00 * (m11 * m22 - m12 * m12)
                   - m01 * (m01 * m22 - m12 * m02)
                   + m02 * (m01 * m12 - m11 * m02))

    # the determinant scales with the third power of the diagonal, use this to
    # find vertices where the least squares fit is ill conditioned
    scale = m00 * m11 * m22
    solvable = found & (np.abs(determinant) > 1e-8 * np.abs(scale))

    if np.any(solvable):
        # only the value at the vertex is needed, which is the first element
        # of the solution of the normal equations
        cofactor_00 = m11[solvable] * m22[solvable] - m12[solvable]**2
        cofactor_01 = m02[solvable] * m12[solvable] - m01[solvable] * m22[solvable]
        cofactor_02 = m01[solvable] * m12[solvable] - m02[solvable] * m11[solvable]

        result[solvable] = (cofactor_00 * b0[solvable]
                            + cofactor_01 * b1[solvable]
                            + cofactor_02 * b2[solvable]) / determinant[solvable]

    return result, found, solvable


def calculate_fault_x(z, fault_angle, x_flt_surface):

    """
    Calculate the x coordinate of a fault at elevation z
    """

    return -z * np.tan(np.deg2rad(90 - fault_angle)) - 0.01 + x_flt_surface


def build_material_properties(x, y, mp, fault_angle, fault_x_surface):

    """
    Calculate the thermal conductivity, density and heat capacity of the rock
    matrix at coordinates x, y, following the layers specified in the model
    parameter file.
    """

    fault_x = calculate_fault_x(y, fault_angle, fault_x_surface)

    K_solid = np.zeros_like(x)
    porosity = np.zeros_like(x)

    try:
        assert len(mp.layer_bottom) == len(mp.K_solids) == len(mp.porosities)
    except AssertionError:
        msg = 'error, number of layers, K_solid and porosities supplied in input file are not equal'
        msg += '\nlayers = %i, K_solid = %i, porosities = %i' \
               % (len(mp.layer_bottom), len(mp.K_solids), len(mp.porosities))
        msg += '\n check your input file'
        raise ValueError(msg)

    left = x <= fault_x
    right = x > fault_x

    for i, layer_bottom_i in enumerate(mp.layer_bottom):
        in_layer_left = left & (y > layer_bottom_i[0])
        in_layer_right = right & (y > layer_bottom_i[1])
        in_layer = in_layer_left | in_layer_right

        K_solid[in_layer] = mp.K_solids[i]
        porosity[in_layer] = mp.porosities[i]

    K_b = K_solid ** (1 - porosity) * mp.K_water ** porosity
    rho_b = porosity * mp.rho_f + (1 - porosity) * mp.rho_s
    c_b = porosity * mp.c_f + (1 - porosity) * mp.c_s

    return K_b, rho_b, c_b


def build_fault_and_aquifer_masks(x, y, mp, surface_level):

    """
    Locate the fault zones, the fault segments and the aquifers at coordinates
    x, y. The masks are returned as float arrays with a value of 1 inside a
    fault or aquifer and 0 elsewhere.
    """

    subsurface = (y <= surface_level).astype(float)

    fault_zones = []
    for fault_x, fault_angle, fault_width, fault_bottom \
            in zip(mp.fault_xs, mp.fault_angles, mp.fault_widths, mp.fault_bottoms):

        fault_left = calculate_fault_x(y, fault_angle, fault_x) - fault_width / 2.0
        fault_right = fault_left + fault_width + 0.02

        fault_zone = (subsurface * (x > fault_left) * (x < fault_right)
                      * (y > fault_bottom))

        fault_zones.append(fault_zone.astype(float))

    fault_segments_all = None
    if mp.fault_segments is not None:
        fault_segments_all = []
        for fault_x, fault_angle, fault_width, fault_bottom, fault_segments \
                in zip(mp.fault_xs, mp.fault_angles, mp.fault_widths,
                       mp.fault_bottoms, mp.fault_segments):

            fault_left = calculate_fault_x(y, fault_angle, fault_x) - fault_width / 2.0
            fault_right = fault_left + fault_width + 0.02

            fault_segments_i = []
            segment_bottoms = [fault_bottom] + list(fault_segments[:-1])
            for segment_top, segment_bottom in zip(fault_segments, segment_bottoms):
                fault_zone_segment = (subsurface * (x > fault_left) * (x < fault_right)
                                      * (y > segment_bottom) * (y <= segment_top))
                fault_segments_i.append(fault_zone_segment.astype(float))

            fault_segments_all.append(fault_segments_i)

    aquifer_locs = []
    for aquifer_top, aquifer_bottom, aquifer_left_bnd, aquifer_right_bnd, aquifer_angle \
            in zip(mp.aquifer_tops, mp.aquifer_bottoms,
                   mp.aquifer_left_bnds, mp.aquifer_right_bnds, mp.aquifer_angles):

        if aquifer_top is None:
            continue

        try:
            assert aquifer_top > aquifer_bottom
        except AssertionError:
            raise AssertionError('error, the aquifer bottom is higher than the top, '
                                 'check input file')

        # the fault of the last iteration bounds the aquifer if no explicit
        # boundary is given, following the escript version of Beo
        aquifer_left_bnd_i = aquifer_left_bnd
        aquifer_right_bnd_i = aquifer_right_bnd

        if aquifer_left_bnd is None:
            aquifer_left_bnd_i = fault_right
        if aquifer_right_bnd is None:
            aquifer_right_bnd_i = fault_left

        aquifer_thickness = aquifer_top - aquifer_bottom
        aquifer_x = x - aquifer_left_bnd_i
        aquifer_top_i = aquifer_top + aquifer_x * np.tan(np.deg2rad(aquifer_angle))
        aquifer_bottom_i = aquifer_top_i - aquifer_thickness

        aquifer_loc = (subsurface * (y < aquifer_top_i) * (y > aquifer_bottom_i)
                       * (x < aquifer_right_bnd_i) * (x > aquifer_left_bnd_i))

        aquifer_loc = aquifer_loc.astype(float)

        # remove the overlap between the aquifer and the fault zone
        if fault_segments_all is not None:
            for fault_segments_i in fault_segments_all:
                for fault_segment_ii in fault_segments_i:
                    aquifer_loc = aquifer_loc * (fault_segment_ii == 0)

        aquifer_locs.append(aquifer_loc.astype(float))

    # cut off the fault zones above the level where they meet an aquifer.
    # this assumes that the aquifer drains the fault, so that there is no flow
    # in the fault above the aquifer. set aquifer_cuts_off_fault to False if
    # the aquifer feeds the fault instead, which is the case for an aquifer at
    # the base of a fault that supplies the water for a spring or a vent
    if getattr(mp, 'aquifer_cuts_off_fault', True) is True:
        for i, fault_zone in enumerate(fault_zones):
            for aquifer_top, aquifer_bottom in zip(mp.aquifer_tops, mp.aquifer_bottoms):
                if aquifer_top is None:
                    continue
                flow_cutoff_level = (aquifer_top + aquifer_bottom) / 2.0
                if np.any(fault_zone * (y > flow_cutoff_level)):
                    fault_zones[i] = fault_zones[i] * (y <= flow_cutoff_level)

    return fault_zones, fault_segments_all, aquifer_locs


def build_flux_field(x, y, fault_zones, fault_segments_all, fault_angles,
                     fault_flux, aquifer_locs, aquifer_fluxes, aquifer_angles,
                     surface_level, surface_buffer=0.1, radial=False):

    """
    Calculate the horizontal and vertical components of the advective fluid
    flux at coordinates x, y for a single time period.

    For a radially symmetric model the radial flux in an aquifer is divided by
    the circumference at that radius, because the flow crosses a cylinder that
    grows with the distance from the symmetry axis. The flux in the fault is
    not scaled, the flow there is vertical and is spread evenly over the cross
    section of the conduit.
    """

    qh = np.zeros_like(x)
    qv = np.zeros_like(x)

    for h, fault_zone, fault_angle, q_fault_zone in \
            zip(itertools.count(), fault_zones, fault_angles, fault_flux):

        if fault_segments_all is not None:
            for n_segment, fault_zone_segment, q_fault_zone_segment in \
                    zip(itertools.count(), fault_segments_all[h], q_fault_zone):
                qh_fault_zone = -q_fault_zone_segment * np.cos(np.deg2rad(fault_angle))
                qv_fault_zone = q_fault_zone_segment * np.sin(np.deg2rad(fault_angle))
                qh += fault_zone_segment * qh_fault_zone
                qv += fault_zone_segment * qv_fault_zone
        else:
            qh_fault_zone = -q_fault_zone * np.cos(np.deg2rad(fault_angle))
            qv_fault_zone = q_fault_zone * np.sin(np.deg2rad(fault_angle))
            qh += fault_zone * qh_fault_zone
            qv += fault_zone * qv_fault_zone

    for aquifer_loc, aquifer_flux, aquifer_angle in \
            zip(aquifer_locs, aquifer_fluxes, aquifer_angles):
        if radial is True:
            # divide the total flow by the circumference at each radius
            aquifer_flux = aquifer_flux / (2.0 * np.pi * np.maximum(x, 1e-6))
        qh += aquifer_loc * aquifer_flux * np.cos(np.deg2rad(aquifer_angle))
        qv += aquifer_loc * aquifer_flux * np.sin(np.deg2rad(aquifer_angle))

    # make sure there is only flow in the subsurface
    subsurface = (surface_level - y) > surface_buffer
    qh = qh * subsurface
    qv = qv * subsurface

    return qh, qv


def build_boundary_flux_source(mesh, boundary_faces, flux):

    """
    Convert a specified heat flux over a set of boundary faces into a volume
    source term for the cells next to those faces, in W/m3.

    The flux can be a single value or an array with a value for each face of
    the mesh. An array is used for radially symmetric models, where the flux
    over each face is weighted by the radius of that face.
    """

    face_areas = np.array(mesh._faceAreas)
    owner_cells = np.array(mesh.faceCellIDs)[0]

    flux = np.asarray(flux, dtype=float)
    if flux.ndim == 0:
        face_flux = float(flux) * face_areas[boundary_faces]
    else:
        face_flux = flux[boundary_faces] * face_areas[boundary_faces]

    source = np.zeros(mesh.numberOfCells)
    np.add.at(source, owner_cells[boundary_faces].astype(int), face_flux)

    return source / np.array(mesh.cellVolumes)


def check_radial_parameters(mp):

    """
    Check that the model parameters are consistent with a radially symmetric
    model. The fault has to be a single vertical fault on the symmetry axis,
    anything else does not describe a body of rotation.
    """

    if len(mp.fault_xs) != 1:
        raise ValueError('error, a radially symmetric model can only contain a '
                         'single fault, %i faults were specified in the model '
                         'parameter file' % len(mp.fault_xs))

    if abs(mp.fault_xs[0]) > 1e-6:
        raise ValueError('error, the fault of a radially symmetric model has to '
                         'be located on the symmetry axis. set fault_xs = [0.0] '
                         'in the model parameter file, instead of [%0.2f]'
                         % mp.fault_xs[0])

    if abs(abs(mp.fault_angles[0]) - 90.0) > 1e-6:
        raise ValueError('error, the fault of a radially symmetric model has to '
                         'be vertical. set fault_angles = [90.0] in the model '
                         'parameter file, instead of [%0.2f]' % mp.fault_angles[0])

    for aquifer_angle in mp.aquifer_angles:
        if aquifer_angle is not None and abs(aquifer_angle) > 1e-6:
            raise ValueError('error, the aquifers of a radially symmetric model '
                             'have to be horizontal, but an aquifer angle of '
                             '%0.2f degrees was specified' % aquifer_angle)


def build_heat_flow_equation(mesh, K_face, rho_c, q_face, rho_f_c_f,
                             convection_term, source, penalty_coeff,
                             penalty_source, solve_as_steady_state, dt):

    """
    Assemble the heat advection and conduction equation.

    The advection term is written in non-conservative form, see the module
    docstring for details. Dirichlet boundary conditions inside the model
    domain are imposed with a penalty term.
    """

    equation = (rho_f_c_f * convection_term(coeff=q_face)
                == DiffusionTerm(coeff=K_face)
                + ImplicitSourceTerm(coeff=rho_f_c_f * q_face.divergence)
                + source
                - ImplicitSourceTerm(coeff=penalty_coeff)
                + penalty_source)

    if solve_as_steady_state is False:
        equation = TransientTerm(coeff=rho_c) + equation

    return equation


def run_model_fipy(mp):

    """
    Set up the mesh, model parameters and boundary conditions and run a single
    model experiment with the fipy backend.

    Returns a dictionary with the modeled temperature and advective flux
    fields, with temperatures interpolated to the mesh vertices so that the
    output matches the output of the escript backend.
    """

    if check_fipy_available() is False:
        raise ImportError('error, the fipy backend requires fipy and the gmsh '
                          'python module. install these with pip install fipy gmsh')

    dt_growth_factor = getattr(mp, 'dt_growth_factor', None)
    dt_max = getattr(mp, 'dt_max', None)
    if (dt_growth_factor is None) != (dt_max is None):
        msg = 'error, dt_growth_factor and dt_max should either both be set '
        msg += 'to grow the timestep, or both left at None to use a constant dt'
        raise ValueError(msg)

    ############################
    # construct the mesh
    ############################
    print('constructing mesh (note, this may take a while....)')

    z_surface = 0
    z_base = -mp.total_depth

    try:
        assert np.min(np.array(mp.target_zs)) > mp.z_fine
    except AssertionError:
        msg = 'error, the lowest value of the model parameter target_z should '
        msg += 'always be higher than z_fine, otherwise the grid algorithm '
        msg += 'does not work.\n'
        mp.z_fine = np.min(np.array(mp.target_zs)) - mp.cellsize_fine
        msg += 'lowering z_fine to %0.1f' % mp.z_fine
        print(msg)

    if mp.add_exhumation is False:
        print('construct static mesh without exhumation')
        mp.exhumation_rate = 0.0
        exhumed_thickness = 0
        elevation_top = z_surface + mp.air_height
    else:
        print('construct dynamic mesh with exhumation')
        exhumed_thickness = mp.exhumation_rate * (np.sum(np.array(mp.durations)) / mp.year)
        exhumation_steps = mp.exhumation_steps

        if exhumed_thickness / exhumation_steps < mp.min_layer_thickness:
            print('warning, exhumation levels would be smaller than %0.2f m'
                  % mp.min_layer_thickness)
            exhumation_steps = int(np.ceil(exhumed_thickness) / mp.min_layer_thickness)
            if exhumation_steps < 1:
                exhumation_steps = 1
            print('reducing exhumation steps to %i' % exhumation_steps)

        if exhumed_thickness != 0:
            mp.target_zs = np.linspace(0, exhumed_thickness, exhumation_steps + 1)

        elevation_top = z_surface + exhumed_thickness + mp.air_height

    radial = getattr(mp, 'radial_symmetry', False)

    if radial is True:
        check_radial_parameters(mp)
        mesh, x_flt, z_flt, z_base = mesh_functions_fipy.setup_radial_mesh_fipy(
            mp.width, mp.fault_widths[0], mp.fault_bottoms[0] - 500.0,
            elevation_top, mp.z_fine, z_base, mp.cellsize,
            mp.cellsize_air, mp.cellsize_surface,
            mp.cellsize_fault_surface, mp.cellsize_fault,
            mp.cellsize_fine, mp.cellsize_base, mp.target_zs,
            mesh_filename=getattr(mp, 'fipy_mesh_filename', 'beo_mesh_radial.msh'))
    else:
        mesh, x_flt, z_flt, z_base = mesh_functions_fipy.setup_mesh_with_exhumation_fipy(
            mp.width, mp.fault_xs[0], mp.fault_widths[0],
            mp.fault_angles[0], mp.fault_bottoms[0] - 500.0,
            elevation_top, mp.z_fine, z_base, mp.cellsize,
            mp.cellsize_air, mp.cellsize_surface,
            mp.cellsize_fault_surface, mp.cellsize_fault,
            mp.cellsize_fine, mp.cellsize_base, mp.target_zs,
            mp.discretize_borehole, mp.borehole_xs, mp.borehole_depths,
            mp.borehole_cellsize,
            mesh_filename=getattr(mp, 'fipy_mesh_filename', 'beo_mesh_fipy.msh'))

    if mp.add_exhumation is False:
        x_flt = np.ones(len(mp.target_zs)) * mp.fault_xs[0]
        z_flt = np.zeros(len(mp.target_zs))

    #########################################
    # coordinates of cells, faces and vertices
    #########################################
    cell_centers = np.array(mesh.cellCenters)
    xc, yc = cell_centers[0], cell_centers[1]
    n_cells = mesh.numberOfCells

    face_centers = np.array(mesh.faceCenters)
    xf, yf = face_centers[0], face_centers[1]

    vertex_coords = np.array(mesh.vertexCoords)
    xyz_array = vertex_coords.T
    vertex_weights = build_cell_to_vertex_weights(mesh)

    exterior = np.array(mesh.exteriorFaces, dtype=bool)

    ###############################
    # boundary conditions
    ###############################
    print('setting up boundary conditions')

    y_tolerance = 1e-6
    top_faces = exterior & (yf > yf.max() - y_tolerance)
    bottom_faces = exterior & (yf < yf.min() + y_tolerance)

    if mp.thermal_gradient is not None:
        bottom_temperature = mp.air_temperature + mp.thermal_gradient * mp.total_depth
        basal_flux_source = np.zeros(n_cells)
        print('using a specified temperature of %0.1f degrees C at the base of '
              'the model domain' % bottom_temperature)
    else:
        bottom_temperature = None
        # for a radially symmetric model the heat flow equation is multiplied
        # by the radius, so the flux over each boundary face is weighted by the
        # radius of that face
        if radial is True:
            basal_heat_flux = mp.basal_heat_flux * xf
        else:
            basal_heat_flux = mp.basal_heat_flux
        basal_flux_source = build_boundary_flux_source(mesh, bottom_faces,
                                                       basal_heat_flux)
        print('using a specified heat flux of %0.3f W/m2 at the base of the '
              'model domain' % mp.basal_heat_flux)

    ###############################
    # material properties
    ###############################
    K_b, rho_b, c_b = build_material_properties(xc, yc, mp, mp.fault_angles[0],
                                                mp.fault_xs[0])

    surface_level = exhumed_thickness
    subsurface = yc <= surface_level

    if mp.variable_K_air is True:
        # base K_air on a small difference between air and land temperature,
        # otherwise the equation for the heat transfer coefficient fails
        import beo_core
        K_air = beo_core.calculate_surface_heat_transfer_coeff(
            mp.rho_air, mp.c_air, mp.ra, mp.dz, mp.air_temperature + 0.5,
            mp.air_temperature)
    else:
        K_air = mp.K_air

    K_var = np.where(subsurface, K_b, K_air)
    c_var = np.where(subsurface, c_b, mp.c_air)
    rho_var = np.where(subsurface, rho_b, mp.rho_air)

    ###############################
    # faults and aquifers
    ###############################
    print('locating fault zones and aquifers')

    fault_zones_c, fault_segments_c, aquifer_locs_c = \
        build_fault_and_aquifer_masks(xc, yc, mp, surface_level)
    fault_zones_f, fault_segments_f, aquifer_locs_f = \
        build_fault_and_aquifer_masks(xf, yf, mp, surface_level)

    for i, fault_zone in enumerate(fault_zones_c):
        if not np.any(fault_zone):
            print('warning, fault zone %i contains no cells, check the fault '
                  'parameters and the cell size in the input file' % i)

    if radial is True:
        # for a radially symmetric model the fault is a cylindrical conduit on
        # the symmetry axis with a radius of half the fault width. the fault
        # fluxes are specified as a total flow in m3/sec and are spread evenly
        # over the cross section of the conduit
        fault_fluxes_m_per_sec = []
        for fault_flux_timeslice in mp.fault_fluxes:
            fts_timeslice = []
            for fault_width, fault_flux in zip(mp.fault_widths, fault_flux_timeslice):
                conduit_area = np.pi * (fault_width / 2.0)**2
                fts_timeslice.append([f / conduit_area for f in fault_flux])
            fault_fluxes_m_per_sec.append(fts_timeslice)

        # the aquifer fluxes are also specified as a total flow in m3/sec. the
        # flow crosses a cylinder with an area of 2 pi r times the aquifer
        # thickness, so here only the division by the thickness is done and the
        # division by the circumference is done in build_flux_field, where the
        # radius of each cell and face is known
        aquifer_fluxes_m_per_sec = []
        for aquifer_flux_timeslice in mp.aquifer_fluxes:
            aquifer_flux_ii = []
            for aquifer_top, aquifer_bottom, aquifer_flux in \
                    zip(mp.aquifer_tops, mp.aquifer_bottoms, aquifer_flux_timeslice):
                if aquifer_top is not None:
                    aquifer_flux_ii.append(aquifer_flux / (aquifer_top - aquifer_bottom))
            aquifer_fluxes_m_per_sec.append(aquifer_flux_ii)
    else:
        # convert fault fluxes from m2/sec to m/sec by dividing by the fault width
        fault_fluxes_m_per_sec = []
        for fault_flux_timeslice in mp.fault_fluxes:
            fts_timeslice = []
            for fault_width, fault_flux in zip(mp.fault_widths, fault_flux_timeslice):
                fts_timeslice.append([f / fault_width for f in fault_flux])
            fault_fluxes_m_per_sec.append(fts_timeslice)

        # convert aquifer fluxes from m2/sec to m/sec
        aquifer_fluxes_m_per_sec = []
        for aquifer_flux_timeslice in mp.aquifer_fluxes:
            aquifer_flux_ii = []
            for aquifer_top, aquifer_bottom, aquifer_flux in \
                    zip(mp.aquifer_tops, mp.aquifer_bottoms, aquifer_flux_timeslice):
                if aquifer_top is not None:
                    aquifer_flux_ii.append(aquifer_flux / (aquifer_top - aquifer_bottom))
            aquifer_fluxes_m_per_sec.append(aquifer_flux_ii)

    ###############################
    # run the heat flow model
    ###############################
    output_time_selector = OutputTimeSelector(
        mp.dt, mp.dt_stored, getattr(mp, 'stored_timesteps', None))
    if dt_growth_factor is not None:
        print('the timestep dt starts at %0.3e s and grows by a factor %0.3f '
              'per timestep up to a maximum of %0.3e s' % (mp.dt, dt_growth_factor, dt_max))
        print('the exact times at which results are stored depend on this '
              'growing timestep, so are not reported in advance')
    else:
        output_time_selector.report_actual_times(mp.dt)

    convection_scheme = getattr(mp, 'fipy_convection_scheme', 'powerlaw')
    convection_term = get_convection_term(convection_scheme)
    print('using the %s scheme for the advection term' % convection_scheme)

    results = model_hydrothermal_temperatures_fipy(
        mesh, mp, xc, yc, xf, yf,
        fault_zones_c, fault_segments_c, aquifer_locs_c,
        fault_zones_f, fault_segments_f, aquifer_locs_f,
        fault_fluxes_m_per_sec, aquifer_fluxes_m_per_sec,
        K_b, rho_b, c_b, K_air, K_var, rho_var, c_var,
        top_faces, bottom_faces, bottom_temperature, basal_flux_source,
        surface_level, vertex_weights, convection_term,
        output_time_selector)

    (runtimes, T_steady, Ts_cell, qh_cell, qv_cell, surface_levels,
     surface_levels_mesh, boiling_temps, exceed_boiling_temps) = results

    ##########################################
    # interpolate the results to the vertices
    ##########################################
    print('interpolating the model results to the mesh vertices')

    T_list = []
    for T_cell, surface_level_i in zip(Ts_cell, surface_levels_mesh):
        T_list.append(_cell_field_to_vertices(T_cell, yc, vertex_coords[1],
                                              surface_level_i, vertex_weights))

    T_array = np.array(T_list)
    T_init_array = _cell_field_to_vertices(T_steady, yc, vertex_coords[1],
                                           surface_levels_mesh[0], vertex_weights)

    Ts = list(T_array)

    xyz_element_array = cell_centers.T
    qh_array = np.array(qh_cell)
    qv_array = np.array(qv_cell)
    q_vectors = [np.vstack((qh, qv)) for qh, qv in zip(qh_array, qv_array)]

    if boiling_temps is not None:
        boiling_temp_array = np.array(
            [_cell_field_to_vertices(bt, yc, vertex_coords[1], sl, vertex_weights)
             for bt, sl in zip(boiling_temps, surface_levels_mesh)])
        exceed_boiling_temp_array = np.array(
            [_cell_field_to_vertices(bt, yc, vertex_coords[1], sl, vertex_weights)
             for bt, sl in zip(exceed_boiling_temps, surface_levels_mesh)])
        xyz_array_exc = xyz_array
    else:
        boiling_temp_array = None
        exceed_boiling_temp_array = None
        xyz_array_exc = None

    return {'runtimes': np.array(runtimes),
            'xyz_array': xyz_array,
            'surface_levels': np.array(surface_levels),
            'x_flt': x_flt,
            'z_flt': z_flt,
            'Ts': Ts,
            'q_vectors': q_vectors,
            'T_init_array': T_init_array,
            'T_array': T_array,
            'boiling_temp_array': boiling_temp_array,
            'xyz_array_exc': xyz_array_exc,
            'exceed_boiling_temp_array': exceed_boiling_temp_array,
            'xyz_element_array': xyz_element_array,
            'qh_array': qh_array,
            'qv_array': qv_array,
            'fault_fluxes_m_per_sec': fault_fluxes_m_per_sec}


def _cell_field_to_vertices(values, y_cells, y_vertices, surface_level,
                            vertex_weights):

    """
    Interpolate a cell centered field to the mesh vertices.

    Cells in the air layer and cells in the subsurface are kept separate, so
    that the strong contrast in thermal conductivity at the land surface does
    not smear out the temperature at the land surface. Vertices at the land
    surface take the value of the subsurface cells.
    """

    subsurface_cells = y_cells <= surface_level
    air_cells = ~subsurface_cells

    values_subsurface = interpolate_to_vertices(values, vertex_weights,
                                                cell_mask=subsurface_cells)

    if not np.any(air_cells):
        return values_subsurface

    values_air = interpolate_to_vertices(values, vertex_weights,
                                         cell_mask=air_cells)

    air_vertices = y_vertices > surface_level + 1e-6

    return np.where(air_vertices, values_air, values_subsurface)


def _surface_temperature(T_cell, y_cells, x_vertices, y_vertices, x_target,
                         surface_level, vertex_weights):

    """
    Find the temperature at the land surface as a function of the x coordinate
    and interpolate this to the coordinates in x_target.
    """

    T_vertex = interpolate_to_vertices(T_cell, vertex_weights,
                                       cell_mask=(y_cells <= surface_level))

    ind = np.abs(y_vertices - surface_level) < 0.01

    # the land surface does not always coincide with a horizontal line of nodes
    # in the mesh, which is the case if track_exact_surface_elev is used. fall
    # back to the closest line of nodes at or below the land surface
    if np.sum(ind) < 2:
        levels = np.unique(y_vertices[y_vertices <= surface_level + 0.01])

        for closest_level in levels[::-1]:
            ind = np.abs(y_vertices - closest_level) < 0.01
            if np.sum(ind) >= 2:
                break

    if np.sum(ind) < 2:
        print('warning, could not find land surface nodes at an elevation of '
              '%0.2f m, using the mean subsurface temperature instead'
              % surface_level)
        subsurface = y_cells <= surface_level
        if not np.any(subsurface):
            subsurface = np.ones_like(y_cells, dtype=bool)
        return np.ones_like(x_target) * np.mean(T_cell[subsurface])

    x_surface = x_vertices[ind]
    T_surface = T_vertex[ind]

    a = np.argsort(x_surface)

    return np.interp(x_target, x_surface[a], T_surface[a])


def model_hydrothermal_temperatures_fipy(mesh, mp, xc, yc, xf, yf,
                                         fault_zones_c, fault_segments_c, aquifer_locs_c,
                                         fault_zones_f, fault_segments_f, aquifer_locs_f,
                                         fault_fluxes, aquifer_fluxes,
                                         K_b, rho_b, c_b, K_air, K_var, rho_var, c_var,
                                         top_faces, bottom_faces, bottom_temperature,
                                         basal_flux_source, surface_level_init,
                                         vertex_weights, convection_term,
                                         output_time_selector,
                                         surface_buffer=0.1,
                                         max_screen_output=500):

    """
    Full single model run of the heat conduction and advection model with the
    fipy backend.

    This is the fipy equivalent of model_hydrothermal_temperatures in
    beo_core.py. All fields are calculated at the cell centers, apart from the
    advective flux, which is calculated at the cell faces as well.
    """

    import beo_core

    fluid_density = 1000.0
    g = 9.81
    atmospheric_pressure = 101325

    c1 = 3.866
    c2 = 25.151
    c3 = 103.28
    c4 = 179.99

    P_buffer = 10.0

    solve_as_steady_state = mp.steady_state
    dt = mp.dt
    dt_initial = mp.dt
    dt_growth_factor = getattr(mp, 'dt_growth_factor', None)
    dt_max = getattr(mp, 'dt_max', None)
    vapour_correction = mp.vapour_correction
    variable_K_air = mp.variable_K_air
    exhumation_rate = mp.exhumation_rate
    target_depths = np.array(mp.target_zs)

    n_cells = mesh.numberOfCells
    cell_volumes = np.array(mesh.cellVolumes)
    cell_size = np.sqrt(cell_volumes)

    vertex_coords = np.array(mesh.vertexCoords)
    x_vertices, y_vertices = vertex_coords[0], vertex_coords[1]

    rho_f_c_f = mp.rho_f * mp.c_f

    ###################################
    # set up the fipy variables and PDE
    ###################################
    T = CellVariable(name='temperature', mesh=mesh, value=mp.air_temperature,
                     hasOld=True)

    # for a radially symmetric model the heat flow equation is multiplied by
    # the radial coordinate, which turns the axisymmetric equation into an
    # ordinary two dimensional equation with all coefficients weighted by the
    # radius. see the module docstring of mesh_functions_fipy for the geometry
    radial = getattr(mp, 'radial_symmetry', False)

    if radial is True:
        print('radially symmetric model, weighting the coefficients of the heat '
              'flow equation with the radius')
        # the faces that lie on the symmetry axis have a radius of exactly
        # zero, which would make the weighted thermal conductivity zero as
        # well. the convection schemes of fipy divide by the diffusion
        # coefficient to find the grid peclet number, so a radius of zero
        # results in a division by zero and a matrix filled with nan. no heat
        # crosses the symmetry axis anyway, so the radius of these faces is
        # replaced by a small value that is negligible compared to the size of
        # the cells next to the axis
        radius_floor = 1e-6 * np.max(xf)
        radial_weight_cell = np.maximum(xc, radius_floor)
        radial_weight_face = np.maximum(xf, radius_floor)
    else:
        radial_weight_cell = np.ones(n_cells)
        radial_weight_face = np.ones(mesh.numberOfFaces)

    K_cell = CellVariable(mesh=mesh, value=K_var)
    rho_c_cell = CellVariable(mesh=mesh, value=radial_weight_cell * rho_var * c_var)
    q_face = FaceVariable(mesh=mesh, rank=1, value=0.0)
    source = CellVariable(mesh=mesh, value=basal_flux_source)
    penalty_coeff = CellVariable(mesh=mesh, value=0.0)
    penalty_source = CellVariable(mesh=mesh, value=0.0)

    # use the harmonic mean of the thermal conductivity of the cells on either
    # side of a face, an arithmetic mean would smooth out the strong contrast
    # in thermal conductivity at the land surface
    K_face = FaceVariable(mesh=mesh, value=radial_weight_face) \
        * K_cell.harmonicFaceValue

    T.constrain(mp.air_temperature, where=FaceVariable(mesh=mesh, value=top_faces))

    if bottom_temperature is not None:
        T.constrain(bottom_temperature,
                    where=FaceVariable(mesh=mesh, value=bottom_faces))

    # scale the penalty term used for specified temperatures inside the model
    # domain to the size of the diagonal of the coefficient matrix
    penalty_scale = 1e6 * np.max(K_var / cell_volumes)
    if solve_as_steady_state is False:
        penalty_scale = max(penalty_scale, 1e6 * np.max(rho_var * c_var / dt))

    ######################################
    # model steady-state temperature field
    ######################################
    print('finding solution for steady state HF')

    equation_steady = (DiffusionTerm(coeff=K_face) + source
                       - ImplicitSourceTerm(coeff=penalty_coeff)
                       + penalty_source == 0)
    equation_steady.solve(var=T, solver=get_solver())

    T_steady = np.array(T).copy()

    print('modeled steady state temperatures, min = %0.1f, max = %0.1f degrees C'
          % (T_steady.min(), T_steady.max()))

    equation = build_heat_flow_equation(mesh, K_face, rho_c_cell, q_face,
                                        rho_f_c_f, convection_term, source,
                                        penalty_coeff, penalty_source,
                                        solve_as_steady_state, dt)

    # reuse the LU factorization of the heat flow equation as long as the
    # coefficients of the equation stay the same
    if solve_as_steady_state is False and getattr(mp, 'cache_lu', True) is True:
        lu_cache = TransientLUCache(
            validate=getattr(mp, 'validate_lu_cache', False))
        print('reusing the LU factorization of the heat flow equation')
    else:
        lu_cache = None

    print('starting transient heat flow calculations')

    runtimes = [0.0]
    t_total = 0
    Ts = [T_steady]
    qh_stored = [np.zeros(n_cells)]
    qv_stored = [np.zeros(n_cells)]
    surface_level = surface_level_init
    surface_level_mesh = surface_level_init
    surface_levels = [surface_level]
    surface_levels_mesh = [surface_level_mesh]

    exceed_max_liquid_T_old = None

    depth = -(yc - surface_level)
    subsurface = yc <= surface_level

    if vapour_correction is True:
        P = np.where(depth > 0, fluid_density * g * depth + atmospheric_pressure,
                     atmospheric_pressure)
        logP = np.log10(P / 1.0e6)
        boiling_temp = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4
        exceed_boiling_temp = (subsurface & (T_steady > boiling_temp)).astype(float)

        boiling_temps = [boiling_temp]
        exceed_boiling_temps = [exceed_boiling_temp]

        if np.any(P < calculate_vapour_pressure(T_steady)):
            print('warning, vapour present at initial steady-state P-T conditions')
    else:
        boiling_temps = None
        exceed_boiling_temps = None

    try:
        assert len(fault_fluxes) == len(mp.durations)
    except AssertionError:
        msg = 'error, the length of the fault_fluxes and durations arrays or list do not match'
        msg += '\ncheck you model parameters file'
        raise ValueError(msg)

    n_time_periods = len(mp.durations)

    for time_period, fault_flux, duration in zip(itertools.count(), fault_fluxes,
                                                 mp.durations):

        print('time period %i of %i' % (time_period + 1, n_time_periods))
        print('duration of %i' % duration)
        print('fluxes in faults: ', fault_flux)

        if time_period < len(aquifer_fluxes):
            aquifer_flux = aquifer_fluxes[time_period]
        else:
            aquifer_flux = []

        # set up the advective flow field in the faults and aquifers, at the
        # cell faces for the PDE and at the cell centers for the model output
        qh_f, qv_f = build_flux_field(xf, yf, fault_zones_f, fault_segments_f,
                                      mp.fault_angles, fault_flux, aquifer_locs_f,
                                      aquifer_flux, mp.aquifer_angles,
                                      surface_level, surface_buffer, radial)
        qh_c, qv_c = build_flux_field(xc, yc, fault_zones_c, fault_segments_c,
                                      mp.fault_angles, fault_flux, aquifer_locs_c,
                                      aquifer_flux, mp.aquifer_angles,
                                      surface_level, surface_buffer, radial)

        print('modeled advective flux:')
        print('\tqh from %0.2e to %0.2e m/sec' % (qh_c.min(), qh_c.max()))
        print('\tqz from %0.2e to %0.2e m/sec' % (qv_c.min(), qv_c.max()))

        # the advective flux is set to zero at the boundary of the model
        # domain, heat is only advected inside the model domain
        exterior = np.array(mesh.exteriorFaces, dtype=bool)
        qh_f[exterior] = 0.0
        qv_f[exterior] = 0.0

        q_face.setValue(np.vstack((qh_f * radial_weight_face,
                                   qv_f * radial_weight_face)))

        if solve_as_steady_state is False and dt_growth_factor is not None:
            # the timestep ramps up again at the start of every time period,
            # so that a period with an abrupt change in fault flux also gets
            # a well resolved start, and not just the very first period
            dt_schedule = build_growing_dt_schedule(dt_initial, dt_growth_factor, dt_max, duration)
            nt = len(dt_schedule)
            print('timestep grows from %0.3e s to a maximum of %0.3e s over %i timesteps'
                  % (dt_schedule[0], max(dt_schedule), nt))
        else:
            dt_schedule = None
            nt = int(duration / dt)

        # calculate the courant number and the grid peclet number. if the
        # timestep is growing, use the largest timestep actually reached in
        # this time period, since the courant number and peclet number
        # increase with dt
        dt_cfl = max(dt_schedule) if dt_schedule is not None else dt
        q_abs = np.sqrt(qh_c**2 + qv_c**2)
        cfl_cond_max = np.max(q_abs * dt_cfl / cell_size)
        print('max. CFL number: %0.2f' % cfl_cond_max)
        if cfl_cond_max > 1:
            print('warning, cfl condition not met, cfl number = %0.2f' % cfl_cond_max)

        Pe = rho_f_c_f * q_abs * cell_size / K_var
        print('max. grid peclet number = %0.2f' % np.max(Pe))

        if solve_as_steady_state is True:
            if variable_K_air is True or vapour_correction is True:
                nt = mp.n_iterations_steady_state
                print('running steady-state model with %i iterations for '
                      'variable K air or vapour correction' % nt)
            else:
                nt = 1
                print('running steady-state model without iterations')

        screen_output_interval = int(np.ceil(nt / max_screen_output))
        if screen_output_interval < 1:
            screen_output_interval = 1
        print('screen output once every %i timesteps' % screen_output_interval)

        print('x' * 10)
        print('starting iterations, total timesteps = %i' % nt)
        start = time.time()

        for t in range(nt):

            if dt_schedule is not None:
                dt = dt_schedule[t]

            # exhumation, calculate the new surface level and snap this to the
            # closest horizontal line of nodes in the mesh
            if exhumation_rate > 0:
                surface_level = surface_level_init - t_total / year * exhumation_rate

                try:
                    surface_level_mesh_id = np.where(target_depths <= surface_level)[0][-1]
                    surface_level_mesh = target_depths[surface_level_mesh_id]
                except IndexError:
                    surface_level_mesh = surface_level

            if exhumation_rate != 0 and solve_as_steady_state is False:

                if getattr(mp, 'track_exact_surface_elev', False) is True:
                    surface_level_active = surface_level
                else:
                    surface_level_active = surface_level_mesh

                subsurface = yc <= surface_level_active

                # switch off the advective flux in the exhumed part of the
                # model domain
                subsurface_f = (surface_level_active - yf) > surface_buffer
                q_face.setValue(np.vstack(
                    (qh_f * subsurface_f * radial_weight_face,
                     qv_f * subsurface_f * radial_weight_face)))
            else:
                surface_level_active = surface_level_mesh

            # calculate the effective thermal conductivity of the air layer
            # from the sensible and latent heat flux at the land surface
            if variable_K_air is True:
                surface_T_int = _surface_temperature(np.array(T), yc, x_vertices,
                                                     y_vertices, xc,
                                                     surface_level_active,
                                                     vertex_weights)
                K_air = beo_core.calculate_surface_heat_transfer_coeff(
                    mp.rho_air, mp.c_air, mp.ra, mp.dz, surface_T_int,
                    mp.air_temperature)

            if variable_K_air is True or exhumation_rate != 0:
                K_var = np.where(subsurface, K_b, K_air)
                c_var = np.where(subsurface, c_b, mp.c_air)
                rho_var = np.where(subsurface, rho_b, mp.rho_air)

                K_cell.setValue(K_var)
                rho_c_cell.setValue(radial_weight_cell * rho_var * c_var)

            # recalculate the fluid pressure and the boiling temperature
            if vapour_correction is True:
                depth = -(yc - surface_level_active)
                P = np.where(depth > 0,
                             fluid_density * g * depth + atmospheric_pressure,
                             atmospheric_pressure)
                logP = np.log10(P / 1.0e6)
                boiling_temp = c1 * logP**3 + c2 * logP**2 + c3 * logP + c4

                T_now = np.array(T)
                if exceed_max_liquid_T_old is None:
                    exceed_boiling_temp = (subsurface & (T_now > boiling_temp))
                else:
                    exceed_boiling_temp = (subsurface & (exceed_max_liquid_T_old == 0)
                                           & (T_now > boiling_temp))

                exceed_boiling_temp = exceed_boiling_temp.astype(float)
                exceed_max_liquid_T_old = exceed_boiling_temp
            else:
                boiling_temp = 1e6
                exceed_boiling_temp = np.zeros(n_cells)

            # specified temperatures inside the model domain, these are the
            # exhumed part of the model domain, which is set to the air
            # temperature, and cells that exceed the boiling temperature
            specified_T_loc = np.zeros(n_cells)
            specified_T = np.zeros(n_cells)

            # the escript backend uses the exact rather than the mesh snapped
            # surface level here, keep the same convention
            exhumed = yc >= surface_level + mp.air_height
            specified_T_loc[exhumed] = 1.0
            specified_T[exhumed] = mp.air_temperature

            if vapour_correction is True:
                boiling = exceed_boiling_temp > 0
                specified_T_loc[boiling] = 1.0
                specified_T[boiling] = boiling_temp[boiling]

            penalty_coeff.setValue(penalty_scale * specified_T_loc)
            penalty_source.setValue(penalty_scale * specified_T_loc * specified_T)

            # solve the heat flow equation for the temperature
            solve_time_start = time.time()
            if solve_as_steady_state is False:
                T.updateOld()
                if lu_cache is not None:
                    lu_cache.solve(equation, T, dt,
                                   [K_cell, rho_c_cell, q_face, penalty_coeff],
                                   [source, penalty_source])
                else:
                    equation.solve(var=T, dt=dt, solver=get_solver())
            else:
                equation.solve(var=T, solver=get_solver())
            solver_time = time.time() - solve_time_start

            t_total += dt
            if solve_as_steady_state is True:
                t_total = 0.0

            # output to screen
            if t % screen_output_interval == 0 or t == nt - 1:

                end = time.time()
                comptime = end - start
                start = end

                print('step %i of %i' % (t + 1, nt))
                print('\truntime = %0.2e yrs' % (t_total / year))
                print('\tcomputational time for one timestep = %0.3f sec'
                      % (comptime / screen_output_interval))
                print('\tcomputational time needed by solver = %0.3f sec' % solver_time)
                print('\tactual surface level ', surface_level)
                if exhumation_rate > 0:
                    print('\tclosest surface in mesh ', surface_level_mesh)
                T_now = np.array(T)
                print('\ttemperature from %0.1f to %0.1f degrees C'
                      % (T_now.min(), T_now.max()))

                if vapour_correction is True:
                    vapour = subsurface & (P - calculate_vapour_pressure(T_now)
                                           + P_buffer < 0)
                    if np.any(vapour):
                        print('\tvapour present in: %0.1f m^2'
                              % np.sum(cell_volumes[vapour]))
                        print('\t\tfrom x = %0.1f to x = %0.1f'
                              % (xc[vapour].min(), xc[vapour].max()))
                        print('\t\tand from y = %0.1f to y = %0.1f'
                              % (yc[vapour].min(), yc[vapour].max()))
                    else:
                        print('\tno vapour present')

                if variable_K_air is True:
                    print('\tcalculated K air from %0.3e to %0.3e W/(m K)'
                          % (np.min(K_air), np.max(K_air)))

            # store output
            if output_time_selector.should_store(t_total, dt):
                Ts.append(np.array(T).copy())
                qh_stored.append(qh_c.copy())
                qv_stored.append(qv_c.copy())
                surface_levels.append(surface_level)
                surface_levels_mesh.append(surface_level_active)

                if vapour_correction is True:
                    boiling_temps.append(boiling_temp.copy())
                    exceed_boiling_temps.append(exceed_boiling_temp.copy())

                runtimes.append(t_total)

        print('final T after advective heating, min = %0.1f, max = %0.1f degrees C'
              % (np.array(T).min(), np.array(T).max()))

        if lu_cache is not None:
            print('LU cache: %i factorizations, %i reused solves'
                  % (lu_cache.n_factorizations, lu_cache.n_reused))
            if lu_cache.validate is True:
                print('LU cache: max relative error in reconstructed rhs = %0.3e'
                      % lu_cache.max_rhs_error)

    return (runtimes, T_steady, Ts, qh_stored, qv_stored, surface_levels,
            surface_levels_mesh, boiling_temps, exceed_boiling_temps)
