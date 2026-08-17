"""
Tests of the fipy backend of Beo.

The tests compare the numerical solution of the heat flow equation with
analytical solutions for conductive heat flow and for advective and
conductive heat flow.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

fipy = pytest.importorskip('fipy')

from fipy import Grid2D, CellVariable, FaceVariable, DiffusionTerm
from fipy import PowerLawConvectionTerm, ImplicitSourceTerm

import lib.beo_fipy as beo_fipy


def build_test_mesh(nx=4, ny=100, dx=10.0, dy=10.0):

    """
    Simple rectangular mesh for a one dimensional geotherm
    """

    return Grid2D(nx=nx, ny=ny, dx=dx, dy=dy)


def test_steady_state_geotherm():

    """
    Conductive heat flow with a specified temperature at the top and a
    specified heat flux at the base should reproduce a linear geotherm
    """

    ny = 100
    dy = 10.0
    mesh = build_test_mesh(ny=ny, dy=dy)

    K = 2.5
    basal_heat_flux = 0.07
    surface_temperature = 10.0

    yf = np.array(mesh.faceCenters)[1]
    exterior = np.array(mesh.exteriorFaces, dtype=bool)
    top_faces = exterior & (yf > yf.max() - 1e-6)
    bottom_faces = exterior & (yf < yf.min() + 1e-6)

    source = beo_fipy.build_boundary_flux_source(mesh, bottom_faces, basal_heat_flux)

    T = CellVariable(mesh=mesh, value=surface_temperature)
    K_cell = CellVariable(mesh=mesh, value=K)
    T.constrain(surface_temperature, where=FaceVariable(mesh=mesh, value=top_faces))

    equation = (DiffusionTerm(coeff=K_cell.harmonicFaceValue)
                + CellVariable(mesh=mesh, value=source) == 0)
    equation.solve(var=T)

    yc = np.array(mesh.cellCenters)[1]
    depth = yc.max() + dy / 2.0 - yc
    T_analytical = surface_temperature + basal_heat_flux / K * depth

    assert np.allclose(np.array(T), T_analytical, atol=1e-6)


def test_advection_conduction_profile():

    """
    Steady state heat flow with a uniform upward fluid flux should reproduce
    the analytical solution for one dimensional advective and conductive heat
    flow, see for example Bredehoeft and Papadopulos (1965)
    """

    ny = 200
    dy = 5.0
    mesh = build_test_mesh(ny=ny, dy=dy)

    K = 2.5
    rho_f_c_f = 1000.0 * 4180.0
    q = 1e-9
    T_top = 10.0
    T_bottom = 100.0
    L = ny * dy

    yf = np.array(mesh.faceCenters)[1]
    exterior = np.array(mesh.exteriorFaces, dtype=bool)
    top_faces = exterior & (yf > yf.max() - 1e-6)
    bottom_faces = exterior & (yf < yf.min() + 1e-6)

    T = CellVariable(mesh=mesh, value=T_top)
    K_cell = CellVariable(mesh=mesh, value=K)
    T.constrain(T_top, where=FaceVariable(mesh=mesh, value=top_faces))
    T.constrain(T_bottom, where=FaceVariable(mesh=mesh, value=bottom_faces))

    # uniform upward flux, set to zero at the model boundary so that heat is
    # only advected inside the model domain
    q_face = FaceVariable(mesh=mesh, rank=1, value=0.0)
    q_values = np.zeros((2, mesh.numberOfFaces))
    q_values[1] = q
    q_values[:, exterior] = 0.0
    q_face.setValue(q_values)

    equation = (rho_f_c_f * PowerLawConvectionTerm(coeff=q_face)
                == DiffusionTerm(coeff=K_cell.harmonicFaceValue)
                + ImplicitSourceTerm(coeff=rho_f_c_f * q_face.divergence))
    equation.solve(var=T)

    yc = np.array(mesh.cellCenters)[1]
    depth = yc.max() + dy / 2.0 - yc

    # analytical solution, with the peclet number based on the model thickness
    beta = rho_f_c_f * q * L / K
    T_analytical = T_top + (T_bottom - T_top) \
        * (1.0 - np.exp(-beta * depth / L)) / (1.0 - np.exp(-beta))

    assert np.max(np.abs(np.array(T) - T_analytical)) < 0.5


def test_non_conservative_advection_term():

    """
    A divergent flux field, as used in Beo where the fluid flux is switched on
    abruptly at the fault walls, should not change a uniform temperature
    field. This is only the case if the advection term is written in
    non-conservative form.
    """

    mesh = build_test_mesh(nx=40, ny=40, dx=10.0, dy=10.0)

    T_init = 50.0
    rho_f_c_f = 1000.0 * 4180.0

    xf = np.array(mesh.faceCenters)[0]
    exterior = np.array(mesh.exteriorFaces, dtype=bool)

    q_values = np.zeros((2, mesh.numberOfFaces))
    fault = (xf > 150.0) & (xf < 250.0)
    q_values[0][fault] = 1e-6
    q_values[1][fault] = 1e-6
    q_values[:, exterior] = 0.0

    q_face = FaceVariable(mesh=mesh, rank=1, value=0.0)
    q_face.setValue(q_values)

    assert np.max(np.abs(np.array(q_face.divergence))) > 0

    K_cell = CellVariable(mesh=mesh, value=2.5)

    T = CellVariable(mesh=mesh, value=T_init, hasOld=True)
    equation = beo_fipy.build_heat_flow_equation(
        mesh, K_cell.harmonicFaceValue,
        CellVariable(mesh=mesh, value=2500.0 * 900.0), q_face, rho_f_c_f,
        PowerLawConvectionTerm, CellVariable(mesh=mesh, value=0.0),
        CellVariable(mesh=mesh, value=0.0), CellVariable(mesh=mesh, value=0.0),
        False, 1e9)

    T.updateOld()
    equation.solve(var=T, dt=1e9)

    assert np.max(np.abs(np.array(T) - T_init)) < 1e-8


def test_cell_to_vertex_interpolation():

    """
    Interpolation of a linear field from the cell centers to the mesh vertices
    should reproduce the analytical values at all vertices, including the
    vertices at the boundary of the model domain, where all surrounding cell
    centers lie on one side of the vertex
    """

    mesh = build_test_mesh(nx=20, ny=20, dx=10.0, dy=10.0)

    cell_centers = np.array(mesh.cellCenters)
    vertex_coords = np.array(mesh.vertexCoords)

    values = 2.0 * cell_centers[0] + 3.0 * cell_centers[1] + 5.0

    weights = beo_fipy.build_cell_to_vertex_weights(mesh)
    values_vertex = beo_fipy.interpolate_to_vertices(values, weights)

    expected = 2.0 * vertex_coords[0] + 3.0 * vertex_coords[1] + 5.0

    assert np.max(np.abs(values_vertex - expected)) < 1e-6

    # interpolation using only the cells in the lower half of the mesh should
    # still be exact, this is used to keep the air layer and the subsurface
    # separate at the land surface
    lower = cell_centers[1] < 100.0
    values_lower = beo_fipy.interpolate_to_vertices(values, weights,
                                                    cell_mask=lower)

    at_or_below = vertex_coords[1] <= 100.0

    assert np.max(np.abs(values_lower[at_or_below] - expected[at_or_below])) < 1e-6


def test_radial_conduction_log_profile():

    """
    Radial conduction between two cylinders should follow a logarithmic
    temperature profile, which is very different from the linear profile of a
    cartesian model. This checks that weighting the coefficients of the heat
    flow equation with the radius really gives an axisymmetric model.
    """

    r1, r2 = 1.0, 100.0
    T1, T2 = 100.0, 10.0
    K = 2.5
    nr = 400

    dr = (r2 - r1) / nr
    mesh = Grid2D(dx=dr, nx=nr, dy=dr, ny=4) + ((r1,), (0.0,))

    rc = np.array(mesh.cellCenters)[0]
    rf = np.array(mesh.faceCenters)[0]
    exterior = np.array(mesh.exteriorFaces, dtype=bool)

    T = CellVariable(mesh=mesh, value=T1)
    T.constrain(T1, where=FaceVariable(mesh=mesh, value=exterior & (rf < r1 + 1e-9)))
    T.constrain(T2, where=FaceVariable(mesh=mesh, value=exterior & (rf > r2 - 1e-9)))

    r_face = FaceVariable(mesh=mesh, value=rf)
    (DiffusionTerm(coeff=r_face * K) == 0).solve(var=T, solver=beo_fipy.get_solver())

    T_analytical = T1 + (T2 - T1) * np.log(rc / r1) / np.log(r2 / r1)

    assert np.max(np.abs(np.array(T) - T_analytical)) < 0.2

    # the same model without the radial weighting is a cartesian model and
    # should give a clearly different, linear, temperature profile
    T_cartesian = CellVariable(mesh=mesh, value=T1)
    T_cartesian.constrain(T1, where=FaceVariable(mesh=mesh,
                                                 value=exterior & (rf < r1 + 1e-9)))
    T_cartesian.constrain(T2, where=FaceVariable(mesh=mesh,
                                                 value=exterior & (rf > r2 - 1e-9)))
    (DiffusionTerm(coeff=K) == 0).solve(var=T_cartesian,
                                        solver=beo_fipy.get_solver())

    assert np.max(np.abs(np.array(T_cartesian) - T_analytical)) > 10.0


def test_radial_solver_actually_solves():

    """
    Weighting the heat flow equation with the radius makes the right hand side
    of the equation much larger for the cells far away from the symmetry axis
    than for the cells close to it. Solvers that judge convergence by the norm
    of the right hand side then decide that the equation is already solved and
    return the temperatures of the previous timestep unchanged.

    This checks that the solver used by the fipy backend does not do that, by
    comparing the solution with a direct solve of the same matrix.
    """

    scipy_sparse = pytest.importorskip('scipy.sparse')
    scipy_linalg = pytest.importorskip('scipy.sparse.linalg')

    from fipy import TransientTerm

    # a wide range of radii, which is what causes the problem
    nr, nz = 100, 40
    mesh = Grid2D(dx=30.0, nx=nr, dy=50.0, ny=nz) + ((0.0,), (-2000.0,))

    rc = np.array(mesh.cellCenters)[0]
    zc = np.array(mesh.cellCenters)[1]
    rf = np.array(mesh.faceCenters)[0]
    zf = np.array(mesh.faceCenters)[1]
    exterior = np.array(mesh.exteriorFaces, dtype=bool)

    K, rho_c, rho_f_c_f = 2.0, 2.5e6, 4.18e6
    dt = 100.0 * 365.25 * 24 * 3600.0

    radius_floor = 1e-6 * np.max(rf)
    weight_cell = np.maximum(rc, radius_floor)
    weight_face = np.maximum(rf, radius_floor)

    T = CellVariable(mesh=mesh, value=10.0 - 0.035 * zc, hasOld=True)
    T.constrain(10.0, where=FaceVariable(mesh=mesh,
                                         value=exterior & (zf > -1e-9)))

    K_face = FaceVariable(mesh=mesh, value=weight_face * K)
    rho_c_cell = CellVariable(mesh=mesh, value=weight_cell * rho_c)

    q_values = np.zeros((2, mesh.numberOfFaces))
    conduit = (rf < 150.0) & (zf < -0.1)
    q_values[1][conduit] = 5e-7
    q_values[:, exterior] = 0.0
    q_face = FaceVariable(mesh=mesh, rank=1, value=0.0)
    q_face.setValue(q_values * weight_face)

    equation = (TransientTerm(coeff=rho_c_cell)
                + rho_f_c_f * PowerLawConvectionTerm(coeff=q_face)
                == DiffusionTerm(coeff=K_face)
                + ImplicitSourceTerm(coeff=rho_f_c_f * q_face.divergence))

    T_start = np.array(T).copy()
    T.updateOld()
    equation.cacheMatrix()
    equation.cacheRHSvector()
    equation.solve(var=T, dt=dt, solver=beo_fipy.get_solver())

    T_fipy = np.array(T).copy()

    # the temperature has to change, a solver that stops immediately returns
    # the temperatures of the previous timestep
    assert np.max(np.abs(T_fipy - T_start)) > 1.0

    # and the solution has to match a direct solve of the same matrix
    matrix = scipy_sparse.csr_matrix(equation.matrix.matrix)
    rhs = np.asarray(equation.RHSvector, dtype=float)
    T_direct = scipy_linalg.spsolve(matrix.tocsc(), rhs)

    assert np.max(np.abs(T_fipy - T_direct)) < 1e-6
