'''
grompylib

Elco Luijendijk, 2013
'''

import time, pdb, os

import numpy as np

try:
    import matplotlib.pyplot as pl
    import matplotlib.patches
    import matplotlib.collections

except:
    print 'warning, cannot import matplotlib'
    make_fig = False

try:
    import scipy.interpolate
except:
    print 'warning, cannot import scipy.interpolate'
    make_fig = False

import esys.escript as es
import esys.pycad as pc
import esys.pycad.gmsh as gmsh
import esys.finley as fl
import esys.escript.linearPDEs

# GIS
try:
    import osgeo.gdal as gdal
except:
    print 'warning: could not find gdal module, cannot load raster files'


def createRectangleMesh(xl, yl, d):
    
    '''
    use gmsh to construct rectangular mesh
    '''
    
    points = [pc.Point(0,0,0),pc.Point(xl,0,0),pc.Point(xl,yl,0),pc.Point(0,yl,0)]

    points.append(points[0])
    
    lines = [pc.Line(p1,p2) for p1,p2 in zip(points[:-1],points[1:])]
    
    curve = pc.CurveLoop(lines)

    surface = pc.PropertySet("surface",pc.PlaneSurface(curve))

    d = gmsh.Design(dim=2, element_size=d)
    
    # add tags
    d.addItems(surface,pc.PropertySet('top',lines[2]))
    
    d.setScriptFileName('rectangle.geo')
    d.setMeshFileName('rectangle.msh')
    
    mesh = fl.MakeDomain(d)

    mesh.write('rectangle.fly')

    return mesh


def solve_gwflow_eq( mypde, mesh, h, K, Q ):
    '''
    solve groundwater flow eq., stable density, incompressible.
    
    returns:
        h   hydraulic head
    '''
    
    
    A = -K * es.kronecker(mesh)
    #A = -K
    Y = -Q
    mypde.setValue( A=A, Y=Y )
    
    h = mypde.getSolution()
    
    return h


def convert_to_array(u):
    
    '''
    return the x,y coordinates and the value of escript variable u as 
    a numpy array 
    '''
    
    coords = u.getFunctionSpace().getX()     
    
    xy = np.array(coords.toListOfTuples())
    
    u_array = np.array(u.toListOfTuples())
    
    assert len(u_array.shape) == 1
        
    
    return xy, u_array


def toRegGrid(u, nx=100, ny=100):
    
    '''
    returns a nx x ny grid representation of the escript object u
    '''

    # get coordinates 
    coords = u.getFunctionSpace().getX()     
    xy, u_array = convert_to_array(u)
    
    # create regular grid
    xi = np.linspace( es.inf(coords[0]), es.sup(coords[0]), nx )
    yi = np.linspace( es.inf(coords[1]), es.sup(coords[1]), ny )
    xg, yg = np.meshgrid(xi,yi)
    
    # interpolate data
    zi = scipy.interpolate.griddata(xy, u_array, (xg,yg), method='cubic')
    
    return xi, yi, zi
    

def readRasterFile(filename):

    '''
    read gdal-compatible raster file and convert to numpy array
    
    Elco Luijendijk
    last edit: 13 feb 2013
    '''
    
    
    if os.path.isfile(filename):
        dataset=gdal.Open(filename,gdal.GA_ReadOnly)
    else:
        print 'error, could not open file %s' %filename
        return None, None, None, None
    
    print '\tnumber of raster bands:',dataset.RasterCount
    inband = dataset.GetRasterBand(1)
    geotransform = dataset.GetGeoTransform()
    nodata = inband.GetNoDataValue()
    origin = [geotransform[0], geotransform[3]]
    pixelsize = [geotransform[1], geotransform[5]]
    
    print '\torigin = (',geotransform[0], ',',geotransform[3],')'
    print '\tpixel Size = (',geotransform[1], ',',geotransform[5],')'
    print '\tdimensions: x= %s, y= %s' %(inband.XSize,inband.YSize)
    print '\tstart reading raster file'
    
    rasterArray = dataset.ReadAsArray()
    
    print '\tfinished reading raster file'
    print '\tmin,max data values = %0.1f - %0.1f' %(rasterArray.min(),rasterArray.max())
    
    return rasterArray, origin , pixelsize, nodata


def drain_function_specified_h( h, zs, specified_h_loc_init, h_buffer=0):
    '''
    drain function. implements drain by setting a specified head boundary 
    condition where hydraulic head exceeds the surface elevation.
    
    '''
    
    # set level of drain (can be above or below surface level)
    drain_level = zs + h_buffer
    
    # recalculate difference hydraulic head and surface level, with extra buffer
    h_drain = h - drain_level
    
    # set specified h at drain locs
    specified_h_loc = specified_h_loc_init + es.wherePositive( h_drain ) + es.whereZero( h_drain)
    
    # correct double counted drain nodes:
    specified_h_loc = es.wherePositive( specified_h_loc )
    
    # set specified hydraulic head values
    specified_h = specified_h_loc_init * drain_level +\
                    es.wherePositive( h_drain ) * drain_level +\
                    es.whereZero( h_drain) * drain_level
    
    return specified_h_loc, specified_h


def drain_function_flux(drain_h, drain_node_old, drain_conductance, 
                        drain_adjustment_factor_upward, 
                        drain_adjustment_factor_downward,
                        dz,
                        Q_drain_old):
    
    '''
    drain flux algorithm
    
    still experimental, drain flux boundary basically doesnt give reasonable
    fluxes or heads for large scale models, still experimenting how to fix this
    '''
    
    
    drain_node = es.wherePositive(drain_h)

    drain_switch = es.whereNegative(drain_node - drain_node_old)

    drain_steady = es.whereZero(drain_node - drain_node_old) * es.whereNonZero(drain_node)

    drain_new = es.wherePositive(drain_node - drain_node_old)

    # find out where new drain nodes are
    new_drain_nodes = es.wherePositive(drain_h) * es.whereNegative(drain_switch)

    # calculate drain flux for new nodes
    Q_drain_new_nodes = drain_new * -es.wherePositive(drain_h) * drain_h * drain_conductance / dz

    # calculate new drain flux fro nodes that switch from drain to non-drain in last iteration
    Q_drain_switch_nodes = drain_switch * Q_drain_old * drain_adjustment_factor_downward

    # 
    Q_drain_old_nodes = drain_steady * Q_drain_old * parameters.drain_adjustment_factor_upward
    
    Q_drain =  Q_drain_old_nodes + Q_drain_switch_nodes + Q_drain_new_nodes
    
    drain_node = drain_node_old + drain_new
    
    return Q_drain, drain_node
    

def filled_cnt_plot(ax, cf_var, cnt_int = 10):
    '''
    '''

    xi_cf, yi_cf, zi_cf = toRegGrid(cf_var)
    c = np.linspace(zi_cf.min(), zi_cf.max(), cnt_int)
    cf = ax.contourf( xi_cf, yi_cf, zi_cf,  cnt_int)
    cn = ax.contour( xi_cf, yi_cf, zi_cf,  cnt_int, colors='black')
    
    return cf


def result_plot( ax, cf_var, c_var, node_var , v,
                cnt_int = 10, scale_factor=5, v_int=1, sc_int = 1,
                transpose=True, log_scale_v=True):
    '''
    function for generic model results fig. escript variables
    
    input parameters:
    cf_var
        variable for filled colour contour
    c_var
        var for contour lines
    node_var    
        variable for scatter plot
    v
        variable for flux arrows
    '''
    
    
    if cf_var != None:
        xi_cf, yi_cf, zi_cf = toRegGrid(cf_var)
        c = np.linspace(zi_cf.min(), zi_cf.max(), cnt_int)
        cf = ax.contourf( xi_cf, yi_cf, zi_cf,  cnt_int)
    else:
        cf = None
    
    if c_var != None:
        xi_c, yi_c, zi_c = toRegGrid(c_var)
        cn = ax.contour( xi_c, yi_c, zi_c,  cnt_int, colors='black')
    else:
        cn = None
        
    if node_var != None:
        xy, node_value = convert_to_array( node_var )
        sc = ax.scatter(    xy[::sc_int,0], xy[::sc_int,1],
                            c=node_value[::sc_int], alpha=0.5, edgecolor='None' )
    else:
        sc = None
    
    if v != None:
        xyv, v_vector = convert_to_array(v)
        v_scalar = (v_vector[:,0]**2 + v_vector[:,1] **2) **0.5
        
        if log_scale_v == True:
            # calculate angle & length arrow
            
            # convert length to log scale
            
            # recalculate vx and vy from length and angle
            pass
        
        scale = np.abs(v_scalar).max() * scale_factor
        qv = ax.quiver( xyv[::v_int,0], xyv[::v_int,1],
                       v_vector[::v_int,0], v_vector[::v_int,1],
                       angles='xy', scale=scale )
    else:
        qv = None
        
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    return cf, cn, sc, qv


def draw_model_elements(ax, escript_var, target_depth=0): 
    
    '''
    create color polygon for each element
    
    code copied from sutra scripts
    
    assumes regular (square) grid
    and creates a rectangle around each grid node
    
    '''
    
    xyz_inp, model_var = convert_to_array(escript_var)
       
    _, mesh_dimensions = xyz_inp.shape
    
    if mesh_dimensions == 3:
        print 'slicing 3D variable'
        model_var = model_var[xyz_inp[:,2]==target_depth]
        xyz = xyz_inp[xyz_inp[:,2]==target_depth][:,:2]
    else:
        xyz = xyz_inp
    
    
    dx = np.unique(xyz[:, 0])[1] - np.unique(xyz[:, 0])[0]
    dy = np.unique(xyz[:, 1])[1] - np.unique(xyz[:, 1])[0]
        
    x_coords = [[x - 0.5*dx, x - 0.5*dx, x + 0.5*dx, x + 0.5*dx] for x in xyz[:, 0]]
    y_coords = [[y - 0.5*dy, y + 0.5*dy, y + 0.5*dy, y - 0.5*dy] for y in xyz[:, 1]]
    
    # create polygons for elements:
    patches = [matplotlib.patches.Polygon(np.array([x_coord, y_coord]).T, True)
                for x_coord, y_coord in zip(x_coords, y_coords)]
    
    # add element values:
    elements = matplotlib.collections.PatchCollection(
                                                patches,
                                                linewidths=0.0,
                                                antialiased = False)
 
    elements.set_array(model_var)
    elements.set_clim(model_var.min(), model_var.max())
          
    elements.set_edgecolor('None')
    elements.set_zorder(0.1)
    ax.add_collection(elements)
    
    ax.set_xlim(xyz[:,0].min()-dx, xyz[:,0].max()+dx )
    ax.set_ylim(xyz[:,1].min()-dy, xyz[:,1].max()+dy )

    return elements
