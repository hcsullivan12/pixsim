import meshio
import numpy 
from pixsim.store import get_result
from pixsim.vector import Scalar

def save_boundary_vtk(ses, mshfile, outname):
    '''
    Save a boundary result into a VTK file.
    '''
    print 'Saving boundary results...'
    res_id = input('Enter the boundary result ID: ')
    result = get_result(ses, None, res_id)
    if result is None:
        print 'No matching results for ID = {}'.format(res_id)
        return
    mesh = meshio.read(mshfile)
    barrs = result.array_data_by_name()

    from tvtk.api import tvtk, write_data
    pd = tvtk.PolyData()
    pd.points = mesh.points
    pd.polys = mesh.cells['triangle']
    pd.cell_data.add_array(mesh.cell_data['triangle']['gmsh:physical'])
    pd.cell_data.get_array(0).name = 'domains'
    for count, name in enumerate(['dirichlet', 'neumann']):
        pd.cell_data.add_array(barrs[name])
        pd.cell_data.get_array(count + 1).name = name
    outname = outname+'.vtk'
    write_data(pd, outname)

def save_raster_vtk(ses, outname):
    '''
    Save a drift result into a VTK file.
    '''
    print 'Saving raster results...'
    res_id = input('Enter the raster result ID: ')
    result = get_result(ses, None, res_id)
    if result is None:
        print 'No matching results for ID = {}'.format(res_id)
        return
    arrs = result.array_data_by_name()
    points = arrs['points'].T

    from tvtk.api import tvtk, write_data

    for thing in ['scalar','gradient']:
        values = arrs[thing]
        npoints = len(points)

        ug = tvtk.UnstructuredGrid()
        point_type = tvtk.Vertex().cell_type

        cell_types = numpy.array([point_type]*npoints)
        cell_array = tvtk.CellArray()
        cells = numpy.array([npoints]+range(npoints))
        cell_array.set_cells(point_type, cells)

        ug.set_cells(1, cell_array)
        if thing == 'scalar':
            ug.points = points
            ug.point_data.scalars = values.reshape(npoints)
            ug.point_data.scalars.name = thing
        else:
            ug.points = points
            field = numpy.asarray([[x,y,z] for x,y,z in zip(values[0,:,:,:].reshape(npoints),values[1,:,:,:].reshape(npoints),values[2,:,:,:].reshape(npoints))])
            ug.point_data.vectors = field
            ug.point_data.vectors.name = thing

        fname = '%s-%s.vtk' % (outname, thing)
        write_data(ug, fname)

def save_step_vtk(ses, outname):
    '''
    Save a step result into a VTK file.
    '''
    print 'Saving stepping results...'
    step_id = input('Enter the step result ID: ')
    step_res = get_result(ses, None, step_id)
    if step_res is None:
        print 'No matching results for ID = {}'.format(step_id)
        return

    # we need the raster data
    vel_res = get_result(ses, None, step_res.parent_id)
    if vel_res is None:
        print 'No matching results for ID = {}'.format(step_res.parent_id)
        return
    rast_res = get_result(ses, None, vel_res.parent_id)
    if rast_res is None:
        print 'No matching results for ID = {}'.format(vel_res.parent_id)
        return
    # get the field and linspaces
    field, linspaces = None, None
    for arr in rast_res.data:
        if arr.typename == 'scalar':
            field = arr.data
        if arr.typename == 'linspace':
            linspaces = arr.data
    assert(field is not None and linspaces is not None)

    field = Scalar(field, linspaces)

    # exporting initial position
    vtxs = [x for x in step_res.data if 'vtx' in x.name]
    points = list()
    pot    = list()
    for vtx in vtxs:
        for pt in vtx.data:
            points.append(pt)
            pot.append(field(pt))
    points = numpy.asarray(points)
    pot = numpy.asarray(pot)

    from tvtk.api import tvtk, write_data
    if len(vtxs):
        ug = tvtk.UnstructuredGrid()    
        point_type = tvtk.Vertex().cell_type
        npoints = len(points)
        cell_types = numpy.array([point_type]*npoints)
        cell_array = tvtk.CellArray()
        cells = numpy.array([npoints]+range(npoints))
        cell_array.set_cells(point_type, cells)

        ug.set_cells(1, cell_array)
        ug.points = points
        ug.point_data.scalars = pot
        ug.point_data.scalars.name = 'potential'

        fname = '%s-%s.vtk' % (outname, 'vtxs')
        write_data(ug, fname)

    # exporting paths
    paths = [x for x in step_res.data if 'path' in x.name]
    points = list()
    pot    = list()
    vel    = list()
    for path in paths:
        for pt in path.data:
            xyz = pt[0:3]
            uvw = pt[3:6]
            points.append(xyz)
            pot.append(field(xyz))
            mag = numpy.sqrt(uvw.dot(uvw))
            vel.append(mag)
    points = numpy.asarray(points)
    pot = numpy.asarray(pot)
    vel = numpy.asarray(vel)

    if len(paths):
        for flv,data in {'potential':pot, 'velocity':vel}.iteritems():
            ug = tvtk.UnstructuredGrid()    
            point_type = tvtk.Vertex().cell_type
            npoints = len(points)
            cell_types = numpy.array([point_type]*npoints)
            cell_array = tvtk.CellArray()
            cells = numpy.array([npoints]+range(npoints))
            cell_array.set_cells(point_type, cells)

            ug.set_cells(1, cell_array)
            ug.points = points
            ug.point_data.scalars = data
            ug.point_data.scalars.name = flv

            fname = '%s-%s-%s.vtk' % (outname, 'paths', flv)
            write_data(ug, fname)

def export(ses, mshfile, save, outname):
    '''
    Determine which exporter to call.
    '''
    if save == 'raster':
        save_raster_vtk(ses, outname)
    elif save == 'boundary':
        save_boundary_vtk(ses, mshfile, outname)
    elif save == 'step':
        save_step_vtk(ses, outname)
    else:
        print 'No methods implented for',save
