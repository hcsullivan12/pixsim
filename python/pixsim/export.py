import meshio
import numpy 

def save_boundary_vtk(result, mshfile, outname):
    '''
    Save a boundary result into a VTK file.
    '''
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

def save_raster_vtk(result, outname):
    '''
    Save a drift result into a VTK file.
    '''
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

def save_step_vtk(result, outname):
    '''
    Save a step result into a VTK file.
    '''

    paths = list(result.data)
    points = list()
    vel    = list()
    for path in paths:
        for pt,v in zip(path.data[:,0:3],path.data[:,3:4]):
            points.append(pt)
            vel.append(v)
    points = numpy.asarray(points)
    vel = numpy.asarray(vel)

    from tvtk.api import tvtk, write_data
    ug = tvtk.UnstructuredGrid()    
    point_type = tvtk.Vertex().cell_type
    npoints = len(points)
    cell_types = numpy.array([point_type]*npoints)
    cell_array = tvtk.CellArray()
    cells = numpy.array([npoints]+range(npoints))
    cell_array.set_cells(point_type, cells)
    
    ug.set_cells(1, cell_array)
    ug.points = points
    ug.point_data.scalars = vel
    ug.point_data.scalars.name = 'velocity'
    
    fname = '%s-%s.vtk' % (outname, 'paths')
    write_data(ug, fname)

def export(result, save, mshfile, outname):
    '''
    Determine which exporter to call.
    '''
    if save == 'raster':
        save_raster_vtk(result, outname)
    elif save == 'boundary':
        save_boundary_vtk(result, mshfile, outname)
    elif save == 'step':
        save_step_vtk(result, outname)
    else:
        print 'No methods implented for',save
