import bempp.api
from bempp.api.file_interfaces import FileReader
import numpy as np
import meshio
import pixsim.potentials as potentials
import pixsim.bem as bem

def solve_boundary(mshfile, **kwds):
    '''
    Solve for the nuemann coefficients using bempp.
    This assumes the surfaces have been defined using GMSH
    and are labled by 'anode', 'cathode', 'walls', etc.
    '''

    # setting bem accuracy knobs
    kwds = bem.knobs(**kwds)

    # import the mesh
    grid = bempp.api.import_grid(mshfile)
    mesh = meshio.read(mshfile)

    # try to extract high level information from mesh
    # e.g. domains, tpc dimensions
    domains = dict()
    for n,v in mesh.field_data.iteritems():
      domains[v[0]] = n
    height        = abs( mesh.points[:,1].max()-mesh.points[:,1].min() ) # cm
    width         = abs( mesh.points[:,2].max()-mesh.points[:,2].min() ) # cm
    drift_length  = abs( mesh.points[:,0].max()-mesh.points[:,0].min() ) # cm
    
    drift_pot = potentials.field_cage(domains, drift_length, **kwds)

    # define operators
    piecewise_lin_space, piecewise_const_space = bem.get_spaces(grid)
    identity, dlp, slp = bem.get_operators(grid)

    # solve
    print 'Solving...'
    dirichlet_fun = bempp.api.GridFunction(piecewise_const_space, fun=drift_pot)
    rhs = dirichlet_fun
    lhs = slp
    #neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)
    sol, info, residuals = bempp.api.linalg.gmres(slp, rhs, tol=1E-6, return_residuals=True, use_strong_form=True)

    from pixsim.models import Array
    return [ Array(typename='scalar', name='dirichlet', data=dirichlet_fun.coefficients), 
             Array(typename='scalar', name='neumann', data=sol.coefficients) ]
