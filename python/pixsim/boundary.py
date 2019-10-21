"""
Solve boundary value problem using BEM++.
"""
import meshio
import bempp.api
import pixsim.potentials as potentials
import pixsim.bem as bem

def boundary(mshfile, method=None, **kwds):
    """Solve for the nuemann coefficients using bempp."""
    # setting bem accuracy knobs
    #kwds = bem.knobs(**kwds)

    # import the mesh
    grid = bempp.api.import_grid(mshfile)
    mesh = meshio.read(mshfile)

    # try to extract high level information from mesh
    # e.g. domains, tpc dimensions
    domains = dict()
    for name, did in mesh.field_data.iteritems():
        domains[did[0]] = name

    # use the specified potential
    drift_pot = None
    if method == 'field_cage':
        drift_pot = potentials.field_cage(domains, **kwds)
    elif method == 'weighting':
        drift_pot = potentials.weighting(**kwds)
    else:
        raise ValueError('No boundary method named', method)

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
    return [Array(typename='scalar', name='dirichlet', data=dirichlet_fun.coefficients),
            Array(typename='scalar', name='neumann', data=sol.coefficients)]
