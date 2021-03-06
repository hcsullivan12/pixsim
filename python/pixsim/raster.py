import bempp.api
import numpy as np
import pixsim.bem as bem

def linear(mshfile, sol,
           linspaces = ( (0,2,100), (-3,3,100), (2,8,100) ), **kwds):
    '''
    Evaluate the potential on a linear grid space.
    '''
    # setting bem accuracy knobs
    kwds = bem.knobs(**kwds)

    # import the mesh
    grid = bempp.api.import_grid(mshfile)

    # convert solution to grid function
    sol = bempp.api.GridFunction(grid, coefficients = sol)

    # define our points
    linspaces = [np.linspace(*ls) for ls in linspaces]
    mgrid = np.meshgrid(*linspaces, indexing='ij')
    points = np.vstack([mgrid[i].ravel() for i in range(3)])

    # evaluate on our space
    print 'Evaluating...'
    piecewise_lin_space, piecewise_const_space = bem.get_spaces(grid)
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space, points)
    u_evaluated = slp_pot * sol

    print 'u_evaluated.shape=',u_evaluated.shape, u_evaluated.T[0], points.T[0]
    u_reshaped = u_evaluated.reshape(mgrid[0].shape)
    print 'u_reshaped.shape=',u_reshaped.shape
    print mgrid[0].shape, points.shape, len(points)
    
    try:
       dxyz = [(ls[1]-ls[0]) for ls in linspaces]
       u_grad = np.asarray(np.gradient(u_reshaped, *dxyz))
       # convert efield to kV to help later
       u_grad /= 1000.
    except:
       print 'Warning: Could not compute gradient'
       u_grad = np.empty([1,1])

    from pixsim.models import Array
    return [ Array(typename='linspace', name='linspaces', data = linspaces),
             Array(typename='scalar',   name='pfield',    data = u_reshaped),
             Array(typename='vector',   name='efield',    data = u_grad),
             Array(typename='points',   name='points',    data = points) ]
