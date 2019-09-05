import bempp.api
import pixsim.operators as operators
import numpy as np

def plot_efield(mshfile, sol,
                efield_linspaces = ( (0,2,80), (-3,3,60), (2,8,60) ),
                plane     = 'xz',
                height    = 0.0,  
                **kwds):
    '''
    Area to plot solution on some domain.
    '''
    linspaces = efield_linspaces
    # import the mesh
    grid = bempp.api.import_grid(mshfile)
    # get the solution
    sol = bempp.api.GridFunction(grid, coefficients = sol)
    linsp = [np.linspace(*ls) for ls in linspaces]
    mgrid = np.meshgrid(*linsp, indexing='ij')
    points = None
    extent = None
    
    import matplotlib
    if plane == 'xy' or plane == 'yx':
        points = np.vstack( (mgrid[0].ravel(),mgrid[1].ravel(),height*np.ones(mgrid[2].size)) )
        extent = (linspaces[1][0], linspaces[1][1], linspaces[0][0], linspaces[0][1], )
    
    if plane == 'xz' or plane == 'zx':
        points = np.vstack( (mgrid[0].ravel(),height*np.ones(mgrid[1].size), mgrid[2].ravel()) )
        extent = (linspaces[2][0], linspaces[2][1], linspaces[0][0], linspaces[0][1])

    if plane == 'yz' or plane == 'zy':
        points = np.vstack( (height*np.ones(mgrid[0].size), mgrid[1].ravel(),mgrid[2].ravel()) )
        extent = (linspaces[2][0], linspaces[2][1], linspaces[1][0], linspaces[1][1])

    # evaluate on our space
    piecewise_lin_space, piecewise_const_space = bem.get_spaces(grid)
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space, points)
    u_evaluated = slp_pot * sol
    print 'u_evaluated.shape=',u_evaluated.shape, u_evaluated.T[0], points.T[0]
    u_reshaped = u_evaluated.reshape(mgrid[0].shape)
    print 'u_reshaped.shape=',u_reshaped.shape

    dxyz = [(ls[1]-ls[0])/(ls[2]-1) for ls in linsp]
    u_grad = np.asarray(np.gradient(u_reshaped, *dxyz))

    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams["image.origin"] = 'lower'

    fig, ax = plt.subplots()
    X,Y,Z = points[0,:],points[1,:],points[2,:]
    U,V,W = u_grad[0,:,:,:].reshape(X.shape[0]),u_grad[1,:,:,:].reshape(X.shape[0]),u_grad[2,:,:,:].reshape(X.shape[0])
    efield = None
    print np.sqrt(U**2+V**2+W**2)
    # resizing 
    X,Y,Z = X[0:-1:5], Y[0:-1:5], Z[0:-1:5]
    U,V,W = U[0:-1:5], V[0:-1:5], W[0:-1:5]
    if plane == 'xy' or plane == 'yx':
        efield = ax.quiver(X, Y, U, V, np.sqrt(U**2+V**2), units='xy' ,scale=10000, color='black')
    if plane == 'xz' or plane == 'zx':
        efield = ax.quiver(X, Z, U, W, np.sqrt(U**2+W**2), units='xy' ,scale=10000, color='black')
    if plane == 'yz' or plane == 'zy':
        efield = ax.quiver(Y, Z, V, W, np.sqrt(W**2+V**2), units='xy' ,scale=10000, color='black')
    #cb = plt.colorbar()
    plt.show()
    

def plot_potential(mshfile, sol,
                  potent_linspaces = ( (0,2,80), (-3,3,60), (2,8,60) ),
                  plane     = 'xz',
                  height    = 0.0,  
                  **kwds):
    '''
    Area to plot solution on some domain.
    '''
    linspaces = potent_linspaces
    # import the mesh
    grid = bempp.api.import_grid(mshfile)
    # get the solution
    sol = bempp.api.GridFunction(grid, coefficients = sol)
    linsp = [np.linspace(*ls) for ls in linspaces]
    mgrid = np.meshgrid(*linsp, indexing='ij')
    points = None
    extent = None
    
    import matplotlib
    if plane == 'xy' or plane == 'yx':
        points = np.vstack( (mgrid[0].ravel(),mgrid[1].ravel(),height*np.ones(mgrid[2].size)) )
        extent = (linspaces[1][0], linspaces[1][1], linspaces[0][0], linspaces[0][1], )
    
    if plane == 'xz' or plane == 'zx':
        points = np.vstack( (mgrid[0].ravel(),height*np.ones(mgrid[1].size), mgrid[2].ravel()) )
        extent = (linspaces[2][0], linspaces[2][1], linspaces[0][0], linspaces[0][1])

    if plane == 'yz' or plane == 'zy':
        points = np.vstack( (height*np.ones(mgrid[0].size), mgrid[1].ravel(),mgrid[2].ravel()) )
        extent = (linspaces[2][0], linspaces[2][1], linspaces[1][0], linspaces[1][1])

    # evaluate on our space
    piecewise_lin_space, piecewise_const_space = bem.get_spaces(grid)
    slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space, points)
    u_evaluated = slp_pot * sol

    matplotlib.rcParams['figure.figsize'] = (15.0, 10.0)
    matplotlib.rcParams["image.origin"] = 'lower'

    import matplotlib.pyplot as plt

    plt.imshow(u_evaluated.T, extent=extent)
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    plt.show()
    