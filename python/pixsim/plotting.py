import bempp.api
import pixsim.bem as bem
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
    

def plot_potential(mshfile, bins, points, pot, height=0.2, draw_potential=1, draw_contours=1, **kwds):
    '''
    Area to plot solution.
    '''
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt

    extents = [ [np.min(points[i]), np.max(points[i])] for i in range(0,3) ]
    npoints = [ len(bins[i]) for i in range(0,3) ]
    grid_x, grid_y, grid_z = np.mgrid[extents[0][0]:extents[0][1]:500j, 
                                      extents[1][0]:extents[1][1]:1j, 
                                      extents[2][0]:extents[2][1]:500j]

    # plotting in xz plane
    points = points.T
    pot = pot.reshape(points.shape[0],1)
    interp = griddata(points, pot, (grid_x, height, grid_z), method='linear')
    plt.subplot(121)
    data = interp[:,0,:,0]
    plt.figure(figsize = (20,20))
    if draw_potential:
        plt.imshow(data,cmap='jet', extent=(extents[2][0],extents[2][1],extents[0][0],extents[0][1]),origin='lower')
    plt.xlabel("z [cm]",fontsize=20)
    plt.ylabel("x [cm]",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    if draw_contours:
        contours = plt.contour(grid_z[:,0,:], grid_x[:,0,:], data, [0.010,0.017,0.025,0.050,0.075,0.1,0.2,0.3], colors='black')
        plt.clabel(contours, inline=True, fontsize=8)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(cax=cax)
    cb.set_label('potential [V]', fontsize=20)
    plt.show()
    
def plot_waveforms(waveforms, **kwds):
    wvfs = np.asarray(waveforms)
    #wvfs = wvfs[0:-1:2]

    import matplotlib.pyplot as plt    
    from scipy.interpolate import make_interp_spline, BSpline
    for w in wvfs:
        x = np.asarray(w[:,0])
        y = np.asarray(w[:,1])
        xnew = np.linspace(x.min(),x.max(), 5*len(x)) 
        spl = make_interp_spline(x, y, k=3)
        smooth = spl(xnew)
        plt.plot(xnew,smooth)

    plt.xlabel("t [us]",fontsize=20)
    plt.ylabel("current [arb]",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()