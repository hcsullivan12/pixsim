import bempp.api

def knobs(hmat_eps = None, hmat_mbs = None, gqo_near = None, gqo_medium = None, gqo_far = None, gqo_ds = None, **kwds):
    '''
    Set the BEM++ knobs that can be set.  This returns filtered kwds
    and will fill in default values if they are not set.

    Solutions may have sub-domains where the potential is
    discontinuous on their borders.  Increasing the near, medium or
    far orders may help produce a more correct and smooth solution.

    http://www.bempp.org/quadrature.html?highlight=global_parameters

    hmat_eps=1e-3 is default.  Increase this if get artifacts in
    calculated Neumann boundary conditions.  Eg, 1e-5.

    hmat_mbs=20 is default.  HMAT minimum block size.  Increase this
    if you get "tearing" in the off-surface potentials.  Eg, 1000.
    Set this just before assembly of the potential operator.

    @todo: add setting of distance scales.
    '''
    import bempp.api
    q = bempp.api.global_parameters.quadrature
    if gqo_ds:
        q.double_singular = gqo_ds
    if gqo_near:
        q.near.single_order = gqo_near
        q.near.double_order = gqo_near
    if gqo_medium:
        q.medium.single_order = gqo_medium
        q.medium.double_order = gqo_medium
    if gqo_far:
        q.far.single_order = gqo_far
        q.far.double_order = gqo_far

    h = bempp.api.global_parameters.hmat
    if hmat_eps:
        h.eps = 1E-5
    if hmat_mbs:
        h.min_block_size = hmat_mbs

    kwds['gqo_ds'] = q.double_singular
    kwds['gqo_near'] = q.near.single_order
    kwds['gqo_medium'] = q.medium.single_order
    kwds['gqo_far'] = q.far.single_order
    kwds['hmat_eps'] = h.eps
    kwds['hmat_mbs'] = h.min_block_size

    print 'Gaussian quadrature orders:',q.double_singular,q.near.single_order,q.medium.single_order,q.far.single_order
    print 'HMAT eps/mbs:',h.eps, h.min_block_size
    return kwds

def get_spaces(grid):
    return bempp.api.function_space(grid, "P", 1), bempp.api.function_space(grid, "DP", 0), 

def get_operators(grid):
    piecewise_lin_space, piecewise_const_space = get_spaces(grid)

    identity = bempp.api.operators.boundary.sparse.identity(piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    dlp = bempp.api.operators.boundary.laplace.double_layer(piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
    slp = bempp.api.operators.boundary.laplace.single_layer(piecewise_const_space, piecewise_const_space, piecewise_const_space)
    return identity, dlp, slp