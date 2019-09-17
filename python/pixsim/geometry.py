import pygmsh
import numpy as np

def make_pixels(pad_spacing   = 0.4, 
                grid_dim      = (6,6),
                pad_offset    = 0.0,
                pad_rmax      = 0.04,
                grid_diameter = 0.04,
                pad_shape     = 'box',
                **kwds):
    '''
    Initialize pixels. 
    '''
    currentZ = 0.5*pad_spacing
    pixels_z, pixels_y = list(), list()
    while len(pixels_z) < grid_dim[0]:
        pixels_z.append(currentZ)
        pixels_z.append(-1*currentZ)
        currentZ += pad_spacing
    currentY = 0.5*pad_spacing
    while len(pixels_y) < grid_dim[1]:
        pixels_y.append(currentY)
        pixels_y.append(-1*currentY)
        currentY += pad_spacing
    pixels_z.sort()
    pixels_y.sort(reverse=True)
    
    pixels = list()
    for y in pixels_y:
        for z in pixels_z:
            pixels.append(np.asarray((pad_offset,y,z)))  

    import pixsim.geo as geoobj
    geo = list()
    for count,cent in enumerate(pixels,1):
        pad = None
        name = 'pixel'+str(count)
        if pad_shape == 'box':
            pad = geoobj.Box(name,0.5*grid_diameter,pad_rmax,pad_rmax,cent)
        elif pad_shape == 'cylinder':
            pad = geoobj.Cylinder(name,pad_rmax,0.5*grid_diameter,cent,np.asarray([1,0,0]))
        elif pad_shape == 'sphere':
            pad = geoobj.Sphere(name,pad_rmax,cent)
        else:
            raise ValueError('Shape not supported: %s' % pad_shape)
        geo.append(pad)
    return geo

def wire_grid(geom, pixels, 
              pad_spacing  = 0.4,
              pad_offset   = 0.0,
              grid_diameter= 0.04,
              grid_dim     = (6,6),
              **kwds):
    '''
    Focusing grid structure made out of wires.
    '''
    nelem = 0
    length = grid_dim[0]*pad_spacing
    sp = [pad_offset,0,-0.5*length]
    center = geom.add_cylinder(sp, [0,0,length], 0.5*grid_diameter)
    for ywire in range(1, int(0.5*grid_dim[0])+1):
        sp = [pad_offset,ywire*pad_spacing,-0.5*length]
        wire = geom.add_cylinder(sp, [0,0,length], 0.5*grid_diameter)
        sp = [pad_offset,-1*ywire*pad_spacing,-0.5*length]
        wire = geom.add_cylinder(sp, [0,0,length], 0.5*grid_diameter)
        nelem += 2
    
    zpos = 0
    ypos = 0.5*pad_spacing
    length = pad_spacing-0.05
    ypos = set([pix.center[1] for pix in pixels])
    print grid_dim, type(grid_dim), kwds
    zpos = set([n*pad_spacing for n in range(0,int(0.5*grid_dim[1])+1)])
    temp = set([-1*z for z in zpos])
    zpos = zpos | temp
    for y in ypos:
        for z in zpos:
            sp = [pad_offset,y-0.5*length, z]
            wire = geom.add_cylinder(sp, [0,length,0], 0.5*grid_diameter)
            nelem += 1
    return geom, nelem

def box_grid_old(active, geom, pixels, 
             tpc_dim = (2,5,5), 
             pad_spacing  = 0.4,
             pad_offset   = 0.0,
             grid_diameter   = 0.04,
             grid_dim     = (6,6),
             **kwds):
    '''
    Focusing grid structure made out of boxes.
    '''
    grid_extent_z = grid_dim[0]*pad_spacing+grid_diameter
    grid_extent_y = grid_dim[1]*pad_spacing+grid_diameter
    grid_corner_z = 0.5*(tpc_dim[2] - grid_extent_z)
    grid_corner_y = 0.5*(tpc_dim[1] - grid_extent_y)
    print grid_corner_z, grid_extent_z
    focus_grid = geom.add_box( [pad_offset, -0.5*tpc_dim[1]+grid_corner_y, grid_corner_z], 
                               [grid_diameter, grid_extent_y, grid_extent_z])
    for x,y,z in pixels:
        grid_off = 0.5*(pad_spacing-grid_diameter)
        box = geom.add_box( [x, y-grid_off, z-grid_off], 
                            [grid_diameter, pad_spacing-grid_diameter, pad_spacing-grid_diameter])
        focus_grid = geom.boolean_difference([focus_grid], [box])
    return geom.boolean_difference([active], [focus_grid])

def construct_pixels_old(active, geom, pixels,
                     pad_rmax   = 0.04,
                     grid_diameter = 0.04,
                     pad_shape  = 'box',
                     **kwds):
    '''
    Placing pixels
    '''
    for x,y,z in pixels:
        if pad_shape == 'box':
            pad = geom.add_box( [x, y-pad_rmax, z-pad_rmax], 
                                [grid_diameter, 2*pad_rmax, 2*pad_rmax] )
            active = geom.boolean_difference([active], [pad])
        elif pad_shape == 'circle':
            pad = geom.add_cylinder( [x, y, z], 
                                     [grid_diameter,0,0],
                                     pad_rmax )
            active = geom.boolean_difference([active], [pad])
        else:
            raise ValueError('Shape not supported: %s' % pad_shape)
    return active

def construct_pixels(geom, pixels,
                     pad_rmax   = 0.04,
                     grid_diameter = 0.04,
                     **kwds):
    '''
    Placing pixels
    '''
    for pix in pixels:
        _,hdim,center,shape = pix.info()
        x,y,z = center
        if shape == 'box':
            pad = geom.add_box( [x-0.5*grid_diameter, y-pad_rmax, z-pad_rmax], 
                                [grid_diameter, 2*pad_rmax, 2*pad_rmax] )
        elif shape == 'circle':
            pad = geom.add_cylinder( [x-0.5*grid_diameter, y, z], 
                                     [grid_diameter,0,0],
                                     pad_rmax )
        elif shape == 'sphere':
            pad = geom.add_ball( center, 
                                 hdim)
    return geom

def add_surface(geofile, bound, values):
    with open(geofile, 'a') as f:
        s = '\n//+\nPhysical Surface(\"'+str(bound)+'\") = {'
        for count,v in enumerate(values, 1):
            if count < len(values):
                s += str(v)+', '
            else:
                s += str(v)+'};'
        f.write(s)

class BoxTpc:
    def __init__(self, geofile, **kwds):
        self.geofile = geofile
        self.kwds = kwds

    def construct_geometry(self, pixels, tpc_dim=(2,2,2), build_pads=0, build_grid=0, **kwds):
        '''
        Builds geometry and writes .geo file for GMSH
        '''
        shape_inc = {'sphere':1, 'cylinder':3, 'box':6}

        # start constructing the geometry
        geom = pygmsh.opencascade.Geometry()
        # active
        active = geom.add_box( [0, -0.5*tpc_dim[1], -0.5*tpc_dim[2]], tpc_dim )
        # define surfaces, not the best way, but better than doing it manually
        bound = {'anode':[1], 'cathode':[2], 'walls':[3,4,5,6]}
        next_surf = 7

        # pixels and grid
        if build_pads:
            geom = construct_pixels(geom, pixels, **kwds)
            for pix in pixels:
                 name,hdim,center,shape = pix.info()
                 bound[name] = [s for s in range(next_surf,next_surf+shape_inc[shape])]
                 assert(len(bound[name]) == shape_inc[shape])
                 next_surf += shape_inc[shape]

        if build_grid:
            geom, nelem = wire_grid(geom, pixels, **kwds)
            bound['grid'] = [s for s in range(next_surf,next_surf+nelem*shape_inc['cylinder']+1)]

        mesh = pygmsh.generate_mesh(geom, geo_filename=self.geofile)
        
        # add the surfaces
        for b,v in bound.iteritems():
            add_surface(self.geofile, b, v)

class FieldCageTpc:
    def __init__(self, geofile):
        self.geofile = geofile

    def construct_geometry(self, tpc_dim=(0.,10.,20.), electrode_thickness=0.5, 
                           build_pads=0, ring_thickness=1.0, nrings=10, **kwds):
        '''
        Builds geometry and writes .geo file for GMSH
        '''
        # initialize pixels
        pixels = define_pixels(**kwds)
        print 'Initialized',len(pixels),'pixels...'
        # start constructing the geometry
        geom = pygmsh.opencascade.Geometry()
        # cathode and anode
        asp = [-0.5*electrode_thickness,0,0]
        csp = [tpc_dim[2]+asp[0],0,0]
        anode   = geom.add_cylinder(asp, [electrode_thickness,0,0], tpc_dim[1])
        cathode = geom.add_cylinder(csp, [electrode_thickness,0,0], tpc_dim[1])
        ## field shaping rings
        # using solid wall
        wsp = [-1*asp[0]+0.1,0,0]
        wl  = csp[0]-0.1
        rin  = geom.add_cylinder(wsp, [wl,0,0], tpc_dim[1]-ring_thickness)
        rout = geom.add_cylinder(wsp, [wl,0,0], tpc_dim[1])
        ring = geom.boolean_difference([rout], [rin])
        #    
        #increment = tpc_dim[2]/(nrings+1)
        #xpos = increment-ring_thickness*0.5
        #for r in range(0, nrings):
        #    rin = geom.add_cylinder([xpos,0,0], [ring_thickness,0,0], tpc_dim[1]-ring_thickness)
        #    rout = geom.add_cylinder([xpos,0,0], [ring_thickness,0,0], tpc_dim[1])
        #    ring = geom.boolean_difference([rout], [rin])
        #    xpos += increment

        # pixels
        if build_pads:
            geom = construct_pixels(geom, pixels, **kwds)
        # grid
        if build_pads:
            geom = construct_grid(geom, pixels, **kwds)

        mesh = pygmsh.generate_mesh(geom, geo_filename=self.geofile)

        from pixsim.models import Result, Array
        return [ Array(typename='points', name='pixels', data=np.asarray(pixels)) ]

class SimpleTpc:
    def __init__(self, geofile):
        self.geofile = geofile

    def construct_geometry(self, pixels, tpc_dim=(2.,4,4.), build_pads=1, **kwds):
        '''
        Builds geometry and writes .geo file for GMSH
        '''
        # start constructing the geometry
        geom = pygmsh.opencascade.Geometry()
        # cathode and anode
        asp = [-0.05,-5,-5]
        csp = [-0.05+tpc_dim[2],-5,-5]
        anode   = geom.add_box(asp,[0.1,10,10])
        cathode = geom.add_box(csp,[0.1,10,10])

        # pixels
        if build_pads:
            geom = construct_pixels(geom, pixels, **kwds)
        # grid
        if build_pads:
            geom = construct_grid_old(geom, pixels, **kwds)
        
        mesh = pygmsh.generate_mesh(geom, geo_filename=self.geofile)

        from pixsim.models import Result, Array
        return [ Array(typename='points', name='pixels', data=np.asarray(pixels)) ]



