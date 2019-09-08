import pygmsh
import numpy as np

class BoxTpc:
    def __init__(self, geofile, **kwds):
        self.geofile = geofile
        self.kwds = kwds

    def construct(self):
        return self.construct_geometry(self.kwds)

    def initPixels(tpc_dim     = (2,5,5),
                   pad_spacing = 0.4, 
                   grid_size   = (6,6),
                   pad_offset  = 0.0):
        '''
        Initialize the positions of the pixels. 
        This builds an array of pixels that is symmetric 
        about the y and z axes.
        '''
        currentZ = 0.5*pad_spacing
        pixels_z, pixels_y = list(), list()
        while len(pixels_z) < grid_size[0]:
            pixels_z.append(currentZ)
            pixels_z.append(-1*currentZ)
            currentZ += pad_spacing
        currentY = 0.5*pad_spacing
        while len(pixels_y) < grid_size[1]:
            pixels_y.append(currentY)
            pixels_y.append(-1*currentY)
            currentY += pad_spacing
        pixels_z.sort()
        pixels_y.sort(reverse=True)
        # translate to 0 < z < length
        for i,p in enumerate(pixels_z):
            pixels_z[i]+=0.5*tpc_dim[2]

        pixels = list()
        for y in pixels_y:
            for z in pixels_z:
                pixels.append(np.asarray((pad_offset,y,z)))
        return pixels

    def construct_pixels(active, geom, pixels,
                         pad_rmax   = 0.04,
                         grid_pitch = 0.04,
                         pad_shape  = 'box'):
        '''
        Placing pixels
        '''
        for x,y,z in pixels:
            if pad_shape == 'box':
                pad = geom.add_box( [x, y-pad_rmax, z-pad_rmax], 
                                    [grid_pitch, 2*pad_rmax, 2*pad_rmax] )
                active = geom.boolean_difference([active], [pad])
            elif pad_shape == 'circle':
                pad = geom.add_cylinder( [x, y, z], 
                                         [grid_pitch,0,0],
                                         pad_rmax )
                active = geom.boolean_difference([active], [pad])
            else:
                raise ValueError('Shape not supported: %s' % pad_shape)
        return active

    def construct_grid(active, geom, pixels, 
                       tpc_dim = (2,5,5), 
                       pad_spacing  = 0.4,
                       pad_offset   = 0.0,
                       grid_pitch   = 0.04,
                       grid_dim     = (6,6)):
        '''
        Constructing focusing grid
        '''
        grid_extent_z = grid_dim[0]*pad_spacing+grid_pitch
        grid_extent_y = grid_dim[1]*pad_spacing+grid_pitch
        grid_corner_z = 0.5*(tpc_dim[2] - grid_extent_z)
        grid_corner_y = 0.5*(tpc_dim[1] - grid_extent_y)
        print grid_corner_z, grid_extent_z
        focus_grid = geom.add_box( [pad_offset, -0.5*tpc_dim[1]+grid_corner_y, grid_corner_z], 
                                   [grid_pitch, grid_extent_y, grid_extent_z])
        for x,y,z in pixels:
            grid_off = 0.5*(pad_spacing-grid_pitch)
            box = geom.add_box( [x, y-grid_off, z-grid_off], 
                                [grid_pitch, pad_spacing-grid_pitch, pad_spacing-grid_pitch])
            focus_grid = geom.boolean_difference([focus_grid], [box])
        return geom.boolean_difference([active], [focus_grid])

    def construct_geometry(tpc_dim=(2.0,10.,10.), build_pads=0):
        '''
        Builds geometry and writes .geo file for GMSH
        '''

        # initialize pixels
        pixels = self.initPixels()
        print 'Initialized',len(pixels),'pixels...'

        # start constructing the geometry
        geom = pygmsh.opencascade.Geometry()
        # active
        active = geom.add_box( [0, -0.5*tpc_dim[1], 0], tpc_dim )
        # pixels
        if build_pads:
            active = self.construct_pixels(active, geom, pixels)
        # grid
        if build_pads:
            active = self.construct_grid(active, geom, pixels, tpc_dim)

        mesh = pygmsh.generate_mesh(geom, geo_filename=self.geofile)

        from pixsim.models import Result, Array
        return [ Array(typename='points', name='pixels', data=np.asarray(pixels)) ]

class FieldCageTpc:
    def __init__(self, geofile):
        self.geofile = geofile

    def initPixels(self, pad_spacing = 0.4, 
                   grid_size   = (6,6),
                   pad_offset  = 0.0,
                   **kwds):
        '''
        Initialize the positions of the pixels. 
        This builds an array of pixels that is symmetric 
        about the y and z axes.
        '''
        currentZ = 0.5*pad_spacing
        pixels_z, pixels_y = list(), list()
        while len(pixels_z) < grid_size[0]:
            pixels_z.append(currentZ)
            pixels_z.append(-1*currentZ)
            currentZ += pad_spacing
        currentY = 0.5*pad_spacing
        while len(pixels_y) < grid_size[1]:
            pixels_y.append(currentY)
            pixels_y.append(-1*currentY)
            currentY += pad_spacing
        pixels_z.sort()
        pixels_y.sort(reverse=True)

        pixels = list()
        for y in pixels_y:
            for z in pixels_z:
                pixels.append(np.asarray((pad_offset,y,z)))
        return pixels

    def construct_pixels(self, geom, pixels,
                         pad_rmax   = 0.04,
                         grid_pitch = 0.04,
                         pad_shape  = 'box',
                         **kwds):
        '''
        Placing pixels
        '''
        for x,y,z in pixels:
            if pad_shape == 'box':
                print x, y-pad_rmax,z-pad_rmax
                pad = geom.add_box( [x, y-pad_rmax, z-pad_rmax], 
                                    [grid_pitch, 2*pad_rmax, 2*pad_rmax] )
            elif pad_shape == 'circle':
                pad = geom.add_cylinder( [x, y, z], 
                                         [grid_pitch,0,0],
                                         pad_rmax )
            else:
                raise ValueError('Shape not supported: %s' % pad_shape)
        return geom

    def construct_grid(self, geom, pixels, 
                       tpc_dim = (2,5,5), 
                       pad_spacing  = 0.4,
                       pad_offset   = 0.0,
                       grid_pitch   = 0.04,
                       grid_dim     = (6,6),
                       **kwds):
        '''
        Constructing focusing grid
        '''
        grid_extent_z = grid_dim[0]*pad_spacing+grid_pitch
        grid_extent_y = grid_dim[1]*pad_spacing+grid_pitch
        grid_corner_z = 0.5*(tpc_dim[2] - grid_extent_z)
        grid_corner_y = 0.5*(tpc_dim[1] - grid_extent_y)
        focus_grid = geom.add_box( [pad_offset, -0.5*tpc_dim[1]+grid_corner_y, grid_corner_z], 
                                   [grid_pitch, grid_extent_y, grid_extent_z])
        for x,y,z in pixels:
            grid_off = 0.5*(pad_spacing-grid_pitch)
            box = geom.add_box( [x, y-grid_off, z-grid_off], 
                                [grid_pitch, pad_spacing-grid_pitch, pad_spacing-grid_pitch])
        return geom

    def construct_geometry(self, tpc_dim=(0.,10.,20.), electrode_thickness=0.5, 
                           build_pads=0, ring_thickness=1.0, nrings=10, **kwds):
        '''
        Builds geometry and writes .geo file for GMSH
        '''
        # initialize pixels
        pixels = self.initPixels(**kwds)
        print 'Initialized',len(pixels),'pixels...'
        # start constructing the geometry
        geom = pygmsh.opencascade.Geometry()
        # cathode and anode
        anode   = geom.add_cylinder([-0.5*electrode_thickness,0,0], [electrode_thickness,0,0], tpc_dim[1])
        cathode = geom.add_cylinder([tpc_dim[2]-0.5*electrode_thickness,0,0], [electrode_thickness,0,0], tpc_dim[1])
        # field shaping rings
        increment = tpc_dim[2]/(nrings+1)
        xpos = increment-ring_thickness*0.5
        for r in range(0, nrings):
            print 0,tpc_dim[2],xpos
            rin = geom.add_cylinder([xpos,0,0], [ring_thickness,0,0], tpc_dim[1]-ring_thickness)
            rout = geom.add_cylinder([xpos,0,0], [ring_thickness,0,0], tpc_dim[1])
            ring = geom.boolean_difference([rout], [rin])
            xpos += increment

        # pixels
        if build_pads:
            geom = self.construct_pixels(geom, pixels, **kwds)
        # grid
        if build_pads:
            geom = self.construct_grid(geom, pixels, **kwds)

        mesh = pygmsh.generate_mesh(geom, geo_filename=self.geofile)

        from pixsim.models import Result, Array
        return [ Array(typename='points', name='pixels', data=np.asarray(pixels)) ]


