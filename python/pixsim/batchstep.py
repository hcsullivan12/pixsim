from larf import units
import larf.util
import larf.drift
from larf.stepping import Stepper, CollectSteps, StuckDetection
from larf.geometry import RayCollision
from larf.models import Array
from larf.bem import knobs
from larf.boundary import result_to_grid_funs
from larf.raster import Points
import larf.points

import numpy
import math

class BatchedGradientFunction(object):
    def __init__(self, scalar_function, lcar=10*units.um):
        self.scalar = scalar_function
        self.lcar = lcar
        self._scalar_points_seen = list()
        self._scalars_seen = list()
        self._gradient_points_seen = list()
        self._gradients_seen = list()

    @property
    def gradient_points(self):
        return numpy.vstack(self._gradient_points_seen)
    @property
    def scalar_points(self):
        return numpy.vstack(self._scalar_points_seen)
    @property
    def scalars(self):
        return numpy.hstack(self._scalars_seen)
    @property
    def gradients(self):
        return numpy.vstack(self._gradients_seen)

    def __call__(self, points):
        '''
        Return gradient of scalar function at each of the given 3-points.
        '''
        npoints = len(points)
        if not npoints:
            return None
        call_with = list()
        indices = [
            [1,1,1],
            [2,1,1],
            [1,2,1],
            [1,1,2],
            [0,1,1],
            [1,0,1],
            [1,1,0],
        ]

        for p in points:

            call_with.append(p)

            for ind in range(3):
                pp = numpy.zeros(3)
                pp[ind] = self.lcar
                pp = p + pp
                call_with.append(pp)

            for ind in range(3):
                pm = numpy.zeros(3)
                pm[ind] = self.lcar
                pm = p - pm
                call_with.append(pm)

        ret = list()

        nsamples = len(indices)
        pots = numpy.asarray(self.scalar(*call_with))
        self._scalars_seen.append(pots)
        self._scalar_points_seen.append(call_with)
        for vec in pots.reshape(npoints, nsamples):
            phi = numpy.zeros((3,3,3))

            for ind,val in zip(indices, vec):
                phi[ind[0],ind[1],ind[2]] = val

            #print len(ret), vec
            #print phi
            grad = numpy.asarray(numpy.gradient(phi, self.lcar, self.lcar, self.lcar))
            ret.append(grad[:,1,1,1])
        grads = numpy.asarray(ret)
        self._gradients_seen.append(grads)
        self._gradient_points_seen.append(points)
        return grads

class BatchedVelocity(object):
    def __init__(self, batched_field, temperature = 89*units.Kelvin):
        self.drift = batched_field
        self.mobility = lambda emag: larf.drift.mobility(emag, temperature)
        self._points_seen = list()
        self._velocity_seen = list()

    @property
    def points(self):
        return numpy.vstack(self._points_seen)
    @property
    def velocities(self):
        return numpy.vstack(self._velocity_seen)

    def __call__(self, points):
        '''
        Return velocity vector at all the the given 3-points.
        '''
        #print 'BatchedVelocity:',len(points), points.shape
        points = numpy.asarray(points)
        points = points[:,:3]   # assure 3-point

        field = self.drift(points)
        emag = numpy.sqrt(field[:,0]**2 + field[:,1]**2 + field[:,2]**2)
        mu = self.mobility(emag)
        mu = mu.reshape((len(field),1))
        velo = mu*field
        self._points_seen.append(points)
        self._velocity_seen.append(velo)
        return velo



def step_maybe(paths, points, stop):
    '''
    For each path and point, add point to path unless it's time to stop.

    Return pair of (paths, points) with which to continue stepping. 
    '''

    next_points = list()
    next_paths = list()

    for path, p2 in zip(paths, points):
        p1 = path[-1]
        if stop(p1, p2):
            print "hit wire at %s, took %d steps starting at %s" % (p2, len(path)-1, path[0])
            continue
        # okay to continue
        path.append(p2)
        next_paths.append(path)
        next_points.append(p2)

    return numpy.asarray(next_points), next_paths


