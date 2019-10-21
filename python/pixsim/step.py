import numpy
import math
from collections import namedtuple
from pixsim.vector import Field

def step_rk4(r, t1, t2, v):
    '''
    Take a step in time from location r at time t1 to time t2 using
    velocity function v via the 4th order Runge-Kutta stepping.

    This is taken from Numerical Recipes.

    @param r: initial location
    @type r: array

    @param t1: initial time
    @type t1: float

    @param t2: time after step
    @type t2: float

    @param v: velocity function
    @type v: callable(position, time)

    @return: tuple(array, None) -- the position after the step and the
        error.

        For rk4 no error is computed.
    '''
    h = t2-t1
    k1 = h * v(t1,         r)
    k2 = h * v(t1 + 0.5*h, r + 0.5*k1)
    k3 = h * v(t1 + 0.5*h, r + 0.5*k2)
    k4 = h * v(t2,         r + k3)
               
    return r + (k1 + 2.0*(k2 + k3) + k4) / 6.0, t2+h

def step_rkck(r, t1, t2, v):
    '''
    Take a step in time from location r at time t1 to time t2 using
    velocity function v via the 4th order adaptive
    Runge-Kutta/Cash+karp stepping.

    The error returned with the step is a vector distance from the
    returned step and a second (6th order) step.  It can be used to
    optimally choose the next t2 to keep the error w/in some bounds.

    This is taken from Numerical Recipes.

    @param r: initial location
    @type r: array

    @param t1: initial time
    @type t1: float

    @param t2: time after step
    @type t2: float

    @param v: velocity function
    @type v: callable(position, time)

    @return: tuple(array, array) -- the position after the step and
        the error.
    '''

    # The Cash/Karp coefficients.  Use None placeholders to match notation
    a  = [None, None, 0.2, 0.3, 0.6, 1.0, 0.875 ]
    cs = [None, 2825.0/27648.0, 0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25]
    c  = [None, 37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0 ]
    b2 = [None, 0.2]
    b3 = [None, 3.0/40.0, 9.0/40.0]
    b4 = [None, 0.3, -0.9, 1.2]
    b5 = [None, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0]
    b6 = [None, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]
    
    h = t2-t1    
    #print '\nONE:', r
    k1 = h * v(t1, r)
    #print 'TWO:', r + b2[1]*k1
    k2 = h * v(t1 + a[2]*h, 
               r + b2[1]*k1)
    #print 'THREE:', r + b3[1]*k1 + b3[2]*k2
    k3 = h * v(t1 + a[3]*h,
               r + b3[1]*k1 + b3[2]*k2)
    #print 'FOUR:', r + b4[1]*k1 + b4[2]*k2 + b4[3]*k3
    k4 = h * v(t1 + a[4]*h,
               r + b4[1]*k1 + b4[2]*k2 + b4[3]*k3)
    #print k5, r + b5[1]*k1 + b5[2]*k2 + b5[3]*k3 + b5[4]*k4
    k5 = h * v(t1 + a[5]*h,
               r + b5[1]*k1 + b5[2]*k2 + b5[3]*k3 + b5[4]*k4)
    #print k6, r + b6[1]*k1 + b6[2]*k2 + b6[3]*k3 + b6[4]*k4 + b6[5]*k5
    k6 = h * v(t1 + a[6]*h,
               r + b6[1]*k1 + b6[2]*k2 + b6[3]*k3 + b6[4]*k4 + b6[5]*k5)

    rnext  = r +  c[1]*k1 +  c[2]*k2 +  c[3]*k3 +  c[4]*k4 +  c[5]*k5 +  c[6]*k6
    rnexts = r + cs[1]*k1 + cs[2]*k2 + cs[3]*k3 + cs[4]*k4 + cs[5]*k5 + cs[6]*k6
    return rnext, rnext - rnexts
    
class BoundPrecision(object):
    '''
    A callable which will provide a relative value that keeps precision in bounds.
    '''

    def __init__(self, prec=0.001, maxrelval=None):
        '''
        Create a bound precision function.

        @param prec: desired precision in same units as error
        @type prec: float

        @param maxrelval: maximum allowed correction
        @type maxrelval: float

        '''
        self.prec = prec
        self.maxrat = maxrelval

    def __call__(self, error):
        '''
        Return a correction to apply to the a subsequent step size to
        bound the precision.

            bp = BoundPrecision(...)
            step *= bp(error)

        @param error: error, eg as vector displacement between two
            attempted steps to the same place.
        @type error: N-array
        @return: correction as ratio of current step size
        @rtype: float

        '''
        err = numpy.max(numpy.abs(error))
        if err == 0.0:
            return 1.0
        rel = self.prec/err
        rel = math.pow(rel, 0.2) # adaptive runge-kutta magic
        if self.maxrat is None:
            return rel
        return min(self.maxrat, rel)

class StopDetection(object):
    def __init__(self, distance=0.1, geom=None, nallowed=1):
        '''
        Create a StopDetection object.

        @param distance: the minimum distance a step is allowed to take before being considered potentially stuck.
        @type distance: float

        @param nallowed: the number of minimum steps in a row allowed.
        @type nallowed: int
        '''
        self.distance2 = distance**2
        self.nallowed = nallowed
        self.nstuck = 0
        self.geom = geom

    def __call__(self, r1, r2):
        # check geometry 
        if self.inside(r2):
            return True
        # are we stepping far enough?
        return self.stuck(r1,r2)
    
    def stuck(self, r1, r2):
        d2 = sum([(a-b)**2 for a,b in zip(r1,r2)])
        if d2 > self.distance2:
            self.nstuck = max(0, self.nstuck-1)
            return False
        self.nstuck += 1
        if self.nstuck <= self.nallowed:
            return False
        return True

    def inside(self, r):
        if self.geom is None:
            print 'Warning: No geometry defined.'
            return False
        for node in self.geom:
            if node.inside(r):
                if 'pixel25' not in node.name:
                    raise Exception
                return True

class CollectSteps(object): 
    '''
    A Stepper visitor that collects all steps and provides a bounds on the precision.
    '''

    Step = namedtuple('Step','t1 r1 v1 t2 r2 error')

    def __init__(self, stop=StopDetection(), bounds = BoundPrecision()):
        '''
        @param stop: a callable returning true if the stepping needs to stop
        @type stop: callable(t1,r1,v1,t2,r2)->bool

        @param bounds: a callable providing a scaling of dt
        @type bounds: callable(error)->float
        '''
        self.steps = list()
        self.bounds = bounds
        self.stop = stop

    def __call__(self, t1, r1, v1, t2, r2, error):
        '''
        Collect one step.

        @param t1: the time at the beginning of the step
        @type t1: float

        @param r1: the position at the beginning of the step
        @type r1: N-array

        @param v1: the speed at the beginning of the step
        @type v1: float

        @param t2: the time at the end of the step
        @type t2: float

        @param r2: the position at the end of the step
        @type r2: N-array

        @param error: and estimate of the error made in the step as a displacement vector
        @type error: N-array

        @return: a duration to use for the subsequent step
        @rtype: float

        @raise StopIteration: if step indicates the stepping is stuck or needs to stop

        '''
        self.steps.append(self.Step(t1,r1,v1,t2,r2,error))
        if self.stop(r1,r2):
            raise StopIteration
        return (t2-t1)*self.bounds(error)

    def __len__(self): return len(self.steps)

    @property
    def points(self):
        '''
        The N+1 step points.

        @type: list of N-array
        '''
        if not self.steps:
            return None
        p = [s.r1 for s in self.steps]
        p.append(self.steps[-1].r2)
        return p

    @property
    def times(self):
        '''
        The N+1 step times.

        @type: list of float
        '''
        if not self.steps:
            return None
        t = [s.t1 for s in self.steps]
        t.append(self.steps[-1].t2)
        return t

    @property
    def vel(self):
        '''
        The N+1 vel.

        @type: list of lists
        '''
        if not self.steps:
            return None
        v = [s.v1 for s in self.steps]
        v.append(self.steps[-1].v1)
        return v

    @property
    def distance(self):
        '''
        N distances between step k and k+1.

        @type: list of float
        '''
        ret = list()
        for s in self.steps:
            diff = s.r2 - s.r1
            dist = math.sqrt(sum([d**2 for d in diff]))
            ret.append(dist)
        return ret

    @property
    def array(self):
        '''
        Return array of shape (N+1,7).  Layout of last dimension is (x,y,z,vx,vy,vz,t).
        '''
        if not self.steps:
            return numpy.ndarray((0,5), dtype=float)
        p = numpy.asarray(self.points)
        t = numpy.asarray(self.times)
        v = numpy.asarray(self.vel)
        return numpy.vstack((p.T,v.T,t.T)).T
        
    def __str__(self):
        if not self.steps:
            return "0 steps"

        t = self.times
        mint = min(t)
        maxt = max(t)
        first = self.steps[0]
        last = self.steps[-1]
        return "%d steps from %f,%s --> %f,%s min:%f, max:%f" % \
            (len(self), first.t1, first.r1, last.t2, last.r2, mint, maxt)

def StepperParams(params, **defaults):
    '''
    Stepper parameters govern how steppers step.

    Not all keyword arguments are used by all steppings.

    @keyword maxiter: maximum number of steps to take (default:100)
    @type maxiter: int

    @keyword dt: initial step in time (default:None - calculate from
        velocity and lcar)
    @type dt: float [units of time]

    @keyword lcar: characteristic length used to set initial time step
        from initial velocity (default:None - rely on dt)
    @type lcar: float [units of distance]

    @return: a namedtuple of keywords
    @rtype: StepperParams namedtuple 
    '''
    p = dict(maxiter=100, dt = None, lcar = None,)
    p.update(**defaults)
    p.update(**params)
    SP = namedtuple('StepperParams','maxiter dt lcar')
    return SP(**p)

class Stepper(object):
    def __init__(self, velo_fun, step_fun = step_rkck, fixed_step=None, **kwds):
        '''
        Create a stepper.
        '''
        self.velo_fun = velo_fun
        self.step_fun = step_fun
        self.fixed_step = fixed_step
        self._defaults = kwds

    def params(self, **kwds):
        return StepperParams(kwds, **self._defaults)

    def __call__(self, time, position, visitor = CollectSteps(), **kwds):
        '''
        Step from starting time and position until end condition met.

        See StepperParams for keyword arguments.

        End conditions:

            - maximum number of iterations is reached,

            - visitor raises StopIteration.
        '''
        p = self.params(**kwds)
        position = numpy.asarray(position)

        vstart = self.velo_fun(time, position)
        vmag = math.sqrt(sum([v**2 for v in vstart]))
        dt = p.lcar / vmag
        if self.fixed_step:
                dt = self.fixed_step
        for count in range(p.maxiter):
            try:
                pnext, error = self.step_fun(position, time, time+dt, self.velo_fun)
            except ValueError:
                print 'Value Error'
                break
            #print 'Pos: (',position[0],',',position[1],',',position[2],')', '  time:', time, 'dt:', dt, '   ', pnext,'  ',error
        
            try:
                vel = self.velo_fun(time, position)
                dtnext = visitor(time, position, vel, time+dt, pnext, error)
            except StopIteration:
                print 'StopIteration'
                break
            if dtnext is None:
                print 'dtnext is None'
                break

            # update for next iteration
            time += dt  
            if self.fixed_step:
                dt = self.fixed_step
            else:
                dt = dtnext
            position = pnext
            #print '\tit =',count, 'position',position
        return visitor

def get_vertices(step_range, step_inc, stepfile):
    '''
    Return list of vertices based on input from stepfile or on a patch
    determined by step_range and step_inc.
    '''
    vtxs = list()

    # stepfile takes precedence 
    if stepfile is not None:
        try:
            f = open(stepfile)
        except:
            print 'Stepfile',stepfile,'does not exist!'

        # expecting list of name x y z 
        while True:
            linevec = f.readline().split()
            if len(linevec)<1:
                break
            name,x,y,z = str(linevec[0]), float(linevec[1]), float(linevec[2]), float(linevec[3])
            vtxs.append([name,(x,y,z)])
        f.close()
        return vtxs

    # otherwise create a patch
    xl, yl, zl = step_range
    xr, yr, zr = xl[1]-xl[0], yl[1]-yl[0], zl[1]-zl[0]
    xs, ys, zs = step_inc

    # we will take x as a constamt
    x = xl[0]
    count = 0
    for y in numpy.linspace(yl[0], yl[1], 1+int(yr/ys)):
        for z in numpy.linspace(zl[0], zl[1], 1+int(zr/zs)):
            name = 'path'+str(count)
            count += 1
            vtxs.append([name,(x,y,z)])
    return vtxs


def step(vfield,
         linspaces,
         geom,
         step_range = ( (1.5,1.5), (-1,1), (-1,1) ),
         step_inc   = (0.1, 0.1, 0.1),
         start_time = 0.0,
         lcar       = 0.01,
         stuck      = 0.01,
         stepfile   = None,
         stepper    = 'rkck', 
         **kwds):
    '''
    Return paths stepped through vfield starting at many points on a
    rectangular patch given by range = ( (xmin, xmax, step), ... ).  
    '''
    # define stepping function
    step_fun = eval("step_%s" % stepper) 

    # get the velocity field
    velo = Field(vfield, linspaces)
    def velocity(notused, r):
        return velo(r)

    # define stepper
    print 'Stepping with',step_fun
    stepper = Stepper(velocity, lcar=lcar, step_fun=step_fun, **kwds)

    # get step vertices
    assert(stepfile is not None)
    vtxs = get_vertices(step_range, step_inc, stepfile)
    
    from pixsim.models import Array
    paths, start = list(), list()
    for vtx in vtxs:
            name  = vtx[0]
            position = vtx[1]
            x,y,z = position
            start.append(position)
            visitor = stepper(start_time, position, CollectSteps(StopDetection(distance=stuck, geom=geom)))
            paths.append( Array(typename='tuples', name='path_from_'+name, data=visitor.array) )
            print 'Stepped',len(visitor.array),' times from (',x,',',y,',',z,') to (',visitor.array[-1][0],',',visitor.array[-1][1],',',visitor.array[-1][2],')'
    
    arr = Array(typename='points', name='vtxs', data=numpy.asarray(start))
    return [arr]+paths


