import ROOT
import pixsim.step as step
import pixsim.vector as vector
import numpy as np
import time

class Node(object):
    """Simple node structure."""
    def __init__(self, data=None, next_node=None):
        self.data = data
        self.next_node = next_node

    def set_next(self, new_next):
        self.next_node = new_next

class LinkedList(object):
    """Simple linked list structure."""
    def __init__(self, head=None):
        self.head = head
        self.size = 0
        if self.head is not None:
            self.size = 1

    def insert(self, data):
        node = Node(data=data, next_node=self.head)
        self.head = node
        self.size += 1

    def __len__(self):
        return self.size

class BatchedStepper_rkck(object):
    """Create a batched RK-C/K stepper using 
        velocity function."""
    def __init__(self, vfield):
        self.velo = vfield
        # The Cash/Karp coefficients.  Use None placeholders to match notation
        self.a  = [None, None, 0.2, 0.3, 0.6, 1.0, 0.875 ]
        self.cs = [None, 2825.0/27648.0, 0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25]
        self.c  = [None, 37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0 ]
        self.b2 = [None, 0.2]
        self.b3 = [None, 3.0/40.0, 9.0/40.0]
        self.b4 = [None, 0.3, -0.9, 1.2]
        self.b5 = [None, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0]
        self.b6 = [None, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]

    def velo_at(self, points, deltat=0.0, deltar=np.zeros(3).T):
        dr = points[:,:3]
        points[:,:3] += deltar
        points[:,3] += deltat
        return self.velo(points[:,:3])

    def __call__(self, dt, points):
        pcp = np.copy(points)
        k1 = dt * self.velo_at(pcp)
        k2 = dt * self.velo_at(pcp, self.a[2]*dt, self.b2[1]*k1)
        k3 = dt * self.velo_at(pcp, self.a[3]*dt, self.b3[1]*k1 + self.b3[2]*k2)
        k4 = dt * self.velo_at(pcp, self.a[4]*dt, self.b4[1]*k1 + self.b4[2]*k2 + self.b4[3]*k3)
        k5 = dt * self.velo_at(pcp, self.a[5]*dt, self.b5[1]*k1 + self.b5[2]*k2 + self.b5[3]*k3 + self.b5[4]*k4)
        k6 = dt * self.velo_at(pcp, self.a[6]*dt, self.b6[1]*k1 + self.b6[2]*k2 + self.b6[3]*k3 + self.b6[4]*k4 + self.b6[5]*k5)

        rnext = np.copy(points)
        rnext[:,:3] += self.c[1]*k1 + self.c[2]*k2 + self.c[3]*k3 + self.c[4]*k4 + self.c[5]*k5 + self.c[6]*k6
        rnext[:,3] += dt
        return rnext

def check_next(did_stop, tips, next_points, stopper):
    for pt in range(0, len(next_points)):
        if did_stop[pt]:
            continue
        the_tip_pt  = tips[pt].head.data
        the_next_pt = next_points[pt]
        if stopper(the_tip_pt[0:3],the_next_pt[0:3]):
            did_stop[pt] = True
        tips[pt].insert(the_next_pt)
    return did_stop, tips

def do_batch_stepping(tips, stepper=None, stopper=None, maxiter=100, fixed_step=0.1, **kwds):
    """Method to apply stepping procedure on
    batch of points."""
    if len(tips) == 0:
        return
    # array to keep track of stopped paths
    did_stop = [False for n in range(0,len(tips))]
    for istep in range(maxiter):
        #print 'Step =',istep
        # do the next step
        points = np.asarray([x.head.data for x in tips])
        next_points = stepper(fixed_step, points)
        assert(len(points) == len(next_points))
        ## check for collision with pixels
        did_stop, tips = check_next(did_stop, tips, next_points, stopper)

def drift(vel, wfield, points, linspaces, geom, entry, 
          backstep = 1.0, 
          stuck    = 0.01,
          **kwds):
    """This method will step the clouds from an initial point
    through the electric field, collect on nearest pixel, and
    calculate induced current on pixel(s).
    
    Taking the path of least resistance, it is assumed that the
    points used are the positions of the charge clouds at the 
    readout plane, so in some sense we've already messed up.
    But we will keep the yz positions and step back to 
    x ~= backstep cm to run the stepping procedure. Either way, 
    I think the error here affects all clouds in the same 
    way and manifests itself as a shift in x or y and z."""

    # method to grab a batch of points
    def get_batch(pt, bck):
        smr_ys, smr_zs = np.asarray(entry.ides_y[pt:bck]), np.asarray(entry.ides_z[pt:bck])
        pix_ys, pix_zs = np.asarray(entry.ides_voxel_y[pt:bck]), np.asarray(entry.ides_voxel_z[pt:bck])
        rel_ys, rel_zs = np.asarray([smr_ys - pix_ys]).T, np.asarray([smr_zs - pix_zs]).T
        smr_xs = np.asarray([backstep*np.ones(len(smr_ys))]).T
        smr_ts = np.asarray([[x/0.1648 for x in entry.ides_voxel_x[pt:bck]]]).T
        rel_pos = np.concatenate( (smr_xs, rel_ys, rel_zs, smr_ts), axis=1 )
        rel_pos = [LinkedList(Node(data=x)) for x in rel_pos]
        return rel_pos

    # get the velocity field
    velo = vector.Field(vel, linspaces)
    def velocity(points):
        return velo.do_many_old(points)

    stepper = BatchedStepper_rkck(velocity)
    stopper = step.StopDetection(geom=geom, distance=stuck)

    # grabbing batches of points
    import sys
    batch_size = 100
    for batch,pt in enumerate(range(0,len(entry.ides_x), batch_size),1):
        bck = pt+batch_size
        if bck >= len(entry.ides_x):
            bck = -1 

        msg = 'Batch = '+str(batch)+' / '+str(len(entry.ides_x)/batch_size)
        sys.stdout.write("\r%s" % msg)
        sys.stdout.flush()
        true_start_xs  = np.asarray(entry.ides_x[pt:bck])
        rel_pos = get_batch(pt, bck)
        do_batch_stepping(rel_pos, stepper=stepper, stopper=stopper, **kwds)
    print '\n'

def sim(vel, wfield, points, linspaces, geom, ntuple=None, **kwds):
    """This method will loop over the events in an ntuple.
    The ntuple, created with larg4, is expected to have
    the following information: 
        1) 3D coordinates of charge clouds close to 
           pixel plane.
        2) 2D coordinates of pixels"""

    rfile = ROOT.TFile.Open(ntuple)
    rtree = rfile.Get('myana/anatree')
    for entry in rtree:
        print 'drifting event',entry.event
        drift(vel, wfield, points, linspaces, geom, entry, **kwds)



    