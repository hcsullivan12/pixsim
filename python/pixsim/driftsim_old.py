import ROOT
import pixsim.step as step
import pixsim.vector as vector
import numpy as np
import random as random
import time
import pixsim.ds as ds

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

class Batch(object):
    """Structure to keep track of tips and ids 
    for a batch of electron clouds from ntuple."""
    def __init__(self, tips, ids):
        self.tips = tips
        self.ids = ids

def check_next(did_stop, exclude, tips, next_points, stopper):
    """Given the current list of steps and the next step,
    check if we should stop here or keep stepping."""
    for pt in range(0, len(next_points)):
        if did_stop[pt]:
            continue
        the_tip_pt  = tips[pt].head.data
        the_next_pt = next_points[pt]
        try:
            didit = stopper(the_tip_pt[0:3],the_next_pt[0:3])
            if didit:
                did_stop[pt] = True
        except:
            exclude[pt] = True

        tips[pt].insert(the_next_pt)
    return did_stop, exclude, tips

def do_batch_stepping(tips, stepper=None, stopper=None, maxiter=300, fixed_step=0.05, **kwds):
    """Method to apply stepping procedure on
    batch of points."""
    if len(tips) == 0:
        return
    # array to keep track of stopped paths and paths that do not make it to center pixel
    did_stop = [False for n in range(0,len(tips))]
    exclude = [False for n in range(0,len(tips))]
    for istep in range(maxiter):
        # do the next step
        points = np.asarray([x.head.data for x in tips])
        next_points = stepper(fixed_step, points)
        assert(len(points) == len(next_points))
        p1 = points[0]
        p2 = next_points[0]
        diff = p1-p2
        #print (diff[0]**2+diff[1]**2+diff[2]**2)**0.5
        
        # check for collision with pixels
        did_stop, exclude, tips = check_next(did_stop, exclude, tips, next_points, stopper)
    return tips, exclude

def round_down(x, fixed_step):
    return np.floor(x / fixed_step) * fixed_step

def round_near(x, base):
    return (base * np.round(x/base)).astype(int)

def do_batch_induction(my_data, exclude, entry, tips, ids, vfield, wfield, linspaces, 
                       fixed_step=0.05,
                       **kwds):
    """This method will calculate the induced current
    on a pixel due to each electron cloud in the batch.

    The wfield is taken to be centered on the pixel that
    collected the cloud."""

    velo = vector.Field(vfield, linspaces)
    def velocity(points):
        return velo.do_many(points)

    wght = vector.Field(wfield, linspaces)
    def weight(points):
        return wght.do_many(points)

    base = int(str(fixed_step)[::-1].find('.'))
    tick_conversion = 10**base

    for eid, tip, excl in zip(ids, tips, exclude):
        if excl:
            continue
        pid = entry.ides_voxel_ch[eid]
        el  = entry.ides_numElectrons[eid]
        # convert electrons to fC
        el *= 1.6 / 10000

        # compute waveform
        # velocity --> cm/us
        # charge --> fC
        # current = q * W * v = fC * 1/cm * cm/us = nA
        xyzts = tip.array()
        ws = weight(xyzts[:,0:3])
        vs = velocity(xyzts[:,0:3])
        wvf = np.einsum('ij,ij->i', vs, ws)
        # converting ticks to integers
        # this is to fix edge cases where a tick is not rounding correctly
        new_ticks = round_near(tick_conversion*xyzts[:,-1], base)

        # clip the end if we have duplicate
        if new_ticks[-1] == new_ticks[-2]:
            rm = new_ticks.pop()
            wvf = wvf[:-1]

        wvf *= -1 * el / abs(fixed_step*np.sum(wvf)) # normalize to numElectrons
        # add entry 
        if pid in my_data.keys():
            old_wvf = my_data[pid]
            for k,v in zip(new_ticks,wvf):
                if k in old_wvf.keys():
                    old_wvf[k] += v
                else:
                    old_wvf[k] = v    
            my_data[pid] = old_wvf
        else:
            new_wvf = {k:v for k,v in zip(new_ticks,wvf)}
            my_data[pid] = new_wvf
    return my_data

def make_rtds(my_data,
              schmitt_time=0.02, 
              threshold=0.1, 
              gain=1.0): 
    """Integrate waveforms directly and return RTDs."""

    pids = [p for p in my_data.keys()]
    base = int(str(fixed_step)[::-1].find('.'))
    tick_conversion = 10**base
    for p in pids:
        wvf = my_data[p]
        ts = sorted(wvf)
        ys = [wvf[t] for t in ts]
        ts = [1.0*t/tick_conversion for t in ts] # converting back to float (us)
        # appending extra zeros
        step_size = ts[1]-ts[0]
        to_append = [s for s in np.arange(ts[-1], ts[-1]+0.5,step_size)]
        ts = np.append(ts, to_append)
        ys = np.append(ys, [0]*len(to_append)) 

        # integrate
        sch_data = {t:0 for t in np.arange(ts[0], ts[-1], schmitt_time)}
        ana_data = sch_data
        charge = 0
        for tick,amp in zip(ts,ys):
            charge += amp*step_size # nA * us = fC
            if charge >= threshold:
                charge = 0
                sch_data[tick] = 0.5
            ana_data[tick] = charge

        sch_ts = sorted(sch_data)
        sch_ys = [sch_data[t] for t in sch_ts]
        ana_ys = [ana_data[t] for t in sch_ts]
        return sch_ts, sch_ys, ana_ys

def make_waveforms(my_data,
                   schmitt_time=0.02,
                   threshold=0.5, # fC
                   gain=10., # V/fC
                   fixed_step=0.05,
                   noise=300, # electrons
                   **kwds):
    """Combine the contents of the data to
    form single waveforms for each pixel."""
    import matplotlib.pyplot as plt
    pids = [p for p in my_data.keys()]
    base = int(str(fixed_step)[::-1].find('.'))
    tick_conversion = 10**base

    # convert electrons to fC
    noise *= 1.6 / 10000

    for p in pids:
        wvf = my_data[p]
        ts = sorted(wvf)
        ys = [wvf[t] + np.random.normal(0, noise) for t in ts]

        # appending extra zeros
        step_size = ts[1]-ts[0]
        to_append = [s for s in np.arange(ts[-1]+step_size, ts[-1]+50,step_size)]
        ts = np.append(ts, to_append)
        ys = np.append(ys, np.random.normal(0, noise, len(to_append))) 

        # integrate
        sch_data = {t:0 for t in np.arange(ts[0], ts[-1], 100*schmitt_time)}
        sch_hits = list()
        ana_data = {t:0 for t in ts}
        charge = 0
        for tick,amp in zip(ts,ys):
            charge += amp*fixed_step # nA * us = fC
            ana_data[tick] = charge*gain
            if charge >= threshold:
                charge = 0
                sch_data[tick] = 2.0
                sch_hits.append(tick)

        sch_ts = sorted(sch_data)
        ana_ts = sorted(ana_data)
        sch_hits.sort()
        sch_ys = [sch_data[t] for t in sch_ts]
        ana_ys = [ana_data[t] for t in ana_ts]

        # reconstructing
        reco_data = {t:0 for t in ts}
        for tick in range(1, len(sch_hits)):
            pre_time = sch_hits[tick-1]
            cur_time = sch_hits[tick]
            delta_t = cur_time-pre_time
            for tbin in np.arange(pre_time, cur_time, step_size):
                reco_data[tbin] = tick_conversion * threshold / delta_t 

        reco_ys = [reco_data[t] for t in ts]

        # converting back to float (us)
        ts     = [1.0*t/tick_conversion for t in ts]
        sch_ts = [1.0*t/tick_conversion for t in sch_ts]
        ana_ts = [1.0*t/tick_conversion for t in ana_ts]
        
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(15, 15))

        ax1.step(ts, ys)
        ax1.step(ts, reco_ys)
        ax2.step(sch_ts, sch_ys)
        ax2.step(ana_ts, ana_ys)

        fig.subplots_adjust(hspace=0)
        plt.xlabel('Time [us]',fontsize=15)
        ax1.set_ylabel('Current [nA]',fontsize=15)
        ax2.set_ylabel('Voltage [V]',fontsize=15)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.show()

def drift(vfield, wfield, points, linspaces, geom, entry, 
          backstep = 0.5, 
          stuck    = 0.01,
          drift_speed = 0.163, # cm/us
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
    way and manifests itself as a shift in x or y and z.
    """

    # method to grab a batch of points
    def get_batch(entry_list):
        smr_ys, smr_zs = np.asarray([entry.ides_y[i] for i in entry_list]), np.asarray([entry.ides_z[i] for i in entry_list])
        pix_ys, pix_zs = np.asarray([entry.ides_voxel_y[i] for i in entry_list]), np.asarray([entry.ides_voxel_z[i] for i in entry_list])
        # NASTY BUG in larg4 nearest channel
        for count in range(0,len(smr_ys)):
            y_diff = smr_ys[count] - pix_ys[count]
            z_diff = smr_zs[count] - pix_zs[count]
            if abs(y_diff) > 1.0 or abs(z_diff) > 1.0:
                print '\nWARNING: Nasty Bug',y_diff,z_diff
                smr_ys[count] = pix_ys[count]
                smr_zs[count] = pix_zs[count]

        rel_ys, rel_zs = np.asarray([smr_ys - pix_ys]).T, np.asarray([smr_zs - pix_zs]).T
        smr_xs = np.asarray([backstep*np.ones(len(smr_ys))]).T
        smr_ts = np.asarray([[round_down(entry.ides_voxel_x[i]/drift_speed, kwds['fixed_step']) for i in entry_list]]).T

        rel_pos = np.concatenate( (smr_xs, rel_ys, rel_zs, smr_ts), axis=1 )
        rel_pos = [ds.LinkedList(ds.Node(data=x)) for x in rel_pos]
        batch = Batch(rel_pos, entry_list)
        return batch

    # get the velocity field
    velo = vector.Field(vfield, linspaces)
    def velocity(points):
        return velo.do_many(points)

    stepper = BatchedStepper_rkck(velocity)
    stopper = step.StopDetection(geom=geom, distance=stuck)

    # grabbing batches of points
    import sys
    my_data = dict()
    entries = [x for x in range(0,len(entry.ides_x))]
    random.shuffle(entries)

    batch_size = 100
    total_batches = int(np.ceil(len(entries)/batch_size))+1
    for bid,index in enumerate(range(0,len(entries), batch_size),1):
        bck = index+batch_size
        if bck >= len(entries):
            bck = len(entries)-1

        msg = 'Batch = '+str(bid)+' / '+str(total_batches)
        sys.stdout.write("\r%s" % msg)
        sys.stdout.flush()

        use_this_batch = [entries[i] for i in range(index,bck+1)]
        batch = get_batch(use_this_batch)

        # do the stepping
        batch.tips, exclude = do_batch_stepping(batch.tips, stepper=stepper, stopper=stopper, **kwds)
        # do the induction
        my_data = do_batch_induction(my_data, exclude, entry, batch.tips, batch.ids, vfield, wfield, linspaces, **kwds)
    print ''
    print len(my_data),'pixels collected charge'
    make_waveforms(my_data, **kwds)
    return my_data
     
def analyze(my_data):
    """Area to study waveforms."""
    return

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
        my_data = drift(vel, wfield, points, linspaces, geom, entry, **kwds)
        analyze(my_data)
        break



    