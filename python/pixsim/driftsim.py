import ROOT
import numpy as np
import sympy as sp
from bisect import bisect_left
from sets import Set
import sys

# default settings, will detect spacing later
_pixel_spacing = 0.4
_eps = np.round(np.sqrt(2*(0.5*_pixel_spacing)**2), 2)
_step_zs, _step_ys = Set(), dict()

_hDiffusion1 = ROOT.TH2D('_hDiffusion1', '_hDiffusion1', 400, 0, 400, 20, 0, 20)
_hDiffusion2 = ROOT.TH2D('_hDiffusion2', '_hDiffusion2', 50, 0, 50, 400, 0, 400)

def is_accurate(dist, eps=_eps):
    """Check against largest distance from pixel center"""
    return dist < eps

def distance2d(pos1, pos2):
    """Calculate 2D distance"""
    return np.sqrt((pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)

def get_nearest(lpos, pos):
    """
    Assumes lpos is sorted. Returns closest value to pos.
    """
    loc = bisect_left(lpos, pos)
    if loc == len(lpos):
        loc = loc-1
    elif loc != 0:
        before = lpos[loc - 1]
        after = lpos[loc]
        if after - pos > pos - before:
           loc = loc-1
    return loc

def nearest_pixel(pos, pixels_y, pixels_z, eps=_eps):
    """
    Find nearest pixel ID. Pixels are numbered
    in increasing order from top-left to right
    and down.

    With this scheme, since the arrays are assumed 
    to be sorted, convert y accordingly. 

    Returns a number in range [0,1,2,...,n-1]
    """
    row = get_nearest(pixels_y, pos[1])
    col = get_nearest(pixels_z, pos[2])
    near_pos = [pos[0], pixels_y[row], pixels_z[col]]
    row = len(pixels_y) - get_nearest(pixels_y, pos[1]) - 1

    dist = distance2d(pos, near_pos) #[(x-y)**2 for x,y in zip(pos, near_pos)]
    #dist = sum(dist)**0.5
    if not is_accurate(dist, eps=eps):
        raise ValueError('Nearest pixel is further than epsilon='+str(eps)+' away ('+str(dist)+').')
    return row * len(pixels_z) + col

def pixel_position(pid, pixels_y, pixels_z):
    """
    Get pixel position from pixel ID.
    Remember, pixels_y is in reverse id 
    order.
    """
    upper_lim = len(pixels_y)*len(pixels_z)-1
    if pid > upper_lim or pid < 0:
        raise ValueError('Pixel ID='+str(pid)+' out of range=[0,'+str(upper_lim)+']')
    row = pid/len(pixels_z)
    col = pid - row * len(pixels_z)
    return [0.0, pixels_y[len(pixels_y)-1-row], pixels_z[col]]

def get_response(rel_pos, rsp_data):
    """Find nearest response function"""
    mxm = float("inf")
    ret_fn, ret_ep = None, None
    for z,ymap in rsp_data.iteritems():
        for y,sat in ymap.iteritems():
            dist = distance2d(rel_pos, [0,y,z])
            if dist < mxm:
                mxm = dist
                ret_fn, ret_ep = {'rsp':sat['rsp'], 'area':sat['area'], 'max':sat['max']}, sat['ep']
    return ret_fn, ret_ep

def get_rotation(x, y):
    """Find rotation necessary to get to 1st quadrant"""
    if x >= 0 and y >= 0:
        return 0
    elif x < 0 and y >= 0:
        return 270
    elif x < 0 and y < 0:
        return 180
    else:
        return 90

def rotate(x, y, theta):
    """Rotate coordinates"""
    c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    return c*x - s*y, s*x + c*y

def check_bipolar(pos0, pixels_y, pixels_z, pix_pos, rsp_data):
    """
    Method to extract information about drift path.
    """
    # find rotation to get to 1st quadrant
    rel_pos = [x-y for x, y in zip(pos0, pix_pos)]
    rot_ang = get_rotation(rel_pos[2], rel_pos[1])
    rot_rel_pos = list(rel_pos)
    rot_rel_pos[2], rot_rel_pos[1] = rotate(rel_pos[2], rel_pos[1], rot_ang)
    if abs(rot_rel_pos[1]) < 0.0001:
        rot_rel_pos[1] = 0.0
    if abs(rot_rel_pos[2]) < 0.0001:
        rot_rel_pos[2] = 0.0
    assert rot_rel_pos[1] >= 0 and rot_rel_pos[2] >= 0, 'Relative positions = '+str(rot_rel_pos[1])+' '+str(rot_rel_pos[2])
    
    # flip over y=z if necessary
    did_flip = False
    if rot_rel_pos[1] > rot_rel_pos[2]:
        rot_rel_pos[1], rot_rel_pos[2] = rot_rel_pos[2], rot_rel_pos[1]
        did_flip = True
    
    # find nearest response function and endpoints
    resp, rel_ep = get_response(rot_rel_pos, rsp_data)
    
    # convert ep back to global ep
    ep = list(rel_ep)
    if did_flip:
        ep[1], ep[2] = ep[2], ep[1]
    ep[2], ep[1] = rotate(ep[2], ep[1], -1*rot_ang)
    ep = [x+y for x, y in zip(ep, pix_pos)]
    
    dist = distance2d(ep, pix_pos)
    is_bip = False
    if not is_accurate(dist):
        is_bip = True

    return is_bip, resp, ep, dist

def integrate_waveform(ts, ys, threshold=0.7, sch_time=0.02):
    """Integrate waveform and handle resets"""

    ys_int, ys_sch, sch_hits = list(ys), list(ys), list()
    charge = 0 

    for tbin in range(0, len(ys_int)):
        charge += ys[tbin] * sch_time
        if charge > threshold:
            charge = 0
            ys_sch[tbin] = 1.0
            sch_hits.append(ts[tbin])
        ys_int[tbin] = charge
    
    return ys_int, ys_sch, sch_hits

def reconstruct(ts, sch_hits, threshold=0.7, sch_time=0.02):
    """Reconstruct waveform from resets"""

    reco_data = {t:0 for t in ts}
    for tbin in range(1, len(sch_hits)):
        pre_time = sch_hits[tbin-1] 
        cur_time = sch_hits[tbin] 
        delta_t = cur_time-pre_time
        for tick in np.arange(pre_time, cur_time, int(100*sch_time)):
            reco_data[tick] = 100. * threshold / delta_t
    
    return reco_data

def analyze(wvfs):
    """Area for analysis"""
    import matplotlib.pyplot as plt
    global _hDiffusion1, _hDiffusion2
    
    for pid, wvf in wvfs.iteritems():
        ts = wvf.keys()
        ts.sort()
        ys = [wvf[t] for t in ts]

        ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys)
        if len(sch_hits) < 2:
            continue
        
        reco_data = reconstruct(ts, sch_hits) 
        reco_ys = [reco_data[t] for t in ts]

        sch_hits.sort()

        #for hit in range(1, len(sch_hits)):
        #    tdiff = (sch_hits[hit] - sch_hits[hit-1]) / 100
        #    xpos = 0.1648 * (sch_hits[hit] / 100)
        #    _hDiffusion2.Fill(tdiff, xpos)

        ts = [t/100. for t in ts]
        thst = ROOT.TH1F('thst'+str(pid), 'thst'+(pid), len(ts), ts[0], ts[-1])
        max_amp = -1
        max_amp_bin = -1

        for hbin in range(0, len(reco_ys)):
            thst.SetBinContent(hbin+1, ys[hbin])
            if reco_ys[hbin] > max_amp:
                max_amp = reco_ys[hbin]
                max_amp_bin = hbin

        time_diff = (sch_hits[-1] - sch_hits[0]) / 100
        fit = ROOT.TF1('fit'+str(pid), 'gaus', ts[0], ts[-1])
        fit.SetParameter(0, max_amp)
        fit.SetParameter(1, ts[max_amp_bin])
        fit.SetParameter(2, 0.25*time_diff)

        thst.Fit(fit, 'QRN')
        _hDiffusion1.Fill(fit.GetParameter(1), np.power(fit.GetParameter(2), 2.))
        #plt.step(ts, ys)
        #plt.step(ts, reco_ys)
        #plt.plot(ts, fit.GetParameter(0)*np.exp(-np.power(np.asarray(ts) - fit.GetParameter(1), 2.) / (2 * np.power(fit.GetParameter(2), 2.))))
        #plt.show()

def best_fit(X, Y):
    """Not used"""

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    
    b = numer / denum
    a = ybar - b * xbar
    
    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))
    
    return a, b
        
    a, b = best_fit(_xs, _ys)

    yfit = [a + b * xi for xi in _xs]
    plt.scatter(_xs, _ys)
    plt.plot(_xs, yfit)
    plt.show()

def plot(wvfs):
    import matplotlib.pyplot as plt
    for pid, wvf in wvfs.iteritems():
        ts = wvf.keys()
        ts.sort()
        ys = [wvf[t] for t in ts]

        ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys) 
        
        reco_data = reconstruct(ts, sch_hits) 
        reco_ys = [reco_data[t] for t in ts]

        fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(15, 15))
        ts = [t/100. for t in ts]
        ax2.step(ts, ys)
        ax2.step(ts, reco_ys)
        ax1.step(ts, ys_int)
        ax1.step(ts, ys_sch)
        plt.xlabel('Time [us]',fontsize=15)
        ax2.set_ylabel('Current [nA]',fontsize=15)
        ax1.set_ylabel('Voltage [V]',fontsize=15)
        #ax2.set_ylim([0,3.5])
        plt.show()

def fill_map(entry, pixels_y, pixels_z, rsp_data):
    """
    Fill pixel to waveform map.

    We begin by guessing the nearest pixel to the drifted yz position.
    Extract the drift information and if the result is bipolar, use
    the end point to find the correct pixel, extract again.
    """

    pid_to_wvf = dict()
    chg_entries = len(entry.ides_y)
    for chg in range(0, chg_entries):
        # getting charge position
        pos0 = np.asarray([0.0, entry.ides_y[chg], entry.ides_z[chg]])
        time_tick = entry.ides_voxel_x[chg]/0.1648 # convert to us
        numEl = entry.ides_numElectrons[chg]
        
        # we will first take nearest pixel
        # if the result is bipolar, redo using correct pixel
        # ignoring charge that missed the pixel plane
        try:
            near_pix_id = nearest_pixel(pos0, pixels_y, pixels_z)
        except:
            continue
        pix_pos = pixel_position(near_pix_id, pixels_y, pixels_z)
        is_bip, resp, ep, dist = check_bipolar(pos0, pixels_y, pixels_z, pix_pos, rsp_data)

        if is_bip:
            try:
                near_pix_id = nearest_pixel(ep, pixels_y, pixels_z)
            except:
                continue
            pix_pos = pixel_position(near_pix_id, pixels_y, pixels_z)
            is_bip, resp, ep, dist = check_bipolar(pos0, pixels_y, pixels_z, pix_pos, rsp_data)

        if is_bip:
            print 'Warning: Could not determine proper pixel ID!'
            continue

        # we have pid, unipolar rsp, chg tick
        to_store = {'el':numEl, 'rsp':resp['rsp'], 'area':resp['area'], 't':time_tick, 'max':resp['max']} 
        if near_pix_id not in pid_to_wvf.keys():
            pid_to_wvf[near_pix_id] = [to_store]
        else:
            pid_to_wvf[near_pix_id].append(to_store)
    return pid_to_wvf

def convolve(pid_to_wvf, rsp_time_bin=0.02):
    """Perform convolution of charge deposits and response function"""

    to_fC_convers = 1.6 / 10000

    pid_count = 0
    tot_pid_count = len(pid_to_wvf)
    for pid, dlist in pid_to_wvf.iteritems():
        pid_count += 1
        msg = 'Convolving '+str(pid_count)+' / '+str(tot_pid_count)
        sys.stdout.write("\r%s" % msg)
        sys.stdout.flush()

        dlist.sort(key=lambda x: x['t'])

        for store in dlist:
            rsp_area = store['area'] * rsp_time_bin
            chg_in_fC = store['el'] * to_fC_convers
            store['rsp_convolved'] = store['rsp'] * chg_in_fC / rsp_area

    print ''
    return pid_to_wvf

def make_waveforms(pid_to_wvf, rsp_time_bin=0.02):
    """Construct waveforms by taking superposition of responses"""

    base = int(100 * rsp_time_bin)
    def round_it(x, base=base):
        return int(base * np.floor(float(x)/base))

    wvfs = dict()

    pid_count = 0
    tot_pid_count = len(pid_to_wvf)
    for pid, dlist in pid_to_wvf.iteritems():
        pid_count += 1
        msg = 'Constructing '+str(pid_count)+' / '+str(tot_pid_count)
        sys.stdout.write("\r%s" % msg)
        sys.stdout.flush()

        wvf = dict()
        
        for store in dlist:
            # start time
            tick = 100 * (store['t'] - rsp_time_bin * (store['max'] + 1))
            tick = round_it(tick)
            
            for tbin in range(0, len(store['rsp_convolved'])):
                if tick not in wvf:
                    wvf[tick] = store['rsp_convolved'][tbin]
                else:
                    wvf[tick] += store['rsp_convolved'][tbin]
                tick += 100 * rsp_time_bin 
                tick = round_it(tick)

        wvfs[pid] = wvf

    print ''
    return wvfs

def calculate(entry, pixels_y, pixels_z, rsp_data):
    """
    Calculate waveforms for each charge cloud. 
    Convolves charge deposits with field/electronics
    response function.
    """
    pid_to_wvf = fill_map(entry, pixels_y, pixels_z, rsp_data)
    pid_to_wvf = convolve(pid_to_wvf)
    wvfs = make_waveforms(pid_to_wvf)
    #plot(wvfs)
    analyze(wvfs)

def get_pixels(entry):
    """Get pixel coordinates from ntuple"""
    pixels_y = [y for y in entry.pixel_y]
    pixels_z = [z for z in entry.pixel_z]
    pixels_y.sort()
    pixels_z.sort()

    global _pixel_spacing, _eps
    _pixel_spacing = np.round(abs(pixels_y[1]-pixels_y[0]), 2)
    _eps = np.round(np.sqrt(2*(0.5*_pixel_spacing)**2), 2)
    return pixels_y, pixels_z

def prepare_data(response):
    """Prepare data structure from responses"""
    rsp_sp, rsp_ep, rsp_fn = None, None, None
    for arr in response.data:
        if arr.name == 'start_points':
            rsp_sp = arr.data
        elif arr.name == 'end_points':
            rsp_ep = arr.data
        elif arr.name == 'responses':
            rsp_fn = arr.data
    assert rsp_sp is not None and rsp_ep is not None and rsp_fn is not None

    rsp_data = dict()
    for sp, ep, fn in zip(rsp_sp, rsp_ep, rsp_fn):
        sy, sz = sp[1], sp[2]
        sat_data = {'ep':ep, 'rsp':fn, 'area':np.sum(fn), 'max':np.argmax(fn)}

        if sz not in rsp_data.keys():
            rsp_data[sz] = {sy:sat_data}
        elif sy not in rsp_data[sz].keys():
            rsp_data[sz][sy] = sat_data
        else:
            raise ValueError('Data already contains yz pair!')
    return rsp_data

def sim(response, ntuple=None, **kwds):
    """
    This method will loop over the events in an ntuple.
    The ntuple, created with larg4, is expected to have
    the following information: 
        1) 3D coordinates of charge clouds close to 
           pixel plane.
        2) 2D coordinates of pixels
    """

    rfile = ROOT.TFile.Open(ntuple)
    rtree = rfile.Get('myana/anatree')

    # prepare response data
    rsp_data = prepare_data(response)

    # Using first entry to fill pixel coordinates
    pixels_y, pixels_z = None, None
    for entry in rtree:
        pixels_y, pixels_z = get_pixels(entry)
        break
    
    for entry in rtree:
        print 'drifting event', entry.event
        calculate(entry, pixels_y, pixels_z, rsp_data)
        #break

    _canvas1 = ROOT.TCanvas('_canvas1', '_canvas1', 800, 800)
    _hDiffusion1.Draw('colz')
    
    wait = input()