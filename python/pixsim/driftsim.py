import ROOT
import numpy as np
import sympy as sp
from bisect import bisect_left

# default settings, will detect spacing later
_pixel_spacing = 0.4
_eps = np.round(np.sqrt(2*(0.5*_pixel_spacing)**2), 2)

def is_accurate(dist, eps=_eps):
    """Check against largest distance from pixel center"""
    #print _eps
    return dist < eps

def distance2d(pos1, pos2):
    """Calculate 2D distance"""
    return np.sqrt((pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)

def get_nearest(pixels_pos, pos):
    """
    Assumes pixels_pos is sorted. Returns closest value to pos.
    """
    loc = bisect_left(pixels_pos, pos)
    if loc == len(pixels_pos):
        loc = loc-1
    elif loc != 0:
        before = pixels_pos[loc - 1]
        after = pixels_pos[loc]
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

def get_response(rel_pos, rsp_sp, rsp_ep, rsp_fn):
    """Find nearest response function"""
    mxm = float("inf")
    nid = -1
    for count, sp in enumerate(rsp_sp):
        dist = distance2d(rel_pos, sp)  #((rel_pos[1]-sp[1])**2 + (rel_pos[2]-sp[2])**2)**0.5
        if dist < mxm:
            mxm = dist
            nid = count
    return rsp_fn[nid], rsp_sp[nid], rsp_ep[nid]

def calculate(entry, pixels_y, pixels_z, rsp_sp, rsp_ep, rsp_fn, backstep=1.0, **kwds):
    """
    Calculate waveforms for each charge cloud. 
    Convolves charge deposits with field/electronics
    response function.

    Finds the nearest pixel to the drifted yz position 
    and uses the relative position to that pixel to 
    lookup a response function. 
    """
    chg_entries = len(entry.ides_y)
    for chg in range(0, chg_entries):
        # getting charge location and pixel info
        pos0 = np.asarray([backstep, entry.ides_y[chg], entry.ides_z[chg]])
        # ignoring charge that missed the pixel plane
        try:
            near_pix_id = nearest_pixel(pos0, pixels_y, pixels_z)
        except:
            continue
        pix_pos = pixel_position(near_pix_id, pixels_y, pixels_z)

        # convert relative pos to 1st quad and reflect over y=z if needed
        did_flip = False
        transf = [np.sign(x-y) for x, y in zip(pos0, pix_pos)]
        rel_pos_abs = [abs(x-y) for x, y in zip(pos0, pix_pos)]

        if rel_pos_abs[1] > rel_pos_abs[2]:
            did_flip = True
            temp = rel_pos_abs[1]
            rel_pos_abs[1] = rel_pos_abs[2]
            rel_pos_abs[2] = temp

        # find nearest response function endpoints
        resp, rel_sp, rel_ep = get_response(rel_pos_abs, rsp_sp, rsp_ep, rsp_fn)

        # convert ep back to global ep
        ep = rel_ep.copy()
        if did_flip:
            temp = ep[1]
            ep[1] = ep[2]
            ep[2] = temp
        ep = [t*x+y for t, x, y in zip(transf, rel_ep, pix_pos)]

        dist = distance2d(ep, pix_pos)
        #print dist
        is_bip = False
        if not is_accurate(dist):
            # this must be a bipolar response since the charge
            # was not collected on the nearest pixel hypothesis
            is_bip = True
            print chg, 'IS BIPOLAR', dist

def get_pixels(entry):
    """Get pixel coordinates from ntuple"""
    pixels_y = [y for y in entry.pixel_y]
    pixels_z = [z for z in entry.pixel_z]
    pixels_y.sort()
    pixels_z.sort()

    global _pixel_spacing, _eps
    _pixel_spacing = np.round(abs(pixels_y[1]-pixels_y[0]), 2)
    _eps = np.round(np.sqrt(2*(0.5*_pixel_spacing)**2), 2)
    print _pixel_spacing, _eps
    return pixels_y, pixels_z

def sim(response, ntuple=None, **kwds):
    """This method will loop over the events in an ntuple.
    The ntuple, created with larg4, is expected to have
    the following information: 
        1) 3D coordinates of charge clouds close to 
           pixel plane.
        2) 2D coordinates of pixels
    """

    rfile = ROOT.TFile.Open(ntuple)
    rtree = rfile.Get('myana/anatree')

    rsp_sp, rsp_ep, rsp_fn = None, None, None
    for arr in response.data:
        if arr.name == 'start_points':
            rsp_sp = arr.data
        elif arr.name == 'end_points':
            rsp_ep = arr.data
        elif arr.name == 'responses':
            rsp_fn = arr.data

    # Using first entry to fill pixel coordinates
    pixels_y, pixels_z = None, None
    for entry in rtree:
        pixels_y, pixels_z = get_pixels(entry)
        break

    for entry in rtree:
        print 'drifting event', entry.event
        calculate(entry, pixels_y, pixels_z, rsp_sp, rsp_ep, rsp_fn, **kwds)
        break