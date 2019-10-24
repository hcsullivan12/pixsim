import ROOT
import numpy as np
import sympy as sp
from bisect import bisect_left
from sets import Set

# default settings, will detect spacing later
_pixel_spacing = 0.4
_eps = np.round(np.sqrt(2*(0.5*_pixel_spacing)**2), 2)
_step_zs, _step_ys = Set(), dict()

def is_accurate(dist, eps=_eps):
    """Check against largest distance from pixel center"""
    #print _eps
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
                ret_fn, ret_ep = sat['rsp'], sat['ep']
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
    assert rot_rel_pos[1] >= 0 and rot_rel_pos[2] >= 0
    
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

def calculate(entry, pixels_y, pixels_z, rsp_data):
    """
    Calculate waveforms for each charge cloud. 
    Convolves charge deposits with field/electronics
    response function.

    We begin by guessing the nearest pixel to the drifted yz position.
    Extract the drift information and if the result is bipolar, use
    the end point to find the correct pixel, extract again.
    """
    chg_entries = len(entry.ides_y)
    for chg in range(0, chg_entries):
        # getting charge position
        pos0 = np.asarray([0.0, entry.ides_y[chg], entry.ides_z[chg]])
        
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
        print resp
        break

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
        sat_data = {'ep':ep, 'rsp':fn}

        if sz not in rsp_data.keys():
            rsp_data[sz] = {sy:sat_data}
        elif sy not in rsp_data[sz].keys():
            rsp_data[sz][sy] = sat_data
        else:
            raise ValueError('Data already contains yz pair!')
    return rsp_data

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
        break