import ROOT
import numpy as np
import sympy as sp
from bisect import bisect_left

def get_nearest(pixel_pos, pos):
    """
    Assumes pixel_pos is sorted. Returns closest value to pos.
    """
    loc = bisect_left(pixel_pos, pos)
    if loc == len(pixel_pos):
        loc = loc-1
    elif loc != 0:
        before = pixel_pos[loc - 1]
        after = pixel_pos[loc]
        if after - pos > pos - before:
           loc = loc-1
    return loc

def nearest_pixel(pos, pixel_y, pixel_z, eps=0.2):
    """
    Find nearest pixel ID. Pixels are numbered
    in increasing order from top-left to right
    and down.

    With this scheme, since the arrays are assumed 
    to be sorted, convert y accordingly. 

    Returns a number in range [0,1,2,...,n-1]
    """
    row = get_nearest(pixel_y, pos[1])
    col = get_nearest(pixel_z, pos[2])
    near_pos = [pos[0], pixel_y[row], pixel_z[col]]
    row = len(pixel_y) - get_nearest(pixel_y, pos[1]) - 1

    dist = [(x-y)**2 for x,y in zip(pos, near_pos)]
    dist = sum(dist)**0.5
    if dist > eps:
        raise ValueError('Nearest pixel is further than epsilon=', eps, 'away.')
    return row * len(pixel_z) + col

def pixel_position(pid, pixel_y, pixel_z):
    """
    Get pixel position from pixel ID.
    Remember, pixel_y is in reverse id 
    order.
    """
    upper_lim = len(pixel_y)*len(pixel_z)-1
    if pid > upper_lim or pid < 0:
        raise ValueError('Pixel ID=', pid, 'out of range=[0,', upper_lim, ']')
    row = pid/len(pixel_z)
    col = pid - row * len(pixel_z)
    return [0.0, pixel_y[len(pixel_y)-1-row], pixel_z[col]]

def calculate(entry, pixels_y, pixels_z, backstep=1.0, **kwds):
    """
    Todo
    """
    chg_entries = len(entry.ides_y)
    for chg in range(0, chg_entries):
        # getting charge location and pixel info
        pos0 = np.asarray([backstep, entry.ides_y[chg], entry.ides_z[chg]])
        near_pix_id = nearest_pixel(pos0, pixels_y, pixels_z)
        pix_pos = pixel_position(near_pix_id, pixel_y, pixel_z)

        # convert relative pos to 1st quad and reflect over y=z if needed
        rel_pos = [abs(x-y) for x, y in zip(pos0, pix_pos)]
        if rel_pos[1] > rel_pos[2]:
            temp = rel_pos[1]
            rel_pos[1] = rel_pos[2]
            rel_pos[2] = temp
        
        # find nearest response function endpoints
        resp, sp, ep = get_response(rel_pos)

def get_pixels(entry):
    """Get pixel coordinates from ntuple"""
    pixels_y = [y for y in entry.pixel_y]
    pixels_z = [z for z in entry.pixel_z]
    pixels_y.sort()
    pixels_z.sort()
    return pixels_y, pixels_z

def sim(ntuple=None, **kwds):
    """This method will loop over the events in an ntuple.
    The ntuple, created with larg4, is expected to have
    the following information: 
        1) 3D coordinates of charge clouds close to 
           pixel plane.
        2) 2D coordinates of pixels
    """

    rfile = ROOT.TFile.Open(ntuple)
    rtree = rfile.Get('myana/anatree')

    # Using first entry to fill pixel coordinates
    pixels_y, pixels_z = None, None
    for entry in rtree:
        pixels_y, pixels_z = get_pixels(entry)
        break

    #for entry in rtree:
    #    print 'drifting event', entry.event
    #    calculate(entry, pixels_y, pixels_z, **kwds)
    #    break



    