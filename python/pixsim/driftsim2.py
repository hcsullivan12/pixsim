import ROOT
import numpy as np
import sympy as sp

def form_waveforms(entry, backstep=1.0, **kwds):
    """
    Todo
    """
    chg_entries = len(entry.ides_y)
    for chg in range(0, chg_entries):
        pos0 = np.asarray([backstep, entry.ides_y[chg], entry.ides_z[chg]])


def get_pixels(entry):
    """Get pixel coordinates from ntuple"""
    pixels_y = [y for y in entry.pixel_y]
    pixels_z = [z for z in entry.pixel_z]
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
    #    form_waveforms(entry, **kwds)
    #    break



    