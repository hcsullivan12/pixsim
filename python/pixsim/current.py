import numpy as np
from pixsim.vector import Field
from pixsim.models import Array

def compute(efield, linspaces, paths, pnames):
    '''
    Apply RS theorem using weighting field and paths to compute current.
    Not worried about speed at the moment, just go brute force.
    '''

    # get the weighting field
    wght = Field(efield, linspaces)
    def weight(notused, r):
        return wght(r)
    
    arr = list()
    for name,path in zip(pnames,paths):
        print 'Computing waveform for path',name

        waveform = list()
        for step in path:
            # should be (x,y,z,u,v,w,t)
            xyz, uvw, time = step[0:3], step[3:6], step[-1] 
            w   = weight(time, xyz)
            print xyz, uvw, time,w,np.dot(uvw,w)
            waveform.append([time, np.dot(uvw,w)])
        arr.append( Array(typename='tuples', name='waveform_'+name, data = np.asarray(waveform)) )

    return arr