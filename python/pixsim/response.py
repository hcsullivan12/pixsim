"""
Calculate total response from field and electronics response.
"""
import sympy as sp
import numpy as np
from pixsim.current import do_norm
from pixsim.models import Array

def response(steps, waveforms, gain=14.0, shape_time=1.0, **kwds):
    if len(waveforms) == 0:
        return list()

    A0 = gain #mV/fC
    tP = shape_time #us

    CT = 1/1.996
    den = tP*CT
    CA = 2.7433/den**4
    p0 = 1.477/den
    pi1, pi2 = 0.598/den, 1.299/den
    pr1, pr2 = 1.417/den, 1.204/den

    # find step width
    wvf = waveforms[0].data
    step = wvf[1][0]-wvf[0][0]

    s, t = sp.symbols('s, t')
    expr = (A0*CA)/((p0+s)*(pi1**2+(pr1+s)**2)*(pi2**2+(pr2+s)**2))
    sol = sp.inverse_laplace_transform(expr, s, t)

    ts = [n for n in np.arange(step, 5*tP, step)]
    ys = [sp.functions.re(sol.subs(t, tick)) for tick in ts]

    start, end, arrays = list(), list(), list()
    steps = [s for s in steps if s.name != 'vtxs']
    assert len(steps) == len(waveforms)
    
    for step, wvf in zip(steps, waveforms):
        if step.name == 'vtxs':
            continue
        print 'Calculating response for', step.name
        wvf = wvf.data[:,1]
        
        # we don't want to normalize bipolar pulses
        # normalize to 1 fC so current is in nA
        #if do_norm(wvf):
        #    wvf *= chg / abs(np.sum(wvf)*step)
        wvf = list(wvf)

        # convolve
        conv = sp.convolution(wvf, ys, dps=2)
        conv = np.asarray([sp.functions.re(x) for x in conv])

        start.append(step.data[0,0:3])
        end.append(step.data[-1,0:3])
        arrays.append(conv)

    return [Array(name='start_points', typename='tuples', data=np.asarray(start)),
            Array(name='end_points', typename='tuples', data=np.asarray(end)),
            Array(name='responses', typename='tuples', data=np.asarray(arrays))]