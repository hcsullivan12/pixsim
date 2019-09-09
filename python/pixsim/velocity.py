#!/usr/bin/env python
'''
Electron drift functionality.
'''

import math
import numpy
from collections import namedtuple

def mobility_function(Emag, temperature = 89):
    '''
    Return the mobility for the given magnitude of the 
    electric field and temperature. Units need to be kV, cm.
    '''

    Trel = temperature / 89
    a0=551.6           # cm2/sec
    a1=7953.7          # cm2/sec/kV
    a2=4440.43         # cm2/sec/kV^3/2
    a3=4.29            # cm2/sec/kV^5/2
    a4=43.63           # cm2/sec/kV^2
    a5=0.2053          # cm2/sec/kV^3

    e2 = Emag*Emag
    e3 = Emag*e2
    e5 = e2*e3
    e52 = math.sqrt(e5)
    e32 = math.sqrt(e3)

    Trel32 = math.sqrt(Trel*Trel*Trel)
    mu = (a0 + a1*Emag +a2*e32 + a3*e52)
    mu /= (1 + (a1/a0)*Emag + a4*e2 + a5*e3) * Trel32
    return mu

mobility = numpy.vectorize(mobility_function)

def drift(potential, linspaces, temperature=89, **kwds):
    '''
    Return an N-field matrix calculated assuming result holds a potential.
    '''
    print 'Drifting...'
    dxyz = [(ls[1]-ls[0]) for ls in linspaces]
    E = numpy.asarray(numpy.gradient(potential, *dxyz))
    # potential is in V and linspaces should be in cm, convert to kV/cm
    E /= 1000.
    Emag = numpy.sqrt(E[0]**2 + E[1]**2 + E[2]**2)
    # mu comes back as cm2/Vs
    mu = mobility(Emag, temperature)
    # convert velocity to cm/us
    vel = mu*E/1000.
    from pixsim.models import Array
    return [ Array(typename='gvector', name='vfield', data=vel) ]
