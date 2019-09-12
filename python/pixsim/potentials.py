

class weighting(object):
    def __init__(self, domain=0, potential=1.0, **kwds):
        self.domain = domain
        self.potential = potential
        print 'weighting for domain=%d potential=%f' % (domain, potential)

    def __call__(self, r, n, index, result):
        result[0] = 0.0
        if index == self.domain:
            result[0] = self.potential

class field_cage(object):
    '''
    Field cage potential. Pass domains and pot,
    applies voltage stepper on walls.
    '''

    def __init__(self, domains=None, drift_length=None, 
                 efield=None, anode_v=None, pad_v=None, 
                 grid_v=None, **kwds):
        self.domains      = domains
        self.drift_length = drift_length
        self.efield       = efield
        self.anode_v      = anode_v
        self.pad_v        = pad_v
        self.grid_v       = grid_v

    def __call__(self, r, n, domain_index, result):
        try:
            name = self.domains[domain_index]
        except:
            raise ValueError('Domain index not found:',domain_index)

        if name == 'walls':
            result[0] = self.anode_v - self.efield * r[0] 
        elif name == 'cathode':
            result[0] = self.anode_v - self.efield * self.drift_length
        elif name == 'anode':
            result[0] = self.anode_v
        elif 'pixel' in name:
            result[0] = self.pad_v
        elif name == 'grid':
            result[0] = self.grid_v
        else:
            raise ValueError('Unknown domain:', name)
