

class field_cage(object):
    '''
    Field cage potential. Pass domains and pot,
    applies voltage stepper on walls.
    '''

    def __init__(self, domains, drift_length, efield=500., anode_v=0., **kwds):
        self.domains      = domains
        self.drift_length = drift_length
        self.efield       = efield
        self.anode_v      = anode_v

    def __call__(self, r, n, domain_index, result):
        try:
            name = self.domains[domain_index]
        except:
            print 'Domain index not found:',domain_index

        if name == 'walls':
            result[0] = self.anode_v - self.efield * r[0] 
        elif name == 'cathode':
            result[0] = self.anode_v - self.efield * self.drift_length
        elif name == 'anode':
            result[0] = self.anode_v
        else:
            raise ValueError('Unknown domain:', name)
