

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
            print 'Domain index not found:',domain_index

        if 'wall' in name:
            result[0] = self.anode_v - self.efield * r[0] 
        elif 'cathode' in name:
            result[0] = self.anode_v - self.efield * self.drift_length
        elif 'anode' in name:
            result[0] = self.anode_v
        elif 'pixel' in name:
            result[0] = self.pad_v
        elif 'grid' in name:
            result[0] = self.grid_v
        else:
            raise ValueError('Unknown domain:', name)
