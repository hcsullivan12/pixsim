import numpy as np

class Box():
    def __init__(self, name, dx, dy, dz, center):
        self.name   = name
        self.hdim   = np.asarray([dx,dy,dz])
        self.center = center

    def inside(self, r):
        dist = abs(r - center)
        ins = True
        for d,hdim in zip(dist,self.hdim):
            if d > hdim:
                ins = False
        return ins

    def shape(self):
        return 'box'

    def info(self):
        return self.name, self.hdim, self.center, self.shape()

class Cylinder():
    def __init__(self, name, r, dz, center, dir):
        self.name   = name
        self.hdim   = np.asarray([r,dz])
        self.center = center
        self.dir    = dir/np.norm(dir)

    def inside(self, r):
        dist = r - self.center
        shadow = numpy.dot(r, self.dir)
        if abs(shadow) > self.height:
            return False
        shadow *= self.dir
        topt = r - shadow
        dist = np.sqrt(sum([t**2 for t in topt]))
        return dist <= self.radius
    
    def shape(self):
        return 'cylinder'
    
    def info(self):
        return self.name, self.hdim, self.center, self.shape()

class Sphere():
    def __init__(self, name, r, center):
        self.name   = name
        self.radius = r
        self.center = center

    def inside(self, r):
        dist = r - shadow
        dist = np.sqrt(sum([x**2 for x in dist]))
        return dist < self.radius

    def shape(self):
        return 'sphere'
    
    def info(self):
        return self.name, self.radius, self.center, self.shape()

