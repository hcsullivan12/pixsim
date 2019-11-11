import ROOT
import math
import matplotlib.pyplot as plt
from matplotlib  import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

pixels_z = []
pixels_y = []

currentZ = 0.2
while (currentZ+0.4) < 250:
    pixels_z.append(currentZ)
    pixels_z.append(-1*currentZ)
    currentZ += 0.4
currentY = 0.2
while (currentY+0.4) < 100:
    pixels_y.append(currentY)
    pixels_y.append(-1*currentY)
    currentY += 0.4

pixels_z.sort()
pixels_y.sort(reverse=True)
for c in range(0, len(pixels_z)):
    pixels_z[c] += 250

npixels = len(pixels_y)*len(pixels_z)
ncol    = len(pixels_z)
nrow    = len(pixels_y)

#for y in pixels_y:
#    print "%.2f" % y

def getPixelPosition(ch):
    row = int(math.floor(ch/ncol))
    col = ch % ncol
    return [pixels_z[col], pixels_y[row]]

f = ROOT.TFile.Open('nu_events.root')
tree = f.Get('anatree/anatree')

x, y, z, en = [], [], [], []

n = 0    
for entry in tree:
    n +=1
    if entry.event!=4: continue
    #continue
    '''
    print'pdg',    entry.nu_pdg[0]
    print'ndau',   entry.nu_ndau[0]
    print'ccnc',   entry.nu_ccnc[0]
    print'mode',   entry.nu_mode[0]
    print'intt',   entry.nu_inttype[0]
    print'startx', entry.nu_StartPointx[0]
    print'starty', entry.nu_StartPointy[0]
    print'startz', entry.nu_StartPointz[0]
    print'endx',   entry.nu_EndPointx[0]
    print'endy',   entry.nu_EndPointy[0]
    print'endz',   entry.nu_EndPointz[0]
    print'vtxx',   entry.nu_Vertexx[0]
    print'vtxy',   entry.nu_Vertexy[0]
    print'vtxz',   entry.nu_Vertexz[0]
    print'startpx', entry.nu_StartPx[0]
    print'startpy', entry.nu_StartPy[0]
    print'startpz', entry.nu_StartPz[0]
    print'endpx',    entry.nu_EndPx[0]
    print'endpy',    entry.nu_EndPy[0]
    print'endpz',    entry.nu_EndPz[0]
    print'ides', len(entry.ides_x)
    '''

    #for c in range(0, entry.geant_list_size):
    #    print entry.PDG[c], entry.StartPointx[c], entry.EndPointx[c]

    for c in range(0, len(entry.ides_x)):
        #print entry.ides_z[c]
        x.append(entry.ides_voxel_x[c])
        y.append(entry.ides_voxel_y[c])
        z.append(entry.ides_voxel_z[c])
        en.append(entry.ides_energy[c])
    #min_z = 10
    #for c in range(0, len(z)):
    #    if np.isnan(y[c]):print 'HEY'
    #    if min_z > abs(x[0]-x[c]):
    #        if abs(x[0]-x[c]) > 0.00001:
    #            min_z = abs(x[0]-x[c])
    #print "%.2f" % min_z
        #print "%.2f" % y[c]
    
    #for c in range(0, len(entry.ides_x)):
    #    print 'ides x', entry.ides_x[c]
    #    print 'ides y', entry.ides_y[c]
    #    print 'ides z', entry.ides_z[c]
    #    print 'ides E', entry.ides_energy[c]
    #    print 'ides tid', entry.ides_tid[c]
    #    print 'ides nel', entry.ides_numElectrons[c]
    #for c in range(0, len(entry.TrackId)):
    #    print 'tid', entry.TrackId[c]
    #    print 'pdg', entry.PDG[c]
    #    print 'starte', entry.StartEnergy[c]
    #    print 'ende', entry.EndEnergy[c]
    #    print 'startpx', entry.StartPx[c]
    #    print 'startpy', entry.StartPy[c]
    #    print 'startpz', entry.StartPz[c]
    #    print 'endpx', entry.EndPx[c]
    #    print 'endpy', entry.EndPy[c]
    #    print 'endpz', entry.EndPz[c]
    #    print 'startx', entry.StartPointx[c]
    #    print 'starty', entry.StartPointy[c]
    #    print 'startz', entry.StartPointz[c]
    #    print 'endx', entry.EndPointx[c]
    #    print 'endy', entry.EndPointy[c]
    #    print 'endz', entry.EndPointz[c]
    #    print 'mother', entry.Mother[c]
    #    #print 'nd', len(entry.Daughters[c])
    #    print 'isPrim', entry.isPrimary[c]
    #    #print 'proc', str(entry.Process[c])




_x = np.asarray(x)
_y = np.asarray(y)
_z = np.asarray(z)
_c = np.asarray(en)

#print _x

#for x in _x:
#	print "%.1f" % x
fig= plt.figure(0)
ax= fig.add_subplot(111, projection='3d')
##ax.patch.set_facecolor('darkblue')
#ax.set_axis_off()
ax.set_ylim([0, 360])
ax.set_xlim([0, 500])
ax.set_zlim([-100, 100])
img = ax.scatter(_z, _x, _y, s=0.1, c=_c, marker='o', cmap=cm.gnuplot, vmin=0.01, vmax=0.1)
ax.set_xlabel('z [cm]')
ax.set_ylabel('y [cm]')
ax.set_zlabel('x [cm]')
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
fig.colorbar(img)

fig= plt.figure(1)
ax= fig.add_subplot(111)
#ax.patch.set_facecolor('darkblue')
#ax.set_axis_off()
ax.set_ylim([0, 360])
ax.set_xlim([0, 500])
#ax.set_zlim([100, 200])
img = ax.scatter(_z, _x, c=_c, s=0.5, marker=',', cmap=cm.gnuplot, vmin=0.01, vmax=0.1)
fig.colorbar(img)

fig= plt.figure(2)
ax= fig.add_subplot(111)
#ax.patch.set_facecolor('darkblue')
#ax.set_axis_off()
ax.set_ylim([-100, 100])
ax.set_xlim([0, 500])
#ax.set_zlim([100, 200])
img = ax.scatter(_z, _y, c=_c, s=0.5, marker=',', cmap=cm.gnuplot, vmin=0.01, vmax=0.1)
fig.colorbar(img)

plt.show()
