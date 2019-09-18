#
# Script for running stepping procedure for muliple domains
#
import argparse
from subprocess import check_call
import os 

raster = [x for x in range(10,26)]
pid    = [8, 9, 10, 11, 14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29]
mapping = {x:y for x,y in zip(raster,pid)}

pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
for ras,pid in mapping.iteritems():
    raster_name = 'weight_raster_'+str(ras)
    step_name   = 'paths_for_pixel'+str(pid)

    com = 'current -s '+step_name+' -r '+raster_name+' -n waveform_rasid'+str(ras)+'_pid'+str(pid)
    
    print '\nRunning for',ras,pid
    dothis = pix+com
    print dothis
    check_call(dothis, shell=True)