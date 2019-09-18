#
# Script for generating weighting fields for muliple domains
#
import argparse
from subprocess import check_call

sources = [n for n in range(10,26)]

pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
for src in sources:
    #par = '-p weight:domain='+str(dom[0])+' '
    com = 'raster -c weight_raster -r '+str(src)+' -n weight_raster_resid_'+str(src)

    print '\nRunning for',str(src)
    dothis = pix+com
    print dothis
    check_call(dothis, shell=True)