#
# Script for running export procedure for muliple domains
#
import argparse
from subprocess import check_call
import os 

ids = [x for x in range(43,59)]

pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
for path in ids:
    com = 'export -s step -r '+str(path)+' -o path_resid_'+str(path)
    
    print '\nRunning for',str(path)
    dothis = pix+com
    print dothis
    check_call(dothis, shell=True)