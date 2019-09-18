#
# Script for running stepping procedure for muliple domains
#
import argparse
from subprocess import check_call
import os 

cwd = os.getcwd()
path = os.path.join(cwd,'vtxs')
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
for file in files:
    filepath = os.path.join(path,file)
    pidname = file.split('_')

    par = '-p step:stepfile='+filepath+' '
    com = 'step -n paths_for_'+pidname[0]
    
    print '\nRunning for',str(filepath)
    dothis = pix+par+com
    print dothis
    check_call(dothis, shell=True)