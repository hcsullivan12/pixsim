#
# Script for generating weighting fields for muliple domains
#
import argparse
from subprocess import check_call

parser = argparse.ArgumentParser(description="Generate weighting fields")
parser.add_argument("-s", "--source", type=str, help="The domain map.", required=True)
args = parser.parse_args()

domains = list()
with open(args.source) as f:
    while True:
        linevec = f.readline().split()
        if len(linevec) < 1:
            break
        dom,name = int(linevec[0]),str(linevec[1])
        cent = [float(linevec[2]),float(linevec[3]),float(linevec[4])]
        domains.append([dom,name,cent])

pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
for dom in domains:
    pos = dom[2]
    if abs(pos[1]) > 0.7 or abs(pos[2]) > 0.7:
        continue
    par = '-p weight:domain='+str(dom[0])+' '
    com = 'boundary -c weight -n weight_'+dom[1]

    print '\nRunning for',dom[1],'domain =',dom[0]
    dothis = pix+par+com
    print dothis
    check_call(dothis, shell=True)