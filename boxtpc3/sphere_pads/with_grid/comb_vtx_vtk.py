import argparse

parser = argparse.ArgumentParser(description="Generate single vtk file from vtk file list.")
parser.add_argument("-s", "--source", type=str, help="File containing list of vtks.", required=True)
parser.add_argument("-t", "--type", type=str, help="Scalar value (potential/speed).", required=True)
args = parser.parse_args()

files, points, pot = list(), list(), list()
with open(args.source) as f:
    while True:
        line = f.readline().split()
        if len(line) < 1:
            break
        files.append(line[0])

for file in files:
    print file
    with open(file) as f:
        line = f.readline().split()
        while 'POINTS' not in line:
            line = f.readline().split()
        line = f.readline().split()
        while 'CELLS' not in line:
            count = 0
            while count < len(line) and len(line) >= 3:
                pos = [float(line[count]),float(line[count+1]),float(line[count+2])]
                points.append(pos)
                count += 3
            line = f.readline().split()
        line = f.readline().split()
        while 'LOOKUP_TABLE' not in line:
            line = f.readline().split()
        while True:
            line = f.readline().split()
            if len(line) < 1:
                break
            for ent in line:
                pot.append(float(ent))

# we have all the points and pot, now write the single file
with open('single_vtk.vtk', 'w') as f:
    f.write('# vtk DataFile Version 4.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(len(points))+' double\n')
    for pt in points:
        f.write(str(pt[0])+' '+str(pt[1])+' '+str(pt[2])+'\n')
    f.write('CELLS 1 '+str(len(points)+1)+'\n')
    out = str(len(points))
    for count, pt in enumerate(points):
        out += ' '+str(count)
    f.write(out+'\n')
    f.write('\n')
    f.write('CELL_TYPES 1\n')
    f.write('1\n')
    f.write('\n')
    f.write('POINT_DATA '+str(len(points))+'\n')
    f.write('SCALARS '+args.type+' double\n')
    f.write('LOOKUP_TABLE default\n')
    for pt in pot:
        f.write(str(pt)+'\n')




