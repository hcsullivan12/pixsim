

start=7
end=42
inc=1
with open('temp.txt', 'w') as f:
	pid = 1
	surf=start
	while surf<=end:
		s = 'Physical Surface(\"pixel'+str(pid)+'\") = {'+str(surf)+'};'
		f.write(s+'\n')
		f.write('//+\n')
		surf += inc
		pid += 1
