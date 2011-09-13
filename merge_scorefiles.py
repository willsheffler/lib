import sys

a = {}
with open(sys.argv[1]) as f:
	for line in f.readlines():
		pdb = line.split()[-1]
		a[pdb] = line.strip()
		
b = {}
with open(sys.argv[2]) as f:
	for line in f.readlines():
		pdb = line.split()[-1]
		b[pdb] = line.strip()
		
for k in a.keys():
	if k in b:
		print a[k].replace(k,"")+ b[k]