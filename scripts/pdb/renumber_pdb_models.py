#!/usr/bin/env python

import sys,os,gzip

for fname in sys.argv[1:]:
	with gzip.open(fname) as fin:
		with open(fname.split(".")[0]+"_nomdl.pdb",'w') as out:
			for line in fin.xreadlines():
				if not line.startswith("ATOM  ") or line.startswith("HETATM"):
					continue
				out.write(line)