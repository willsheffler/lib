#!/usr/bin/python

import sys,os

def line2pdbid(l):
	for x in l.split('_'):
		if len(x) == 4: return x

seenit = set()
for line in reversed(sys.stdin.readlines()):
	fn = line.split()[-1]
	prot = line2pdbid(fn)
	if prot in seenit: continue
	seenit.add(prot)
	print fn
