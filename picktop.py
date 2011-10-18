#!/usr/bin/python

import sys,os

def line2pdbid(l):
	for x in l.split('_'):
		if len(x) == 4: return x

m = {}
for line in sys.stdin.readlines():
	fn = line.split()[-1]
	prot = line2pdbid(fn[0:4])
	print prot
	if prot not in m: m[prot] = []
	m[prot].append( fn )

key = m.keys()
key.sort()
for k in key:
	for i in m[k]:
		print i

