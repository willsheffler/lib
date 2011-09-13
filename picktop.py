#!/usr/bin/python

import sys,os

m = {}
for line in sys.stdin.readlines():
	fn = line.split()[-1]
	prot = fn[0:4]
	if prot not in m: m[prot] = []
	m[prot].append( fn )

key = m.keys()
key.sort()
for k in key:
	for i in m[k]:
		print i

