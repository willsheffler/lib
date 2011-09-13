#!/usr/bin/python

import sys,os

seenit = set()
for line in reversed(sys.stdin.readlines()):
	fn = line.split()[-1]
	prot = fn[0:4]
	if prot in seenit: continue
	seenit.add(prot)
	print fn
