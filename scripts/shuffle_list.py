#!/usr/bin/python
import random,sys

f = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
try:
	list = [x.strip() for x in f.readlines()]
	random.shuffle(list)
	for l in list: print l
except IOError:
	pass
