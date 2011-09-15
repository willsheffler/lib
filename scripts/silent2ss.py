#!/usr/bin/python

import sys,os,re

if len(sys.argv) != 2:
    print "usage: %s <silent mode file>"%sys.argv[0]
    sys.exit(0)

lines = open(sys.argv[1]).readlines()
ss = [re.match("\s+(\d+)\s+(.).*\sS.(\d\d\d\d).\d\d\d\d.*",l) for l in lines]
id = [m.group(3) for m in ss if m is not None]
num = [int(m.group(1)) for m in ss if m is not None]
ss = [m.group(2) for m in ss if m is not None]

nres = max(num)
for ii in range(len(num)/nres):
    print id[nres*ii],
    for jj in range(nres):
        if ss[nres*ii+jj] is 'L':
            print 0,
        if ss[nres*ii+jj] is 'E':
            print 1,
        if ss[nres*ii+jj] is 'H':
            print 2,
    print
