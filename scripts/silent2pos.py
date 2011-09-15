#!/usr/bin/python
from math import sqrt

import sys,os,re

if len(sys.argv) != 2:
    print "usage: %s <silent mode file>"%sys.argv[0]
    sys.exit(0)

#   1 L     0.000  132.598  178.785   24.914   54.621   53.432 S_0001_5441

lines = open(sys.argv[1]).readlines()
ss = [re.match("\s+(\d+)\s+.\s+(\S+)\s+(\S+)\s+(\S+)\s+.*\sS.(\d\d\d\d).(\d\d\d\d).*",l) for l in lines]
num = [int(m.group(1)) for m in ss if m is not None]
x = [float(m.group(2)) for m in ss if m is not None]
y = [float(m.group(3)) for m in ss if m is not None]
z = [float(m.group(4)) for m in ss if m is not None]
id1 = [m.group(5) for m in ss if m is not None]
id2 = [m.group(6) for m in ss if m is not None]

nres = max(num)
if 0:
    for ii in range(len(num)/nres):
        print id1[nres*ii]+id2[nres*ii],
        for jj in range(nres):
            print x[nres*ii+jj],
        print
        print id1[nres*ii]+id2[nres*ii],
        for jj in range(nres):
            print y[nres*ii+jj],
        print
        print id1[nres*ii]+id2[nres*ii],
        for jj in range(nres):
            print z[nres*ii+jj],
        print
