#!/usr/bin/python
from math import sqrt
from sys import stderr

import sys,os,re

if len(sys.argv) != 2:
    print "usage: %s <silent mode file>"%sys.argv[0]
    sys.exit(0)

lines = open(sys.argv[1]).readlines()
ss = [re.match("\s+(\d+)\s+.\s+(\S+)\s+(\S+)\s+(\S+)\s+.*\sS.(\d\d\d\d).(\d\d\d\d).*",l) for l in lines]
num = [int(m.group(1)) for m in ss if m is not None]
x = [float(m.group(2)) for m in ss if m is not None]
y = [float(m.group(3)) for m in ss if m is not None]
z = [float(m.group(4)) for m in ss if m is not None]
id1 = [m.group(5) for m in ss if m is not None]
id2 = [m.group(6) for m in ss if m is not None]

nres = max(num)
stderr.write(`nres`+' residues\n')
for ii in (0,1):#range(len(num)/nres):
    stderr.write(`ii`+'\n')
    for jj in range(nres):
        stderr.write(`x[nres*ii+jj]`+' '+
                     `y[nres*ii+jj]`+' '+
                     `z[nres*ii+jj]`+'\n')
        for kk in range(jj):
            print (x[nres*ii+jj]-x[nres*ii+kk])**2 + \
                  (y[nres*ii+jj]-y[nres*ii+kk])**2 + \
                  (z[nres*ii+jj]-z[nres*ii+kk])**2,
        print '9 '+'0 '*(nres-jj-1)
        
        
