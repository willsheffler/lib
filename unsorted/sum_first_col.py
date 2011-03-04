#!/usr/bin/env python

import sys,os

s = 0
for line in sys.stdin.readlines():
   s += float(line.split()[0])
   print line,
print "total",s
