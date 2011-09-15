#!/usr/bin/python
from math import sqrt
from sys import stderr

import sys,os,re

if len(sys.argv) != 2:
    print "usage: %OAs <silent mode file>"%sys.argv[0]
    sys.exit(0)

lines = open(sys.argv[1]).readlines()
ii = 0
for l in lines:
    tmp = re.match("(\s+\d+\s+.\s+\S+\s+\S+\s+\S+\s+.*\sS.)\d\d\d\d(.\d\d\d\d.*)",l)
    if tmp is not None:
        print tmp.group(1)+"%0.4i"%(ii)+tmp.group(2)
    else:
        tmp = re.match("(.*S.)\d\d\d\d(.\d\d\d\d.*)",l)
        if tmp is not None:
            ii += 1
            print tmp.group(1)+"%0.4i"%(ii)+tmp.group(2)
        else:
            print l[:-1]
