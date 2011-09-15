#!/usr/bin/python

import sys,os,re

if len(sys.argv) != 2:
    print "usage: %s <dssp file>"%sys.argv[0]
    sys.exit(0)

lines = open(sys.argv[1]).readlines()
name = [re.match("\s*\d*\s*\d*\s.\s\S\s\s(.).*",l).group(1) for l in lines]

for n in name:
    print n,
print
