#!/usr/bin/env python
import sys,os

to = sys.argv[1]

for line in os.popen("svn status").readlines():
	if line[0] in "MA?":
		cmd = "rsync -a %s %s"%(line[1:].strip(),to)
		print cmd
	else: print "unknown:",line,
