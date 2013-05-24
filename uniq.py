#!/usr/bin/env python
import sys
from optparse import OptionParser

start = int(sys.argv[1] if len(sys.argv) > 1 else 0)
end   = int(sys.argv[2] if len(sys.argv) > 2 else 0)

seenit = set()
for l in sys.stdin.xreadlines():
	if not end:
	      h = l.split()[start]
	else: h = l[start:end]
	if h in seenit: continue
	seenit.add(h)
	print l,
