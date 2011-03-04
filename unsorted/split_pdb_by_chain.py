#!/usr/bin/python
import sys,os,sets

for fn in sys.argv[1:]:
   print "splitting by chain:",fn
   lall = [l for l in open(fn).readlines() if l.startswith("ATOM  ")]
   call = sets.Set([l[21] for l in lall])
   for c in call:
      cfn = fn[:-4]+c+".pdb"
      o = open(cfn,'w')
      for l in lall:
         if l[21] == c:
            o.write(l)
      o.close()
