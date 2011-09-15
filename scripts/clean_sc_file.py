import sys,os

files = sys.argv[1:]

for f in files:
   out = None
   N = 0
   for line in open(f).xreadlines():
      if N is 0:
         N = len(line.split())
         outname = f+".clean"
         if f.endswith(".sc"): outname = f[:-3]+".clean.sc"
         print "writing",outname
         out = open(outname,"w")
         out.write(line)
      else:
         if len(line.split()) == N:
            out.write(line)
   out.close()
