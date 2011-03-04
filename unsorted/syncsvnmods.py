#!/usr/bin/env python
import sys,os

host = "dig"
rdir = os.getcwd().replace("/Users/sheffler/","")
if len(sys.argv) > 1: host = sys.argv[1]
if len(sys.argv) > 2: rdir = sys.argv[2]

print "syncing cdw to host '%s' dir '%s'"%(host,rdir)

files = {}
for l in os.popen("svn status").readlines():
	if l.endswith("~"): continue
	fn = l.split()[1]
	dn = os.path.dirname(fn)
	if not dn: dn = "."
	fn = os.path.basename(fn)
	if dn not in files: files[dn] = list()
	files[dn].append(fn)

for dn in files:
	fns = (" ").join([dn+"/"+x for x in files[dn]])
	cmd = "rsync -avz %s %s:%s/%s/"%(fns,host,rdir,dn)
	print cmd
	os.popen(cmd)
