#!/usr/bin/env python

import sys, os, shutil

if(len(sys.argv) > 1): target = sys.argv[1].rstrip('/')
else:	               target = os.path.expanduser('~/Dropbox/dev/misc')

for line in sys.stdin.readlines():
	if not line.strip(): continue
	fn = line.split()[-1].strip()
	if (fn.startswith("cmake/") or
		fn.endswith(".dylib") or
		fn.endswith(".os") or
		fn.endswith(".so") or
		fn.endswith(".o") or
		fn.endswith(".pyc") or
		fn.endswith("~") or
		os.path.split(fn)[1].startswith("._") or
		fn.count("/.svn/")
	):
		continue
	if os.path.isdir(fn): fn += '/'
	newfn = target+"/"+fn
	diffs = 1
	if not os.path.exists(os.path.dirname(newfn)):
		os.makedirs(os.path.dirname(newfn))
	else:
		diffs = int(os.popen("diff -rq %s %s|wc -l"%(fn,newfn)).read()) if os.path.exists(newfn) else -1
	if diffs:
		# print "copying %s to %s %i" % (fn,newfn,diffs)
		cmd = "rsync --exclude .svn -a %s %s"%(fn,newfn)
		print cmd
		os.system(cmd)
