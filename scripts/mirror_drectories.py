import sys,os,shutil

list = [x.strip() for x in open(sys.argv[1]).readlines()]
oldloc = sys.argv[2]+'/'
newloc = sys.argv[3]+'/'

for ii in range(len(list)):
    fname = list[ii]
    if not os.path.exists(newloc+os.path.dirname(fname)):
	os.system("mkdir -p "+newloc+os.path.dirname(fname))
    for suffix in (".pdb",".dssp",".pdb.gz",".dssp.gz"):
	if os.path.exists(oldloc+fname+suffix) and not os.path.exists( newloc+fname+suffix ):
	    print float(ii)/float(len(list)),
	    print "copy", oldloc+fname+suffix, newloc+fname+suffix 
	    shutil.copyfile( oldloc+fname+suffix, newloc+fname+suffix )


