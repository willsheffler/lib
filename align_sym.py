import pickle,sys
from LA import *


def xform(x):
	v = Vec( 0.707107, 0.707107, 1.000000)/2.0
	x = rotation_matrix(v,180.0)*x
	return x/1.414214
	
with open(sys.argv[1]) as o: 
	for line in o:
		if not line.startswith("xyz"): print line,
		else:
			s = line.split()			
			x = Vec( *s[2].split(",") )
			y = Vec( *s[3].split(",") )
			x,y = xform(x),xform(y)
			print s[0],s[1],str(x).replace(", ",","),str(y).replace(", ",","),s[4]

