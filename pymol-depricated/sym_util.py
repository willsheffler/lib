import sys,string,re,gzip,itertools,os
sys.path.append(os.path.expanduser("~/lib"))
sys.path.append(os.path.expanduser("~/lib/pymol"))
from LA import *
from pymol_util import *
from pymol import cmd

def c2axis(sele,alignsele=None,chains=["A","B"]):
	if alignsele is None: alignsele = sele
        #cmd.create('tmp98367598',sele)        
        #sele = 'tmp98367598'
	cmd.remove(sele+" and resn HOH")
	trans(sele,-com(alignsele))
	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
	print len(a),len(b)
	if len(a) != len(b) or len(a) == 0:
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord) ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
        #cmd.delete('tmp98367598')
	return axis

def alignc2(sele,alignsele=None,tgtaxis=Vec(0,0,1),chains=["A","B"]):
	if alignsele is None: alignsele = sele
	axis = c2axis(sele,alignsele,chains)
	# print "axis of rotation:",axis
	alignaxis(sele,tgtaxis,axis,Vec(0,0,0))
	return True

def c3axis(sele,alignsele=None,chains=["A","B","C"]):
	if alignsele is None: alignsele = sele
        #cmd.create('tmp98367598',sele)        
        #sele = 'tmp98367598'
	cmd.remove(sele+" and resn HOH")
	cen = com(alignsele)
	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
	# print "subunit lengths:",len(a),len(b),len(c)
	if len(a) != len(b) or len(a) != len(c) or len(a) == 0:
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord)+Vec(c[i].coord) - 3*cen ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
        #cmd.delete('tmp98367598')
	return axis

def alignc3(sele,alignsele=None,tgtaxis=Vec(0,0,1),chains=["A","B","C"]):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	axis = c3axis(sele,alignsele,chains)
	# print "axis of rotation:",axis
	alignaxis(sele,tgtaxis,axis,Vec(0,0,0))
	return True

def c4axis(sele,alignsele=None,chains=["A","B","C","D"]):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	trans(sele,-com(alignsele))
	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
	d = cmd.get_model(alignsele+" and chain "+chains[3]+" and name CA").atom
	# print "subunit lengths:",len(a),len(b),len(c)
	if len(a) != len(b) or len(a) != len(c) or len(a) == 0 or len(d) != len(a):
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord)+Vec(c[i].coord)+Vec(d[i].coord) ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
	return axis

def alignc4(sele,alignsele=None,tgtaxis=Vec(0,0,1),chains=["A","B","C","D"]):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	axis = c3axis(sele,alignsele,chains)
	# print "axis of rotation:",axis
	alignaxis(sele,tgtaxis,axis,Vec(0,0,0))
	return True

def c5axis(sele,alignsele=None,chains=["A","B","C","D","E"]):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	trans(sele,-com(alignsele))
	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
	d = cmd.get_model(alignsele+" and chain "+chains[3]+" and name CA").atom
	e = cmd.get_model(alignsele+" and chain "+chains[4]+" and name CA").atom
	# print "subunit lengths:",len(a),len(b),len(c)
	if len(a) != len(b) or len(a) != len(c) or len(a) == 0:
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord)+Vec(c[i].coord)+Vec(d[i].coord)+Vec(e[i].coord) ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
	return axis

def alignc5(sele,alignsele=None,tgtaxis=Vec(0,0,1),chains=["A","B","C","D","E"]):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	axis = c5axis(sele=sele,alignsele=alignsele,chains=chains)
	print "axis of rotation:",axis,"to",tgtaxis
	alignaxis(sele,tgtaxis,axis,Vec(0,0,0))
	return True

def myint(s):
   i = len(s)
   while i > 0 and not s[:i].isdigit(): i -= 1
   if not i: return None
   return int(s[:i])

def homogenizechains(sel1,sel2):
   cmd.remove("hydro")
   cmd.remove("resn HOH")
   cmd.remove("(HET and not resn MSE+CSW)")
   a = cmd.get_model("%s and name ca"%(sel1))
   b = cmd.get_model("%s and name ca"%(sel2))
   sa = "".join([name1[x.resn] for x in a.atom])
   sb = "".join([name1[x.resn] for x in b.atom])
   if sa==sb: return True
   ra = [myint(x.resi) for x in a.atom]
   rb = [myint(x.resi) for x in b.atom]
#   if max(ra) - min(ra) + 1 != len(ra): print "missing residue numbers",max(ra),min(ra),len(ra)
#   if max(rb) - min(rb) + 1 != len(rb): print "missing residue numbers",rb
   mla,mua,mlb,mub = lcs(sa,sb)
   bla,bua,blb,bub = lcs(sa[  :mla],sb[  :mlb])
   ala,aua,alb,aub = lcs(sa[mua+1:],sb[mub+1:])
   ra = ra[mla:(mua+1)]
   rb = rb[mlb:(mub+1)]
   if len(ra[bla:(bua+1)]) > 10:
      ra = ra[bla:(bua+1)] + ra[mla:(mua+1)] + ra[ala:(aua+1)]
      rb = rb[blb:(bub+1)] + rb[mlb:(mub+1)] + rb[alb:(aub+1)]
   if len(ra[ala:(aua+1)]) > 10:
      ra += ra[ala:(aua+1)]
      rb += rb[alb:(aub+1)]      
   for c,i in getres("%s"%(sel1)):
      if not i in ra: cmd.remove("%s and resi %i"%(sel1,i))
   for c,i in getres("%s"%(sel2)):
      if not i in rb: cmd.remove("%s and resi %i"%(sel2,i))   
   return False

def pickandfixchains(N,sel="all"):
   # find chains 
   # homogenize all pairs until fixed
   cc = []
   for c in getchain(sel):
      cc.append((-cmd.select("%s and chain %s and name CA"%(sel,c)),c))
   cc.sort()
   chains = [x[1] for x in cc[:N]]
   done = False
   count = 0
   while not done:
      if count > 10: break
      count += 1
      done = True;
      random.shuffle(chains)
      for i in range(1,len(chains)):
         done = done and homogenizechains(sel,chains[0],chains[i])
   print chains
   if N is 2: alignc2(sel,"name ca",chains=chains)
   if N is 3: alignc3(sel,"name ca",chains=chains)
   if N is 4: alignc4(sel,"name ca",chains=chains)
   if N is 5: alignc5(sel,"name ca",chains=chains)
   chains.sort()
   return chains[0]
   

def processhomomers():
   o = open("log",'w')
   for n in (2,3,4,5):
      files = glob.glob("c%ipdb/*.pdb.gz"%n)
      random.shuffle(files)
      for f in files:
         o.write(f+"\n")
         o.flush()
         cmd.delete("all")
         try:
            cmd.load(f)
            c = pickandfixchains(n)
            cmd.save("c%ia/"%n+f[3:-3],"chain %s"%c)
         except:
            print "fail on",f
   o.close()


def mki213(N, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("i213_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = com(sel)
	c3 = Vec(0, 0, 0)
	cmd.create( 'tmp80345769083457', sel)
	a2 = c2axis('tmp80345769083457')
	cmd.delete( 'tmp80345769083457')
	a3 = Vec(0, 0, 1)
	cmd.create('base80345769083457', sel+" and chain A and visible")
	seenit = []
	R2 = [rotation_matrix(a2, 0), rotation_matrix(a2, 180), ]
	R3 = [rotation_matrix(a3, 0), rotation_matrix(a3, 120), rotation_matrix(a3, 240), ]
	C = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	print a2, c2, a3, c3
	for i21 in range(2):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(2 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(2 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(2 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(2 if N > 8 else 1):
										test = Vec(0, 0, 0)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "i213_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457 and name n+ca+c')
										rot(n, a2, i21*180.0, c2)
										rot(n, a3, i32*120.0, c3)
										rot(n, a2, i22*180.0, c2)
										rot(n, a3, i33*120.0, c3)
										rot(n, a2, i23*180.0, c2)
										rot(n, a3, i34*120.0, c3)
										rot(n, a2, i24*180.0, c2)
										rot(n, a3, i35*120.0, c3)
										rot(n, a2, i25*180.0, c2)
	print len(seenit)
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def viewi213(sel = "all"):
	cmd.hide('ev')
	cmd.show('rib')
	mki213(sel)
	cmd.show('car', 'not i213*')
	cmd.hide('rib', 'not i213*')
	cmd.show('lines', '(byres (%s and not i213* and chain A) within 7.0 of (%s and not i213* and chain B))'%(sel, sel))
	cmd.show('lines', '(byres (%s and not i213* and chain B) within 7.0 of (%s and not i213* and chain A))'%(sel, sel))









def mkp23(N, R=43.5, i=0, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("p23_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = Vec(0, 0, 0)
	c3 = Vec(R,R,-R)
	cmd.create( 'tmp80345769083457', sel)
	cmd.delete( 'tmp80345769083457')
	a3 = [Vec(0,0,0),Vec(1,1,1),Vec(-1,-1,-1)]
	a2 = [Vec(0,0,0),Vec(1,0,0),Vec(0,1,0),Vec(0,0,1)]
	cmd.create('base80345769083457', sel+" and visible")
	seenit = []
	R2 = [rotation_matrix(a2[1],  0), # hack
	      rotation_matrix(a2[1],180),
	      rotation_matrix(a2[2],180),
	      rotation_matrix(a2[3],180) ]
	R3 = [rotation_matrix(a3[1],  0), # hack!
	      rotation_matrix(a3[1],120),
	      rotation_matrix(a3[2],120), ]
	C = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	print a2, c2, a3, c3
	for i21 in range(4):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(4 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(4 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(4 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(4 if N > 8 else 1):
										test = Vec(1, 1, 1)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "p23_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457')
										if i21 > 0: rot(n, a2[i21], 180.0, c2)
										if i32 > 0: rot(n, a3[i32], 120.0, c3)
										if i22 > 0: rot(n, a2[i22], 180.0, c2)
										if i33 > 0: rot(n, a3[i33], 120.0, c3)
										if i23 > 0: rot(n, a2[i23], 180.0, c2)
										if i34 > 0: rot(n, a3[i34], 120.0, c3)
										if i24 > 0: rot(n, a2[i24], 180.0, c2)
										if i35 > 0: rot(n, a3[i35], 120.0, c3)
										if i25 > 0: rot(n, a2[i25], 180.0, c2)
	print "seen:",len(seenit)
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def selbycomp(trn=0):
	cmd.select("TRI1","TRI and chain A+B+C")
	cmd.select("TRI2","TRI and chain D+E+F")
	cmd.select("TRI3","TRI and chain G+H+I")
	cmd.select("TRI4","TRI and chain J+K+L")
	cmd.select("TRI5","TRI and chain M+N+O")
	cmd.select("TRI6","TRI and chain P+Q+R")
	cmd.select("TRI7","TRI and chain S+T+U")
	cmd.select("TRI8","TRI and chain V+W+X")
	cmd.select("DIM1","DIM and chain A+D")
	cmd.select("DIM2","DIM and chain B+G")
	cmd.select("DIM3","DIM and chain C+J")
	cmd.select("DIM4","DIM and chain E+U")
	cmd.select("DIM5","DIM and chain F+R")
	cmd.select("DIM6","DIM and chain H+T")
	cmd.select("DIM7","DIM and chain I+O")
	cmd.select("DIM8","DIM and chain K+Q")
	cmd.select("DIM9","DIM and chain L+N")
	cmd.select("DIM10","DIM and chain M+V")
	cmd.select("DIM11","DIM and chain P+W")
	cmd.select("DIM12","DIM and chain X+S")
	cmd.delete("LINE*")
	cmd.delete("serf*")

	cmd.do("""alter all, b=50
	alter all, q=1
	set gaussian_resolution,8""")
	ISO="""map_new map%s, gaussian, 2, %s, 10
	isosurface surf%s, map%s"""


	for i in range(1, 9): 
		cmd.do(ISO%(("TRI%i"%i,)*4))
		cmd.color(COLORS[i-1],"surfTRI%i"%i)
		c = com("TRI%i"%i)
		# trans("TRI%i"%i,trn*c.normalized())
		obj = [
			CYLINDER, 
		   	0.0, 0.0, 0.0,
		   	1.6*c.x, 1.6*c.y, 1.6*c.z,
			1.5,
			0.1,0.1,0.1,0.1,0.1,0.1,
		]                                                                                            
		cmd.load_cgo(obj,'LINETRI%i'%i)
	for i in range(1,13): 
		cmd.do(ISO%(("DIM%i"%i,)*4))
		cmd.color(COLORS[i+7],"surfDIM%i"%i)
		c = com("DIM%i"%i)
		# trans("DIM%i"%i,trn*com("DIM%i"%i).normalized())
		obj = [
			CYLINDER, 
		   	0.0, 0.0, 0.0,
		   	1.3*c.x, 1.3*c.y, 1.3*c.z,
			1.0,
			0,0,1,0,0,1
		]                                                                                            
		cmd.load_cgo(obj,'LINEDIM%i'%i)






















