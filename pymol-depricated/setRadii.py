import pymol
from pymol import cmd
import sys,os,random
from math import floor


#PyMOL>print dir(cmd.get_model('r1').atom[0])
#['__cmp__', '__doc__', '__getattr__', '__module__', 'b', 'bohr', 'chain', 'coord', 'defaults', 'flags', 'formal_charge', 'get_implicit_valence', 'get_mass', 'get_number', 'get_signature', 'has', 'hetatm', 'id', 'in_same_residue', 'index', 'name', 'new_in_residue', 'numeric_type', 'partial_charge', 'q', 'resi', 'resi_number', 'resn', 'segi', 'ss', 'stereo', 'symbol', 'vdw']

rainbow = ["0xFF0000","0xFF4C00","0xFF9900","0xFFE500","0xCCFF00","0x80FF00","0x33FF00","0x00FF19","0x00FF66","0x00FFB2","0x00FFFF","0x00B3FF","0x0066FF","0x0019FF","0x3300FF","0x8000FF","0xCC00FF","0xFF00E6","0xFF0099","0xFF004D" ]

def useRosettaRadii():
	cmd.alter("element C", "vdw=2.00")
	cmd.alter("element N", "vdw=1.75")
	cmd.alter("element O", "vdw=1.55")
	cmd.alter("element H", "vdw=1.00")
	cmd.alter("element P", "vdw=1.90")
	cmd.set("sphere_scale", 1.0)
	
cmd.extend('useRosettaRadii', useRosettaRadii)

def useOccColors(sel="all"):
	colors = rainbow
	random.shuffle(colors)
	for ii in range(len(colors)):
		cmd.color( colors[ii] ,"%s and q=%i"%(sel,ii))

def useTempColors(sel="all"):
	for a in cmd.get_model(sel).atom:
		q = a.b
		c = intcolors[ int(floor(q))%len(intcolors) ]
		cmd.color( c ,"%s and resi %s and name %s"%(sel,a.resi,a.name))

def useOccRadii(sel="all"):
	for a in cmd.get_model(sel).atom:
		q = a.q
		if q >= 3:
			print "shrik radius"
			q <- 0.1
		cmd.alter("%s and resi %s and name %s"%(sel,a.resi,a.name),"vdw=%f"%(q))



def useTempRadii(sel="all"):
	for a in cmd.get_model(sel).atom:
		bfac = a.b
		if bfac >= 3:
			print "shrik radius"
			bfac <- 0.1
		cmd.alter("%s and resi %s and name %s"%(sel,a.resi,a.name),"vdw=%f"%(bfac))


def loadPackingPDB(file,name=None):
	if name is None:
		cmd.load(file)
		name = os.path.basename(file)[:-4]
	else:
		cmd.load(file,name)
	cmd.hide('all')
	cmd.do("useRosettaRadii")
	cmd.select("cavities","resi 500-999 and "+name)
	cmd.select("protein","resi 0-499 and "+name)
	cmd.do("useTempRadii cavities")
	cmd.show('spheres',"cavities")
	cmd.do("useOccColors cavities")
	cmd.show("lines","protein")
	cmd.show("cartoon","protein")
	cmd.color("white","protein")
	#cmd.delete("cavities")
	#cmd.delete("protein")

cmd.extend("loadPackingPDB",loadPackingPDB)		   
cmd.extend("useOccRadii",useOccRadii)
cmd.extend("useOccColors",useOccColors)
cmd.extend("useTempRadii",useTempRadii)
cmd.extend("useTempColors",useTempColors)



def which(v,x):
	r = []
	for ii in range(len(v)):
		if v[ii]==x:
			r.append(ii)
	return r

## cmd.reinitialize()
## cmd.do("")

## wd = "/users/sheffler/project/modelpicking/pdbs/"
## pdbs = os.listdir(wd)
## pdbs = [wd+p for p in pdbs if p.startswith('cf') and p.endswith('.pdb')]

## cmd.create('n','none')
## cmd.load(wd+'1hz6_0001.pdb','n')
## #cmd.hide('lines')
## #cmd.show("cartoon")

## for ii in range(50):#len(pdbs)):
## 	print ii
## 	p = pdbs[ii]
## 	cmd.create('d'+`ii`,'none')
## 	cmd.load(p,'d'+`ii`)
## 	cmd.hide('all')
## 	#cmd.show("cartoon")

## cmd.do("dss")
## cmd.show('lines','name ca or name n or name c')
