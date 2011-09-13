
import sys,os
from LA import *

class PDBLine(object):
    def __init__(self,l):
        self.line = l.strip()
        if not l.startswith("ATOM ") and not l.startswith("HETATM "):
            self.mknull()
        else:
            "0123456789 11 14 17 20 23 26 29 32 35 38 41 44 47 50 53 56 59 62 65"
            "|||||||||| |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | "
            "ATOM      1  N   ARG A   1       9.545   1.916 -42.900  1.00  0.00"
            self.ano   = int(l[7:12])
            self.an    = l[12:16]
            self.rn    = l[17:20]
            self.chain = l[21]
            self.ri    = int(l[23:28])
            self.xyz   = Vec(float(l[30:38]),float(l[39:46]),float(l[47:54]))
    def mknull(self):
        self.ano=0;self.an="";self.rn="";self.chain="";self.ri=0;self.x=0;self.y=0;self.z=0
    def dist(self,othr):
        return self.xyz.distance(othr.xyz)
    def __str__(self):
        return str(self.an)+": "+self.line

class Rsd(object):
    def __init__(self,atoms=[]):
        self.map = {}
        self.atoms = []
        for atom in atoms: self.addAtom(atom)
    def addAtom(self,pline):
        self.atoms.append(pline)
        self.map[pline.an.strip()] = pline
    def atom(self,an):
        return self.map[an]
    def has(self,an):
        return an in self.map
    def bbmatch(self,othr):
        if not self.has( "N") or not othr.has( "N"): return False
        if not self.has("CA") or not othr.has("CA"): return False
        if not self.has( "C") or not othr.has( "C"): return False
        if self.atom( "N").dist(othr.atom( "N")) > 0.01: return False
        if self.atom("CA").dist(othr.atom("CA")) > 0.01: return False
        if self.atom( "C").dist(othr.atom( "C")) > 0.01: return False
        return True
    def mutate(self,newr):
        self.map = newr.map
        self.atomr = newr.atoms
    def dump(self,o):
        for a in self.atoms:
            o.write(a.line+"\n")
    def dumpala(self,o):
        for a in self.atoms:
            if a.an.strip() in ("N","CA","C","O","CB"):
                o.write(a.line.replace(a.rn,"ALA")+"\n")
    def __str__(self):
        s = ""
        for a in self.atoms:
            s += str(a)+"\n"
        return s

class PDB(object):
    def __init__(self,plines):
        self.res = {}
        ri = None
        for pl in plines:
            if pl.ri != ri:
                ri = pl.ri
                self.res[pl.ri] = Rsd([pl])
            else:
                self.res[pl.ri].addAtom(pl)
    def __str__(self):
        s = ""
        for i,r in enumerate(self.res):
            s += str(i)+"\n"
            s += str(r)+"\n"
        return s
    def find_bbmatch(self,rsd):
        for r in self.res:
            if r.bbmatch(rsd): return r
        return None

class EDremark(object):
    def __init__(self,l):
        "0123456789 11 14 17 20 23 26 29 32 35 38 41 44 47 50 53 56 59 62 65"
        "|||||||||| |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | "
        "REMARK 666 MATCH TEMPLATE X CTP    0 MATCH MOTIF A ARG  204  1  1"
        self.cn1 = l[26]
        self.rn1 = l[28:31]
        self.ri1 = int(l[33:36])
        self.cn2 = l[49]
        self.rn2 = l[51:54]
        self.ri2 = int(l[56:60])
        def __str__(self):
            return str(self.cn1)+" "+str(self.rn1)+" "+str(self.ri1)+" "+str(self.cn2)+" "+str(self.rn2)+" "+str(self.ri2)

def pdb_from_file(fn):
    f = open(fn)
    lines = f.readlines()
    plines = []
    for l in lines:
        if l.startswith("ENDMDL"): break
        if l.startswith("ATOM") or l.startswith("HETATM"):
            plines.append( PDBLine(l) )
    p = PDB(plines)
    p.ed = [EDremark(l) for l in lines if l.startswith("REMARK 666")]
    f.close()
    return p

# count = 0
# for fn in sys.argv[2:]:
#     print count
#     count += 1
#     template = pdb_from_file(sys.argv[1])
#     p = pdb_from_file(fn)
#     o = open(fn+"_w_scaff.pdb","w")
#     for r in template.res:
#         m = p.find_bbmatch(r)
#         if m: m.dump(o)
#         else: r.dumpala(o)
#     o.close()

for fn in sys.argv[1:]:
 try:
  p = pdb_from_file(fn)
  for minfo in p.ed:
   if minfo.rn1=="CTP":
       c2 = p.res[1].atom('C2').xyz
       c5 = p.res[1].atom('C5').xyz
       c7 = p.res[1].atom('C7').xyz
       c8 = p.res[1].atom('C8').xyz
       c9 = p.res[1].atom('C9').xyz
       o1 = p.res[1].atom('O1').xyz
       o2 = p.res[1].atom('O2').xyz
       o3 = p.res[1].atom('O3').xyz
       o4 = p.res[1].atom('O4').xyz
       o5 = p.res[1].atom('O5').xyz
       o6 = p.res[1].atom('O6').xyz
       if minfo.rn2=="ARG":
           n1 = p.res[minfo.ri2].atom('NH1').xyz
           n2 = p.res[minfo.ri2].atom('NH2').xyz
           ne = p.res[minfo.ri2].atom('NE' ).xyz
           cz = p.res[minfo.ri2].atom('CZ' ).xyz
           d1 = min(o3.d(n1),o3.d(n2),o3.d(ne),o5.d(n1),o5.d(n2),o5.d(ne))
           d2 = min(o4.d(n1),o4.d(n2),o4.d(ne),o6.d(n1),o6.d(n2),o6.d(ne))
           print 28,d1
           print 28,d2
           a1 = max([angle(c2,c7,cz),angle(c8,c9,cz)])
           a2 = max([angle(c7,cz,x) for x in (n1,n2,ne)]+[angle(c9,cz,x) for x in (n1,n2,ne)])
           print 180,a1
           print 180,a2
       if minfo.rn2=="LYS":
           nz = p.res[minfo.ri2].atom('NZ').xyz
           d1 = o2.d(nz)
           print 30,d1
       if minfo.rn2 in ("GLU","ASP"):
           if minfo.rn2=="GLU":
               dec  = p.res[minfo.ri2].atom('CB').xyz
               deo1 = p.res[minfo.ri2].atom('OE1').xyz
               deo2 = p.res[minfo.ri2].atom('OE2').xyz
           if minfo.rn2=="ASP":
               dec  = p.res[minfo.ri2].atom('CG').xyz
               deo1 = p.res[minfo.ri2].atom('OD1').xyz
               deo2 = p.res[minfo.ri2].atom('OD2').xyz
           d1 = min([o1.d(deo1),o1.d(deo2)])
           print 30,d1
           if abs(angle(dec,deo1,o1)-120.0) < abs(angle(dec,deo2,o1)-120.0):
               a1 = angle(dec,deo1,o1)
           else: a1 = angle(dec,deo2,o1)
           print 120,a1
           if abs(angle(deo1, o1,c5)-108.5) < abs(angle(deo2, o1,c5)-108.5):
               a2 = angle(deo1,o1,c5)
           else: a2 = angle(deo2,o1,c5)
           print 108,a2

   else:
       n1,n2,ne,o1,o2 = None,None,None,None,None
       if minfo.rn1=="ARG":
           n1 = p.res[minfo.ri1].atom('NH1').xyz
           n2 = p.res[minfo.ri1].atom('NH2').xyz
           ne = p.res[minfo.ri1].atom('NE' ).xyz
           cz = p.res[minfo.ri1].atom('CZ' ).xyz
       else: pass#print "??????????"
       if   minfo.rn2=="ASP":
           o1 = p.res[minfo.ri2].atom('OD1').xyz
           o2 = p.res[minfo.ri2].atom('OD2').xyz
           c1 = p.res[minfo.ri2].atom('CB').xyz
           c2 = p.res[minfo.ri2].atom('CG').xyz
       elif minfo.rn2=="GLU":
           o1 = p.res[minfo.ri2].atom('OE1').xyz
           o2 = p.res[minfo.ri2].atom('OE2').xyz
           c1 = p.res[minfo.ri2].atom('CG').xyz
           c2 = p.res[minfo.ri2].atom('CD').xyz
       else: pass#print "??????????"
       d1 = min([o1.d(n1),o1.d(n2),o1.d(ne)])
       d2 = min([o2.d(n1),o2.d(n2),o2.d(ne)])
       print 30,d1
       print 30,d2
       a1 = angle(c1,c2,cz)
       a2 = max([angle(c2,cz,x) for x in (n1,n2,ne)])
       print 180,a1
       print 180,a2
 except:
     pass#print 'fail!'

