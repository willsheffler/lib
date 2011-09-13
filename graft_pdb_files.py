import sys,os
from LA import *
from pyPDB import *

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

