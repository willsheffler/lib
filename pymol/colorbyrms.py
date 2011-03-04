import pymol
from math import *

def resrms(sel1,sel2):
    r1 = cmd.get_model(sel1).get_residues()
    r2 = cmd.get_model(sel2).get_residues()
    assert len(r1)==len(r2)
    rrms = []
    for ii in range(1,len(r1)+1):
        s1 = sel1+' and resi '+`ii`
        s2 = sel2+' and resi '+`ii`
        rms = cmd.rms_cur( s1 , s2 )
        rrms.append( rms )
        print rms, type(rms)
        cmd.do("alter %s , b=%f"%(s1,(max(0.0,min(rms,10)))))
    return resrms


def ballrelrms(sel1,sel2):
    r1 = cmd.get_model(sel1).get_residues()
    r2 = cmd.get_model(sel2).get_residues()
    assert len(r1)==len(r2)
    brms = []
    for ii in range(1,len(r1)+1):
        s1 = sel1+' within 10 of (name ca and resi '+`ii`+')'
        s2 = sel2
        rms = cmd.rms( s1 , s2 )
        brms.append( rms )
        print rms, type(rms)
        cmd.do("alter %s , b=%f"%(sel1 + " and resi "+`ii`,sqrt(max(0.5,min(rms,10)))))
    return brms

cmd.extend('resrms',resrms)
cmd.extend('ballrms',ballrelrms)
