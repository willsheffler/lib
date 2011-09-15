#!/usr/bin/python

import sys,os

if len(sys.argv) != 3:
    print "usage:pdb_starting_res_index.py <pdb dir> <pdb list>"
    sys.exit()

res3 = ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")


pdbdir  = sys.argv[1]
pdblist = sys.argv[2]

for pdbfilename in open(pdblist).readlines():
    #print "'%s'"%pdbfilename
    if len(pdbfilename) < 4:
        continue
    chain_inc = 1 # for now assume chain_inc
    pdbfile = open(pdbdir+'/'+pdbfilename[:-6]+".pdb")
    chain = pdbfilename[-6]
    resname = []
    resnum = []
    lastres = ""
    foundchain = 0
    for line in pdbfile:
        if line[0:4] == "ATOM":
            # is 21st always chain?
            # will there always be a CA?
            if chain == '_' or chain == line[21]:
                foundchain = 1
                #print "'%s'"%line[22:27],
                x = line[22:27]
                aa = line[17:20]
                if  x != lastres and aa in res3:
                    resname.append(aa)
                    resnum.append(x)
                    lastres = x
        elif line[:6] == "ENDMDL":
            break
        elif  foundchain and line[:4] == "TER ":
            break
#    print len(resnum)

#    print pdbfilename[:-5],len(resnum)                               
#    for i in range(len(resnum)-1):
#        if resnum[i] == resnum[i+1]:
#            print resnum[i],
        
    for i in range(len(resnum)):
        print "%s\t%i\t%s\t%s"%(pdbfilename[:-5],i+1,resnum[i],resname[i])



    pdbfile.close()
        

