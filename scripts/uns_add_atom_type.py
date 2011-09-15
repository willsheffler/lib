#!/usr/bin/python
import sys,os

if len(sys.argv) != 2:
    print "usage: uns_add_tome_type.py <uns_list>"
    sys.exit()

res1tores3 = {
    "A":"ALA",
    "C":"CYS",
    "D":"ASP",
    "E":"GLU",
    "F":"PHE",
    "G":"GLY",
    "H":"HIS",
    "I":"ILE",
    "K":"LYS",
    "L":"LEU",
    "M":"MET",
    "N":"ASN",
    "P":"PRO",
    "Q":"GLN",
    "R":"ARG",
    "S":"SER",
    "T":"THR",
    "V":"VAL",
    "W":"TRP",
    "Y":"TYR"
}

hbond_atom_types = {
('C','HG')  :22,
('D','OD1') :15,
('D','OD2') :15,
('E','OE1') :15,
('E','OE2') :15,
('H','HE2') :22,
('H','ND1') :8,
('K','1HZ') :22,
('K','2HZ') :22,
('K','3HZ') :22,
('N','1HD2'):22,
('N','2HD2'):22,
('N','OD1') :14,
('Q','1HE2'):22,
('Q','2HE2'):22,
('Q','OE1') :14,
('R','1HH1'):22,
('R','1HH2'):22,
('R','2HH1'):22,
('R','2HH2'):22,
('R','HE')  :22,
('S','HG')  :22,
('S','OG')  :13,
('T','HG1') :22,
('T','OG1') :13,
('W','HE1') :22,
('Y','HH')  :22,
('Y','OH')  :13  
}

print "tag pdb aa res atomnum score atomname sasa14 sasa10 sasa7 atomtype"
f = open(sys.argv[1])
for line in f.xreadlines():
    x = line.split()
    res = x[3]
    atom = x[7]
    if    atom is 'O': type = 20
    elif  atom is 'H': type = 25
    else: type = hbond_atom_types[(res,atom)]
    #print "%s\t%s\t%s"%(res,atom,type)
    print line[3:17],res1tores3[line[17]],line[18:-1],'  ',type
