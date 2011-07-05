import sys,os

r = set()
seen = set()
f = open(sys.argv[1])
for l in f.readlines():
    s = l.strip().split()
    if s[0]=="REMARK" and s[1]=="666" and s[2]=="MATCH" and s[3]=="TEMPLATE":
        r.add(int(s[6 ]))
        r.add(int(s[11]))
    elif (s[0]=="ATOM" or s[0]=="HETATM") and l[21]!=" ":
        i = int(l[22:27])
        if i in r or l[21]=="X":
            print l.strip()
        seen.add(i)
f.close()


