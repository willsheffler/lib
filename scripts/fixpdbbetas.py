import sys,os

for fname in sys.argv[1:]:
    print fname
    pdb = open(fname)
    lines = pdb.readlines()
    pdb.close()
    pdb = open(fname,'w')
    for line in lines:
        if line[:5] == "ATOM ":
            pdb.write(line[:60]+ "  0.00\n")
        else:
            pdb.write(line)
    pdb.close()
