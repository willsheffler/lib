import sys,os

for fname in sys.argv[1:]:
    print fname
    f = open(fname,'r+')
    begin = f.tell()
    print begin
    ii = 0
    line = f.readline()
    while line != "":
        end = f.tell()
        if line[30] not in ('E','L','H'):
            f.seek(begin+30)
            f.write('L')
            f.seek(end)
        begin = f.tell()        
        line = f.readline()
        ii += 1
        print ii
        #break
    f.close()
