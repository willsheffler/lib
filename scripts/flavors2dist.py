#!/usr/bin/python

import sys,os,re

lines = open(sys.argv[1]).readlines()
features = {}
for l in lines:
    m = re.match(".*VALUES.\s+(\S+).*\[(.*)\].*NAT DIST\s+(\S+)\s+SIZE.\s+(\S+).*",l)
    feat = [int(s[s.find(':')+1:]) for s in m.group(1).split(',')]
    feat.sort()
    #print
    #print m.group(2)
    flav = [float(t) for s in m.group(2).split(',') for t in s.strip()[1:-1].split()]
    #print flav
    if not tuple(feat) in features.keys():
        features[tuple(feat)] = []
    features[tuple(feat)].append(flav)
    
for k in features.keys():
    l = len(features[k])
    for ii in range(l):
        for jj in range(l):
            f1 = features[k][ii]
            f2 = features[k][jj]
            d = 0
            for kk in range(len(f1)):
                v1 = f1[kk]
                v2 = f2[kk]
                while v1 < 0: v1 += 360
                while v2 < 0: v2 += 360
                while v1 >= 360: v1 -= 360
                while v2 >= 360: v2 -= 360
                d += min(abs(v1-v2),360-abs(v1-v2))
            d = d/len(f1)
            print '"'+`k`+'"',ii+1,jj+1,d
