#!/usr/bin/env python

import sys,os,math

def sub(X,Y):
   return X[0]-Y[0], X[1]-Y[1], X[2]-Y[2]

def cross(X,Y):
   return X[1]*Y[2]-X[2]*Y[1], X[2]*Y[0]-X[0]*Y[2], X[0]*Y[1]-X[1]*Y[0]


def norm(X):
   s = math.sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2])
   return X[0]/s,X[1]/s,X[2]/s

a = None
b = None
c = None

for l in open(sys.argv[1]).xreadlines():
   if not l.startswith("ATOM "): continue
   # print "'"+l[30:38]+"'", "'"+l[38:46]+"'", "'"+l[46:54]+"'"
   if   a is None: a = float(l[30:38]),float(l[38:46]),float(l[46:54])
   elif b is None: b = float(l[30:38]),float(l[38:46]),float(l[46:54])
   elif c is None: c = float(l[30:38]),float(l[38:46]),float(l[46:54])
   else: break
   
X = norm(sub(b,a))
Y = norm(cross(X,sub(c,b)))
C = a
print X[0],X[1],X[2], " ", Y[0],Y[1],Y[2], " ", C[0],C[1],C[2]
