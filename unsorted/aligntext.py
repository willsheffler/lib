#!/usr/bin/env python
from numpy import *
import sys,re

Aij = 0 # come from A and match (no gap)
Bij = 1 # come from B and match (no gap)
Cij = 2 # come from C and match (no gap)
Ai  = 3 # come from A and gap in j
Bi  = 4 # come from B and gap in j
Ci  = 5 # come from C and gap in j
Aj  = 6 # come from A and gap in i
Bj  = 7 # come from B and gap in i
Cj  = 8 # come from C and gap in i

Infinity     = 999999

def which_min(V):
  """return min in V and index of min"""
  m = Infinity
  a = -1
  for i,v in enumerate(V):
    if v < m:
      m = v
      a = i
  return m,a

#TODO somehow make a starting ' ' equivaluent to a gap!!!!!

def StartGapCost(s):
  if s == ' ':
    return 2
  #  if s1 in ('/', ':', '.', ',', '(', ')'):
  #    return 1
  return 10

def ContinueGapCost(s):
  if s == ' ':
    return 1
  #  if s1 in ('/', ':', '.', ',', '(', ')'):
  #    return 1
  return 2

def MatchCost(s1,s2):
  if s1==s2:
    #if s1 in "(){},;'\"[]#".split():
    #  return -5
    #else:
      return -1
  else:
    return 2

def backtrace(S,i,j,X,Y):
  """output alignment given matrix"""
#  print i,j
  si,sj = '',''
  if   i>0 and S[i,j] == S[i-1,j]+ContinueGapCost(X[i-1]):
    tmpi,tmpj = align(S,i-1, j ,X,Y)
    si = tmpi + X[i-1]
    sj = tmpj + '_'  
#    print "'"+si+"'"    
#    print "'"+sj+"'"    
  elif j>0 and S[i,j] == S[i,j-1]+ContinueGapCost(Y[j-1]):
    tmpi,tmpj = align(S, i ,j-1,X,Y)
    si = tmpi + '_'
    sj = tmpj + Y[j-1]
#    print "'"+si+"'"    
#    print "'"+sj+"'"    
  elif i>0 and j>0 and S[i,j]==S[i-1,j-1]+MatchCost(X[i-1],Y[j-1]):
    tmpi,tmpj = align(S,i-1,j-1,X,Y)
    si = tmpi+X[i-1] 
    sj = tmpj+Y[j-1]
#    print "'"+si+"'"    
#    print "'"+sj+"'"    
  elif i==0 and j==0:
    pass
  else:
    print "ERROR!!!"
#    print i,j
#    print S
  return si,sj

def alignment(X, Y):
  m = len(X)
  n = len(Y)
  S = zeros((m+1,n+1))         # c[i,j] = optimal alignment cost
                             # for X_i and Y_j
  for i in range(m + 1):     # fill in 0th row and column of c
    S[i,0]        = i * ContinueGapCost(X[i-1]) # with gap costs
  for j in range(n + 1):
    S[0,j]        = j * ContinueGapCost(Y[j-1])

  for i in range(1, m + 1):
    for j in range(1, n + 1):
      if X[i - 1] == Y[j - 1]:        # initialize c[i,j] to cost of
         S[i  ,j] = S[i-1,j-1]        # (mis)matching x_i and y_j
      else:                           # plus the optimal cost of 
        mmc = MatchCost(X[i-1],Y[j-1]) 
        S[i,j]    = S[i-1,j-1] + mmc  # aligning X_i-1 and Y_j-1
      gp = ContinueGapCost(X[i-1])
      if gp + S[i-1,j] < S[i,j]:      # check if aligning x_i with a
        S[i,j]    = gp + S[i-1,j]     # gap is better
      else:
        gp = ContinueGapCost(Y[j-1])
        if gp + S[i,j-1] < S[i,j]:      # check if aligning y_j with a
          S[i,j]    = gp + S[i,j-1]     # gap is better

  return S

def affineMatricies(s1,s2,offdiag=Infinity,localstart=False,fixedstart=True):
  m = len(s1)
  n = len(s2)
  A  = zeros((m+1,n+1)) + Infinity
  B  = zeros((m+1,n+1)) + Infinity # gap in i
  C  = zeros((m+1,n+1)) + Infinity # gap in j
  Ta = zeros((m+1,n+1)) + Infinity # traceback A
  Tb = zeros((m+1,n+1)) + Infinity # traceback B
  Tc = zeros((m+1,n+1)) + Infinity # traceback C

  A[0,0] = 0
  B[0,0] = StartGapCost(s1[0])
  C[0,0] = StartGapCost(s2[0])
  for i in range(1,m + 1):     # fill in 0th row and column of c
    A[i,0] = Infinity
    B[i,0] = Infinity
    if localstart:
      C[i,0] = 1 # can start anywhere if local alignment
    elif fixedstart:
      C[i,0] = Infinity
    else:
      C[i,0] = C[i-1,0]+ ContinueGapCost(s1[i-1]) # with gap costs
    Tc[i,0] = Ci
  for j in range(1,n + 1):
    A[0,j] = Infinity
    if localstart:
      B[0,j] = 1
    elif fixedstart:
      B[0,j] = Infinity
    else:
      B[0,j] = B[0,j-1] + ContinueGapCost(s2[j-1])
    C[0,j] = Infinity    
    Tb[0,j] = Bj
    
  for i in range(1, m + 1):
    for j in range(max(1,i-offdiag-1), min(n+1,i+offdiag+2)):
        
      mtc  = MatchCost(s1[i-1],s2[j-1])
      gci = ContinueGapCost(s1[i-1])
      gcj = ContinueGapCost(s2[j-1])
      A[i,j],Ta[i,j] = which_min(( A[i-1,j-1]+mtc, B[i-1,j-1]+mtc, C[i-1,j-1]+mtc, 
                                   Infinity,       B[i-1, j ]+gci, C[i-1, j ]+gci    ,
                                   Infinity,       B[ i ,j-1]+gcj, C[ i ,j-1]+gcj  ) )

                                
      BcostFromA = A[ i ,j-1] + StartGapCost(s2[j-1]) + ContinueGapCost(s2[j-1])
      BcostFromB = B[ i ,j-1] +                     0 + ContinueGapCost(s2[j-1])
      BcostFromC = C[ i ,j-1] + StartGapCost(s2[j-1]) + ContinueGapCost(s2[j-1])

      B [i,j],tmp = which_min((BcostFromA,BcostFromB,BcostFromC))
      Tb[i,j] = 6+tmp #HACK! this depends on ordering of path ENUM
      
      CcostFromA = A[i-1, j ] + StartGapCost(s1[i-1]) + ContinueGapCost(s1[i-1])
      CcostFromB = B[i-1, j ] + StartGapCost(s1[i-1]) + ContinueGapCost(s1[i-1])
      CcostFromC = C[i-1, j ] +                     0 + ContinueGapCost(s1[i-1])

      C[i,j],tmp = which_min((CcostFromA,CcostFromB,CcostFromC))
      Tc[i,j] = 3+tmp #HACK! depends on traceback ENUM

  return A,B,C,Ta,Tb,Tc

def affineAlignment(s1,s2,A,B,C,Ta,Tb,Tc,offdiag=Infinity,localend=True):
  m,n = len(s1),len(s2)
  T = (Ta,Tb,Tc)
    
  mn,tmp = which_min((A[m,n],B[m,n],C[m,n]))
  tmp = T[tmp][m,n]
  arr,step = tmp%3,tmp//3
  i,j = m,n

  if localend:
    aim,bim,cim = min(A[:,n]),min(B[:,n]),min(C[:,n])
    ajm,bjm,cjm = min(A[m,:]),min(B[m,:]),min(C[m,:])
    jm = min(ajm,bjm,cjm)
    im = min(aim,bim,cim)
    # print aim,bim,cim,im
    # print ajm,bjm,cjm,jm
    if jm < im: # gap in i at end
      if   ajm == jm:
        S,j = which_min(A[m,:])
        arr = 0
      elif bjm == jm:
        S,j = which_min(B[m,:])
        arr = 1
      elif cjm == jm:
        S,j = which_min(C[m,:])
        arr = 2
    else: # gap in j at end
      if   aim == im:
        S,i = which_min(A[:,n])
        arr = 0
      elif bim == im:
        S,i = which_min(B[:,n])
        arr = 1
      elif cim == im:
        S,i = which_min(C[:,n])
        arr = 2
    tmp = T[arr][i,j]
    arr,step = tmp%3,tmp//3

  Pi,Pj = s1[i:m],s2[j:n]
  Si,Sj = [],[]
    
  while i>0 or j>0:
#    print i,j,arr,step,(A,B,C)[arr][i,j]
#    print abs(i-j)
#    if abs(i-j) >= offdiag:
#      print "OUT OF BOUNDS!"  
#      return None
    if step == 0 or step == 1:
      Si.append(s1[i-1])
      i -= 1
    else:
      Si.append(' ')
    if step == 0 or step == 2:
      Sj.append(s2[j-1])
      j -= 1
    else:
      Sj.append(' ')
    # print "'"+`arr`+"'"
    tmp = T[int(arr)][i,j]
    arr,step = tmp%3,tmp//3

  Si.reverse()
  Sj.reverse()
  return ''.join(Si)+Pi,''.join(Sj)+Pj
  
  

def align2D(seqs,localstart=False,localend=True,offdiag=Infinity):
  if len(seqs) < 2:
    return seqs
  if len(seqs) is 2:
    s1 = seqs[0]
    s2 = seqs[1]
    A,B,C,Ta,Tb,Tc    = affineMatricies(s1,s2,offdiag=offdiag,localstart=localstart)
    # print A
    # print
    # print B
    # print
    # print C
    # print
    # print Ta
    # print 
    # print Tb
    # print 
    # print Tc
    # print
    aligned1,aligned2 = affineAlignment(s1,s2,A,B,C,Ta,Tb,Tc,offdiag=offdiag,localend=localend)
  
  return (aligned1,aligned2)


def alignMultiple(strs,refstr=0,localstart=False,localend=False,offdiag=Infinity):
  if type(strs) is type(''):
    strs = strs.split("\n")
  indent = len(re.findall("^(\s*)\S",strs[refstr])[0])
  strs = [s.strip() for s in strs]
  refstr = strs[refstr]
  newstrs = [refstr]
  # for ii,s in enumerate(strs[1:]):
  #   refstr,other = align2D( (refstr,s) ,localstart=localstart,localend=localend,offdiag=Infinity )
  #   newstrs[0] = refstr
  #   newstrs.append(other)
  
  for ii in range(1,len(strs)):
    strs[ii-1:ii+1] = align2D( (strs[ii-1],strs[ii]) ,localstart=localstart,localend=localend,offdiag=Infinity )

  r =  range(1,len(strs)); 
  r.reverse()
  for ii in r:
    strs[ii-1:ii+1] = align2D( (strs[ii-1],strs[ii]) ,localstart=localstart,localend=localend,offdiag=Infinity )

  return [" "*indent + s for s in strs]

def main():

  
  teststr = """  rosetta::conformation::protein::ProteinKey const & rosetta::conformation::protein::Protein::key() const [member function]
  rosetta::conformation::MoleculeKey const & rosetta::conformation::Molecule::key() const [member function]
  rosetta::conformation::PolymerKey const & rosetta::conformation::Polymer::key() const [member function]
  rosetta::conformation::StructureKey const & rosetta::conformation::Structure::key() const [member function]
  """
  s = teststr
  # s = sys.stdin.read()
  s = [ x for x in s.split("\n") if x.strip() ]
  if len(s) > 10: sys.exit()
  print "\n".join( alignMultiple( s ) )

  sys.exit()
  
  input = sys.stdin.read().split('\n')
  if input[-1] == '':
    input = input[:-1]
  s1 = input[0]
  s2 = input[-1]
  #print s1 , s2
  indent1,indent2 = 0,0
  if not s1.startswith(' '):
    indent1 = len(re.findall("^(\s*)\S",s1)[0])
  if not s2.startswith(' '):
    indent2 = len(re.findall("^(\s*)\S",s2)[0])
  s1 = "#"*indent1+re.sub(' ',' ',s1[indent1:])
  s2 = "#"*indent2+re.sub(' ',' ',s2[indent2:])
  #sys.exit()
  localstart = False
  localend =   False
  if abs(len(s1)-len(s2)) > 30:
    localstart = True
    localend   = True
  a1,a2 = align2D( (s1,s2), localstart=localstart,localend=localend,offdiag=Infinity )
  print " "*indent1+a1[indent1:]
  for i in input[1:-1]:
    if i is not "":
      print i
  print " "*indent2+a2[indent2:]
    
  if False:#__name__ == "__main__":
  
    seqs = [ "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:703: error: candidates are: numeric::xyzVector< <template-parameter-1-1> >& numeric::xyzVector< <template-parameter-1-1> >::normalize_or_zero() [with T = float]",
        "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:810: error:                 numeric::xyzVector< <template-parameter-1-1> >& numeric::xyzVector< <template-parameter-1-1> >::normalize_or_zero(const T&) [with T = float]",
        "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:689: error: candidates are: void numeric::xyzVector< <template-parameter-1-1> >::normalized(numeric::xyzVector< <template-parameter-1-1> >&) const [with T = float]", 
        "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:796: error:                 void numeric::xyzVector< <template-parameter-1-1> >::normalized(const T&, numeric::xyzVector< <template-parameter-1-1> >&) const [with T = float]",
        "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:888: error:                 numeric::xyzVector< <template-parameter-1-1> > numeric::xyzVector< <template-parameter-1-1> >::normalized() const [with T = float]",    # 
            "/Users/sheffler/svn/branches/python-bindings/librosetta/src/numeric/xyzVector.hh:930: error:                 numeric::xyzVector< <template-parameter-1-1> > numeric::xyzVector< <template-parameter-1-1> >::normalized(const T&) const [with T = float]", 
      "examples.cc:161: error: conflicting declaration 'numeric::xyzVector_float& (numeric::xyzVector<float>::* xyzVector_float_normalized_1)(numeric::xyzVector_float&)'" ]
  
  #  seqs = ('bc','ab')
  
  #  seqs = ["lized(const T&)",
  #          "lized() const"]
  
    for i in range(len(seqs)):
      for j in range(i):
        print i,j      
        a1,a2 = align2D( (seqs[i],seqs[j]), localstart=False,localend=False,offdiag=Infinity )
        # print '"'+a1+'"'
        #   print '"'+a2+'"'
        print a1
        print a2
        print "done",i,j
      
if __name__ == '__main__':
  main()