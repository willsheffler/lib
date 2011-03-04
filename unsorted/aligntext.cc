#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cstdlib>
#include <iostream>
using namespace boost::python;
using namespace std;

#define Infinity 999999999

struct Grid {
  int *array;
  int const dim1;
  int const dim2;
  Grid(int const d1, int const d2,int const init=0) : dim1(d1), dim2(d2) {
    array = new int[d1*d2];
    for(int i = 0; i < d1*d2; ++i)
      array[i] = init;
  }
  inline int & operator()(int const i, int const j) {
//    assert i > 0 && j > 0
    return array[i*dim1+j];
  }
  inline int & operator[](int const i) {
    return array[i];
  }
  inline void operator=(int const set_to_this) {
    for( int ii = 0; ii < dim1*dim2; ii += 1 ) {
      array[ii] = set_to_this;
    }
  }
};


inline int StartGapCost(char const s) {
  if( s == ' ')
    return 2;
//  if s1 in ('/', ':', '.', ',', '(', ')'):
  //    return 1
  return 10;
}

inline int ContinueGapCost(char const s){
  if (s == ' ')
    return 1;
  return 2;
}

inline int MatchCost(char const s1, char const s2) {
  if( s1==s2 )
//    if s1 in "(){},;'\"[]#".split()
//      return -5;
//    else
      return -1;
  else
    return 2;
}


int const Aij = 0; // come from A and match (no gap)
int const Bij = 1; // come from B and match (no gap)
int const Cij = 2; // come from C and match (no gap)
int const Ai  = 3; // come from A and gap in j
int const Bi  = 4; // come from B and gap in j
int const Ci  = 5; // come from C and gap in j
int const Aj  = 6; // come from A and gap in i
int const Bj  = 7; // come from B and gap in i
int const Cj  = 8; // come from C and gap in i

inline int min      ( int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8 ) {

  return min(min(min(min(i0,i1),min(i2,i3))),min(min(min(i4,i5),min(i6,i7)),i8));


  return min(
             min(
                 min(
                     min(
                         i0
                         ,i1
                         )
                     ,min(
                         i2
                         ,i3
                         )
                     )
                 )
             ,min(
                 min(
                     min(
                         i4
                         ,i5
                         )
                     ,min(
                         i6
                         ,i7
                         )
                     )
                 )
             ,i8
             );
  

}

//print '\n'.join(["  if( i%i <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i%i;"%(i,i) for i in range(9)])
inline int which_min( int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8 ) {
  if( i0 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i0;
  if( i1 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i1;
  if( i2 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i2;
  if( i3 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i3;
  if( i4 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i4;
  if( i5 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i5;
  if( i6 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i6;
  if( i7 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i7;
  if( i8 <= min( i0,i1,i2,i3,i4,i5,i6,i7,i8 )) return i8;  
}

void affineMatricies( string const s1
                    , string const s2
                    ,  int const offdiag=999999999
                    , bool const localstart=false
                    , bool const fixedstart=false
                    , Grid &  A, Grid &  B, Grid &  C
                    , Grid & Ta, Grid & Tb, Grid & Tc) 
{
  m = s1.size();
  n = s2.size();
  A  = Infinity;
  B  = Infinity; // # gap in i
  C  = Infinity; // # gap in j
  Ta = Infinity; // # traceback A
  Tb = Infinity; // # traceback B
  Tc = Infinity; // # traceback C

  A[0,0] = 0
  B[0,0] = StartGapCost(s1[0])
  C[0,0] = StartGapCost(s2[0])
  for( int ii = 1; ii <= m; ii += 1 ) { // fill in 0th row and column of c
    A[ii,0] = Infinity;
    B[ii,0] = Infinity;
    if( localstart )
      C[ii,0] = 1; // can start anywhere if local alignment
    else if( fixedstart )
      C[ii,0] = Infinity;
    else
      C[ii,0] = C[ii-1,0]+ ContinueGapCost(s1[ii-1]); // # with gap costs
    Tc[ii,0] = Ci;
  }
//  for j in range(1,n + 1):
  for( int ii = 1; ii <= n; ii += 1 ) {
    A[0,j] = Infinity;
    if( localstart )
      B[0,j] = 1;
    else if( fixedstart )
      B[0,j] = Infinity;
    else 
      B[0,j] = B[0,j-1] + ContinueGapCost(s2[j-1]);
    C[0,j] = Infinity;
    Tb[0,j] = Bj;
  }
//  for i in range(1, m + 1):
  for( int ii = 1; ii <= m; ii += 1 ) {
//    for j in range(max(1,i-offdiag-1), min(n+1,i+offdiag+2)):
    for( int jj = max(1,ii-offdiag-1); jj < min(n+1,i+offdiag+2); jj += 1 ) {
      mtc  = MatchCost(s1[i-1],s2[j-1]);
      gci = ContinueGapCost(s1[i-1]);
      gcj = ContinueGapCost(s2[j-1]);
      A[i,j],Ta[i,j] = which_min(  A[i-1,j-1]+mtc , B[i-1,j-1]+mtc , C[i-1,j-1]+mtc 
                                ,  Infinity       , B[i-1, j ]+gci , C[i-1, j ]+gci    
                                ,  Infinity       , B[ i ,j-1]+gcj , C[ i ,j-1]+gcj   
                                );
                                                                
      BcostFromA = A[ i ,j-1] + StartGapCost(s2[j-1]) + ContinueGapCost(s2[j-1]);
      BcostFromB = B[ i ,j-1] +                     0 + ContinueGapCost(s2[j-1]);
      BcostFromC = C[ i ,j-1] + StartGapCost(s2[j-1]) + ContinueGapCost(s2[j-1]);

      B [i,j],tmp = which_min((BcostFromA,BcostFromB,BcostFromC));
      Tb[i,j] = 6+tmp; //HACK! this depends on ordering of path ENUM
      
      CcostFromA = A[i-1, j ] + StartGapCost(s1[i-1]) + ContinueGapCost(s1[i-1]);
      CcostFromB = B[i-1, j ] + StartGapCost(s1[i-1]) + ContinueGapCost(s1[i-1]);
      CcostFromC = C[i-1, j ] +                     0 + ContinueGapCost(s1[i-1]);

      C[i,j],tmp = which_min((CcostFromA,CcostFromB,CcostFromC));
      Tc[i,j] = 3+tmp;// #HACK! depends on traceback ENUM
    }
  }
}

void align2d(string const s1, string const s2) {
  Grid g(s1.size(),s2.size());
  
  cout << "hello" << endl;

}

BOOST_PYTHON_MODULE(hello) {

}

