#!/usr/bin/python

import sys,os,time

def make_header(fname,brief=None,author='Will Sheffler (willsheffler@gmail.com)'):
  tplt = """// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   %(fname)s
/// @brief  %(brief)s
/// @author %(author)s
/// @date   %(date)s
///
"""
  date = time.ctime()
  if not brief: brief = fname
  return tplt%vars() + "%(ifndef)s\n\n%(nsbegin)s\n\n%(BODY)s\n\n%(nsend)s\n\n"
  
def get_ifndef(fname,fwd=False):
  ifndef = """
#ifndef %(symbol)s_HH
#define %(symbol)s_HH
"""
  symbol = fname.replace("/","_")
  if fwd: symbol += '_FWD'
  return ifndef%vars()
  
def main():
  """make a .hh and .cc file. assumes is a class if beings in cap letter"""
  miniloc = "./"
  if "MINIROSETTA" in os.environ: miniloc = os.environ['MINIROSETTA']+'/src/'
  fname = sys.argv[1]
  namespaces = fname.split("/")[:-1]
  brief = " ".join(sys.argv[2:])
  nsbegin,nsend = "",""
  for n in namespaces: nsbegin += "namespace "+n+" {\n"
  for n in reversed(namespaces): nsend += "} // end namespace "+n+'\n'
  HHBODY = ""
  CCBODY = ""
  cname = fname.split('/')[-1]
  if cname[0] == cname[0].upper() and not cname.count("."):
    HHBODY = """class %(cname)s {\n\npublic:\n\n  %(cname)s();\n\n}; // end class %(cname)s\n\n"""%vars()
    CCBODY = """%(cname)s::%(cname)s( ) {\n  \n}"""%vars()
    fwd = open(miniloc+fname+'.fwd.hh','w')
    BODY = "// forward\nclass "+cname+";\n"
    ifndef = get_ifndef(fname,fwd=True)
    fwd.write(make_header(fname,brief)%vars()+"\n\n#endif\n\n")
    fwd.close()
  ifndef = get_ifndef(fname)
  BODY = HHBODY
  hh = open(miniloc+fname+".hh",'w')
  hh.write(make_header(fname+".hh",brief)%vars()+"\n\n#endif\n\n")
  hh.close()
  BODY = CCBODY
  ifndef = '\n#include "'+fname+'.hh"\n'
  cc = open(miniloc+fname+".cc",'w')
  cc.write(make_header(fname+".cc",brief)%vars())  
  cc.close()
  
if __name__ == '__main__':
  main()