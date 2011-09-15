#!/usr/bin/env python

import sys,os,re
sys.path.append("/work/sheffler/scripts")
from tree_parse import ParseTree,Node,StringNode,GroupNode


def fun_body_re():
  sp0 = r"\s*"
  sp = r"\s+"
  blk = r"(?:%s)"
  const = r"""(?:const)?"""
  virt = r"""(?:virtual)?"""
  typere = blk%r""".*?"""
  name = blk%r"""[a-zA-Z_][-a-zA-Z_0-9]*"""
  args = blk%r"""[^()]*"""
  funre = ".*?" + virt + sp0 + "(?:"+blk%typere+sp+blk%name+"|operator[^()]*?)" + sp + "\(" +  args + "\)" + sp0 + const + sp0 # + blk%(":.*?")+"?"
  r = re.compile(funre,re.VERBOSE|re.M|re.DOTALL)
  return r


#  {\s*return\s.min_type_;\s*}

def process_string_node(node):
  s = node.string
  return s
  

def replace_for_fun_body(nodes,funchk):
  # print "replace_for_fun_body"
  if not nodes: return ""
  assert isinstance(nodes[0],StringNode)
  ret = process_string_node(nodes[0])
  for ii in range(1,len(nodes)):
    # print ii
    if isinstance(nodes[ii],GroupNode) and isinstance(nodes[ii-1],StringNode):
      # print "group preeceeded by string",ii,len(nodes)
      s = nodes[ii-1].string
      s = re.split(r"[;#]",s)[-1]
      s = re.sub(r"//.*|public:|private:|protected:","",s).strip()
      s = re.sub(r"/\*.*?\*/","",s)
      s = s.replace("inline","")
      # if s.count("class"): continue
      # if s.count("namespace"): continue
      s = s.replace("("," ( ")
      s = s.replace(")"," ) ")      
      s = re.sub(r"\s+"," ",s)      
      if funchk(s) and not re.search(r"\b(?:enum|class|struct|namespace)\b",s):
        # print "removing fun body"
        # print s
        # b = nodes[ii+1].get_re()
        # print "------------------------------------------------------"
        # print b
        # print "======================================================"
        ret += ";\n\n"
      else:
        # print "recursing"
        # print nodes[ii]
        # print "======================="
        ret += "\n{\n"+replace_for_fun_body(nodes[ii].children, funchk )+"\n}\n"
    elif isinstance(nodes[ii],GroupNode): # probably shouldn't happen
      assert False
      # print "elif is group",ii,len(nodes)
      # ret += replace_for_fun_body(nodes[ii].children,funchk)
    elif isinstance(nodes[ii],StringNode):  
      # print "elif is string",ii,len(nodes)
      ret += process_string_node(nodes[ii])
    else:
      assert False
  return ret

def strip_fun_body(cppstr):
  funre = fun_body_re()
  cppstr = re.sub("//.*","",cppstr)
  p = ParseTree(cppstr,pairchars=[("{","}")],separators=["NOT_A_SEPARATOR"])
  cppstr = replace_for_fun_body(p.root.children, funre.match )
  r = re.search('\)\s*:.*?[{;]',cppstr,re.DOTALL)
  while r:
    # print "==========="
    # print cppstr[r.start(0):r.end(0)]
    # print "-----"
    cppstr = cppstr[:r.start(0)+1] + cppstr[r.end(0)-1:]
    r = re.search('\)\s*:.*?[;{]',cppstr,re.DOTALL)
  # for b in replace_for_fun_body(p.root.children, funre.match):
  #   # print "======================================"
  #   # print b
  #   l = len(cppstr)
  #   cppstr = re.sub("\s*"+b,";",cppstr)
  #   if len(cppstr) - l > 1000:
  #     print
  #     print b
  #     sys.exit()
  return cppstr


def main():
   fname = sys.argv[1]#"/Users/sheffler/svn/branches/mini/src/core/scoring/EnergyMap.hh"
   o = open(fname)
   s = o.read()
   try:
      print strip_fun_body( s )
   except:
      print s
   o.close()

if __name__ == '__main__':
  main()
