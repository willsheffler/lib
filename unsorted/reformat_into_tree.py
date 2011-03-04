#!/usr/bin/python
import sys,re,string

I = lambda x: x

def INDENTFUNC(d,sep="             "):
  return sep*(d)

def INDENTFUNC2(d):
  return "_INDENT_%i_"%(d)


class DepthTracker:
  
  def __init__(self, pairs = ( ('(',')'), ("[","]"), ("{","}") )  ):
    self.pairs = pairs
    self.pairbegin = tuple( [p[0] for p in pairs ] )
    self.pairend   = tuple( [p[1] for p in pairs ] )
    self.depths = {}
    for p in self.pairbegin:
      self.depths[p] = 0
    self.pairendtobegin = {}
    for p in self.pairs:
      self.pairendtobegin[p[1]] = p[0]
  
  def startPair(self,c):
    if c in self.pairbegin:
      self.depths[c] += 1
      return True
    else:
      return False
  
  def endPair(self,c):
    if c in self.pairend:
      begin = self.pairendtobegin[c]
      assert self.depths[begin] > 0
      self.depths[begin] -= 1
      return True
    else:
      return False
  
  def depth(self):
    return sum(self.depths.values())
  


def flatten(x):
   result = []
   for el in x:
       #if isinstance(el, (list, tuple)):
       if hasattr(el, "__iter__"):
           result.extend(flatten(el))
       else:
           result.append(el)
   return result


class ParseTree(object):
  """holes a tree of Nodes which represent structured text"""
  def __init__(self, stringtoparse=None, pairchars=None, separators=None ):
    self.stringtoparse = stringtoparse
    self.pairchars = pairchars
    if not self.pairchars:
      self.pairchars = (  ('(',')'), ("[","]"), ("{","}"), ("''","''"), ("<",">")  )
    if not separators:
      self.separators = (",","=",":")
    self.pairbegin = tuple( [p[0] for p in self.pairchars ] )
    self.pairend   = tuple( [p[1] for p in self.pairchars ] )
    self.charstacks = {}
    for s,e in self.pairchars:
      self.charstacks[s] = []
    self.pairendtobegin = {}
    for s,e in self.pairchars:
      self.pairendtobegin[e] = s
    self.pairbegintoend = {}
    for s,e in self.pairchars:
      self.pairbegintoend[s] = e
    self.root = None
    if stringtoparse:
      self.parse(stringtoparse)
  
  def partner_of_token(self,t):
    if t in self.pairbegin: return self.pairbegintoend[t]
    if t in self.pairend:   return self.pairendtobegin[t]
    return None
  
  def tokenize(self,to_tokenize):
    tokens = [ "["+ "][".join(list(x)) + "]" for x in flatten(self.pairchars)+list(self.separators)]
    r = re.compile("(%s)"%(  '|'.join(tokens) ))
#    return filter( lambda x: len(x), map( string.strip, r.split(string) ) )
    return filter( I, map( string.strip, r.split(to_tokenize) ) )
  
  def push_node(self,token,node):
    self.charstacks[token].append(node)
  
  def pop_node(self,token):
    return self.charstacks[self.partner_of_token(token)].pop()
  
  def parse(self,string):
    self.root = GroupNode(None)
    isstart = lambda t: t in self.pairbegin
    isend   = lambda t: t in self.pairend    
    istoken = lambda t: isstart(t) or isend(t)
    tokens = self.tokenize(string)
    n = self.root
    while tokens:
      t,tokens = tokens[0],tokens[1:]
      if isstart(t):
        n = GroupNode(n,t,self.partner_of_token(t))
        #self.push_node(t,n)
      elif isend(t):
        n = n.parent
        #self.pop_node(t)
      elif t in self.separators:
        SeparatorNode(n,t)
      else:
        StringNode(n,t)
  
  def __str__(self):
    return str(self.root)


class Node(object):
  """holds one entry in the tree"""
  def __init__(self, parent, FOO=INDENTFUNC2 ):
    self.indentfunc = FOO
    self.parent = parent
    if parent:
      self.parent.add_child(self)
      self.depth = parent.depth + 1
    else:
      self.depth = 0
  
  def indent(self,extra=0):
    return self.indentfunc( self.depth+extra )
  
  def msg(self,*args):
    print self.indent()," ".join(args),`self`
  
  
class StringNode(Node):
  """node that holds a string"""
  def __init__(self, parent, string ):
    self.parent = parent
    self.string = string
    Node.__init__(self,parent)
  
  def __str__(self):
    #self.msg("leaf_node '"+self.string+"' visited!")
    return ""+self.string+""
  

class SeparatorNode(Node):
  """node that holds a string"""
  def __init__(self, parent, sepchar ):
    self.parent = parent
    self.sepchar = sepchar
    Node.__init__(self,parent)
    self.depth = parent.depth
  
  def __str__(self):
    return '\n' +self.indent() + self.sepchar

class GroupNode(Node):
  """internal node holds child strings and groups"""  
  def __init__(self, parent, pairchar='', pairend='', in_between_children="\n", in_between_grouping_chars='\n' ):
    Node.__init__(self,parent)
    self.parent = parent
    self.pairchar = pairchar
    self.pairend  = pairend
    self.children = []
    self.in_between_children = in_between_children
    self.in_between_grouping_chars = in_between_grouping_chars
  
  def add_child(self,node):
    self.children.append(node)
  
  def __str__(self):
    #self.msg("entered:  '"+self.pairchar+"'")
    
    def f(node):
      x = str(node)
#      if type(node) is type(StringNode(None,"")):
#        return x
      if len( re.sub("__INDENT_\d+_","",x.replace(" ","")) ) < 0:
        return re.sub("__INDENT_\d+_","",re.sub("\s+"," ",re.sub("^\s+","",x).replace(",","__COMMA__") ))
      else:
        return x

    child_strs = map( f , self.children )
    joinchar=""
    middle = joinchar.join( child_strs )
    
    #self.msg("finished: '"+self.pairend+"'")
    prefix = self.pairchar + self.indent()
    suffix = '\n' + self.indent()+self.pairend
    
    s = prefix + middle + suffix
    
    # pattern = re.compile( r""" ^                     # line start
    #                            (  \s+  )         # beginning spaces indent
    #                            (  [^\s] [^,]* )  [,]   # start of line up to comma 
    #                            (  .*?   )            # rest of line
    #                            $                    # end of line
    #                        """
    #                      , re.VERBOSE|re.MULTILINE 
    #                      )
    # 
    #s = pattern.sub( r"\1|\2|\3", s )
    return s.replace("__COMMA__",",")
  
  


def process_indents_tabs(s): 
  """docstring for process_indents"""
  for i in range(10):
    s = s.replace("_INDENT_%i_"%i,"        "*i)

def process_indents_clever(s,orig_indent=None): 
  """docstring for process_indents"""
#  return s
  lines = [ (int(x[0]), x[1]) for x in re.findall("_INDENT_(\d+)_(.*?)(?=\n|_INDENT_)",s) ]
  indents = [None]*99
  if orig_indent is None:
    orig_indent = "CODE TREE: "
  prevline = orig_indent
  for indent,line in lines:
    if not indents[indent]:
      if indent > 0:
        indents[indent] = indents[indent-1]
      else:
        indents[indent] = ""
      indents[indent] += " " * len(prevline)
    prevline = line
  r = ""
  prevline = ""
  for indent,line in lines:
    prefix = "\n" + indents[indent]
    if prevline.endswith("("):
      prefix = ""
    l = prefix + line
    r += l
    prevline = l
#  lines = [ re.match("^_INDENT_(\d+)_(.*)$",x).groups() for x in s.split("_INDENT_") ]
#  print lines
  return r


def main():
  """print a pretty formatting of tree-structured text"""
  
  in_string = sys.stdin.read()
  
  orig_indent = re.match("^\s+",in_string)
  if orig_indent:
    orig_indent = orig_indent.group(0)
  else:
    orig_indent = ""
    
  assert in_string.strip().count("\n") == 0

  p = ParseTree(in_string)

  s = str(p).strip()

  s = process_indents_clever(s,orig_indent=orig_indent)

  s = s.replace("__COMMA__","CMA")
  
  print s



       #  print re.sub("(^\s+)([^\s,][^,]*),(.*)$","XXXXXXXXXX", str(p) )

       # pattern = re.compile( r""" ^                     # line start
       #                            (  \s+  )         # beginning spaces indent
       #                            (  [^\s] [^,]* )  [,]   # start of line up to comma 
       #                            (  .*?   )            # rest of line
       #                            $                    # end of line
       #                        """
       #                      , re.VERBOSE|re.MULTILINE 
       #                      )
       # 


# def main():
#   """docstring for main"""
#   if len(sys.argv) > 1:
#     s = sys.argv[1]
#   else:
#     s = sys.stdin.read()
#   pairs = (  ('(',')'), ("[","]"), ("{","}"), ("''","''"), ("<",">")  )  
#   d = DepthTracker(pairs)
#   ii = 0
#   def spacer():
#     return "\t"*d.depth()
#   
#   while ii < len(s):      
#     c = s[ii]
#     if   d.startPair(c):
#       print
#       print spacer(), c#, "start pair"
#       print spacer()+'\t',
#     elif d.endPair(c):
#       print
#       print spacer()+'\t', c#, 'endpair'
#       print spacer(),
#     else:
#       sys.stdout.write(c)      
#     
#     ii += 1
#   print

if __name__ == '__main__':
 main()
