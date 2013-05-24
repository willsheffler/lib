#! /usr/bin/python
# Copyright (c) 2003 Robert L. Campbell

"""
data2bfactor: contains the functions 
   data2b_atom(mol='',data_file=''),  
   data2b_res(mol='',data_file=''),
   data2q_atom(mol='',data_file='')  and 
   data2q_res(mol='',data_file='')
"""

import os,sys,time,re

comment = re.compile('^\s*$|^\s*#')

def atom_data_extract(data_file):
  """
  Read the specified 'by-atom' data file and extract the data from it
  and store it in parallel dictionaries specifying the data
  and residue names (both with keys of chain and residue number and atom name).
  The data file can contain comment lines starting with '#' (on lines by themselves).
  These comment lines are ignored.
  """
  bdat = {}
  residue = {}
  chain = ''

  data_lines = file(data_file).readlines()

  for line in data_lines:
# ignore comment lines (beginning with a '#') or blank lines
    if not comment.match(line):
      words = line.split()

# check number of columns of data
      if len(words) == 5:
        chain = words[0]
        resnum = int(words[1])
        resname = words[2]
        atname = words[3]
        if chain == '-':
          chain = ''
        data = float(words[4])
      elif len(words) == 4:
        resnum = int(words[0])
        resname = words[1]
        atname = words[2]
        data = float(words[3])
      else:
        sys.stderr.write("Error in reading data files -- check number of columns")
        sys.exit(1)

      if bdat.has_key(chain):
        if bdat[chain].has_key(resnum):
          bdat[chain][resnum][atname] = data
        else:
          bdat[chain][resnum] = {atname:data}
      else:
        bdat[chain] = {resnum:{atname:data}}

      if residue.has_key(chain):
        if residue[chain].has_key(resnum):
          residue[chain][resnum][atname] = resname
        else:
          residue[chain][resnum] = {atname:resname}
      else:
        residue[chain] = {resnum:{atname:resname}}


  return bdat,residue

def residue_data_extract(data_file):
  """
  Read the specified 'by-residue' data file and extract the data from it
  and store it in parallel dictionaries specifying the data
  and residue names (both with keys of chain and residue number).
  The data file can contain comment lines starting with '#' (on lines by themselves).
  These comment lines are ignored.
  """
  bdat = {}
  residue = {}
  chain = ''

  data_lines = file(data_file).readlines()

  for line in data_lines:
# ignore comment lines (beginning with a '#') or blank lines
    if not comment.match(line):
      words = line.split()

# check number of columns of data
      if len(words) == 4:
        chain = words[0]
        resnum = int(words[1])
        resname = words[2]
        if chain == '-':
          chain = ''
        data = float(words[3])
      elif len(words) == 3:
        resnum = int(words[0])
        resname = words[1]
        data = float(words[2])
      elif len(words) == 2:
        resnum = int(words[0])
        data = float(words[1])
        resname = ''
      else:
        sys.stderr.write("Error in reading data files -- check number of columns")
        sys.exit(1)

      if bdat.has_key(chain):
        bdat[chain][resnum] = data
      else:
        bdat[chain] = {resnum:data}

      if residue.has_key(chain):
        residue[chain][resnum] = resname
      else:
        residue[chain] = {resnum:resname}

  return bdat,residue

def data2b_atom(mol='',data_file=''):
  """
  usage: data2b_atom <mol>, <data_file>

  where <mol> is the molecular object whose B-factor data you wish to modify
  and <data_file> is a file contain the data (one value for each atom)
  The format of <data_file> should be:

     chain resnum resname name data
  or
     resnum resname name data

  (i.e. "chain" is optional if all atoms are in one chain). 
  Lines beginning with '#' are ignored as comments.  
  """
#  call the function to extract the data per atom from 'data_file'
#  alter 'mol' with it.

  from pymol import cmd

  # read the data file and extract the 
  # data to the b_dict and res_dict dictionaries
  b_dict,res_dict = atom_data_extract(data_file)

  # make lists of chains and residues within the chains and sort them
  chain_list = b_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = b_dict[chain].keys()
    res_list.sort()
    for resnum in res_list:
        atom_list = b_dict[chain][resnum].keys()
    # now do the alteration of the B-factor data
        for at in atom_list:
          cmd.do("alter /%s//%s/%s/%s/,b=%f" % (mol,chain,resnum,at,b_dict[chain][resnum][at]))

    # force a rebuild of the drawing
    cmd.rebuild()

def data2b_res(mol='',data_file=''):
  """
  usage: data2b_res <mol>, <data_file>

  where <mol> is the molecular object whose B-factor data you wish to modify
  and <data_file> is a file contain the data (one value for each residue)
  The format of <data_file> should be:

     chain resnum resname data
  or
     resnum resname data

  (i.e. "chain" is optional). Lines beginning with '#' are ignored as comments.  
  """

#call the function to extract the data per residue from 'data_file'
#alter 'mol' with it.

  from pymol import cmd

  # read the data file and extract the 
  # data to the b_dict and res_dict dictionaries
  b_dict,res_dict = residue_data_extract(data_file)

  # make lists of chains and residues within the chains and sort them
  chain_list = b_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = b_dict[chain].keys()
    res_list.sort()
    # now do the alteration of the B-factor data
    for resnum in res_list:
      cmd.do("alter /%s//%s/%s/,b=%f" % (mol,chain,resnum,b_dict[chain][resnum]))

    # force a rebuild of the drawing
    cmd.rebuild()

def data2q_atom(mol='',data_file=''):
  """
  usage: data2q_atom <mol>, <data_file>

  where <mol> is the molecular object whose occupancy data you wish to modify
  and <data_file> is a file contain the data (one value for each atom)
  The format of <data_file> should be:

     chain resnum resname name data
  or
     resnum resname name data

  (i.e. "chain" is optional). Lines beginning with '#' are ignored as comments.  
  """
#  call the function to extract the data per atom from 'data_file'
#  alter 'mol' with it.

  from pymol import cmd

  # read the data file and extract the 
  # data to the q_dict and res_dict dictionaries
  q_dict,res_dict = atom_data_extract(data_file)

  # make lists of chains and residues within the chains and sort them
  chain_list = q_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = q_dict[chain].keys()
    res_list.sort()
    for resnum in res_list:
        atom_list = q_dict[chain][resnum].keys()
    # now do the alteration of the occupancy (q) data
        for at in atom_list:
          cmd.do("alter /%s//%s/%s/%s/,q=%f" % (mol,chain,resnum,at,q_dict[chain][resnum][at]))

    # force a rebuild of the drawing
    cmd.rebuild()

def data2q_res(mol='',data_file=''):
  """
  usage: data2q_res <mol>, <data_file>

  where <mol> is the molecular object whose occupancy data you wish to modify
  and <data_file> is a file contain the data (one value for each residue)
  The format of <data_file> should be:

     chain resnum resname data
  or
     resnum resname data

  (i.e. "chain" is optional). Lines beginning with '#' are ignored as comments.  
  """
#  call the function to extract the data per residue from 'data_file'
#  alter 'mol' with it.

  from pymol import cmd

  # read the data file and extract the 
  # data to the q_dict and res_dict dictionaries
  q_dict,res_dict = residue_data_extract(data_file)

  # make lists of chains and residues within the chains and sort them
  chain_list = q_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = q_dict[chain].keys()
    res_list.sort()
    # now do the alteration of the occupancy (q) data
    for resnum in res_list:
      cmd.do("alter /%s//%s/%s/,q=%f" % (mol,chain,resnum,q_dict[chain][resnum]))

    # force a rebuild of the drawing
    cmd.rebuild()

cmd.extend('data2b_res',data2b_res)
cmd.extend('data2b_atom',data2b_atom)
cmd.extend('data2q_res',data2q_res)
cmd.extend('data2q_atom',data2q_atom)
###########################################################################################
# for testing purposes:
# if calling this as a program on its own, read the pdb_file name from
# the command line and run residue_data_extract on it. (does not require
# importing cmd from pymol

if __name__ == '__main__':
  pdb_file = sys.argv[1]
  b_dict,res_dict = residue_data_extract(pdb_file)
  chain_list = b_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = b_dict[chain].keys()
    res_list.sort()
    for resnum in res_list:
      print "b-factors %s %s %s %s  new B='%s'" % (pdb_file,chain,res_dict[chain][resnum],resnum,b_dict[chain][resnum])


