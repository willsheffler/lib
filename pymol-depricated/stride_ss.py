#! /usr/bin/python2.2

import os,sys,time

stride_2_pymol = {
        'C':'L',
        'B':'B',
        'b':'B',
        'E':'S',
        'T':'T',
        'G':'G',
        'H':'H',
        }


def stride_extract(pdb_file):
  """
  This function reads the specified PDB format file and runs stride on it.
  (stride must be found in your $PATH).  It extracts the secondary 
  structure data from the standard output of stride and stores it in 
  parallel dictionaries specifying the secondary structure and residue 
  names (both with keys of chain and residue number).
  """
  secondary = {}
  residue = {}

  # run stride (must be on your PATH) and read its output
  stride_cmd = "stride " + pdb_file
  stride_output = os.popen(stride_cmd).readlines()

  for line in stride_output:
    words = line.split()
    if words[0] == "ASG":
      resname = words[1]
      chain = words[2]
      if chain == '-':
        chain = ''
      resnum = int(words[3])
      sec_short = words[5]
      sec_long = words[6]
      if secondary.has_key(chain):
        secondary[chain][resnum] = sec_short
      else:
        secondary[chain] = {resnum:sec_short}
      if residue.has_key(chain):
        residue[chain][resnum] = resname
      else:
        residue[chain] = {resnum:resname}

  return secondary,residue

def stride2pymol(mol='',sel=''):
  """
  usage: stride2pymol object, [selection]

  where the object name is required, but the selection is optional 
  (defaults to all of the object).  
  
  e.g. stride2pymol prot, i. 1-300 and c. a
  """
  #Save the specified molecule (and optional selection sel) to a temporary 
  #PDB format file and call stride_extract on it

  from pymol import cmd

  # map stride's secondary structure alphabet onto pymol's
  # strides definitions:
  # G = 3/10 helix
  # B = Bridge
  # C = coil
  # T = turn
  # H = helix
  # E = strand/sheet

  # create tmpfile from the selection name and the current time in seconds since the epoch
  ext = str(int(time.time()))
  tmpfile=mol + '.pdb_' + ext
  if sel:
    # make sure that selection forces the molecule name to be included in the selection
    # (i.e. is 'and'-ed with the rest of the selection)
    sel = mol + ' and ' + sel
    cmd.save(tmpfile,sel,1,'pdb')
  else:
    cmd.save(tmpfile,mol,1,'pdb')
  # run stride on the saved pdb file and read the standard output to extract the 
  # data to the sec_dict and res_dict dictionaries
  sec_dict,res_dict = stride_extract(tmpfile)
  os.remove(tmpfile)

  # make lists of chains and residues within the chains and sort them
  chain_list = sec_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = sec_dict[chain].keys()
    res_list.sort()
    # now do the alteration of the secondary structure.
    for resnum in res_list:
      cmd.do("alter /%s//%s/%s,ss='%s'" % (mol,chain,resnum,stride_2_pymol[sec_dict[chain][resnum]]))

    # make sure that cartoon automatic is set for this molecule
    cmd.cartoon('automatic',mol)
    # force a rebuild of the drawing
    cmd.rebuild()

###########################################################################################
# for testing purposes:
# if calling this as a program on its own, read the pdb_file name from the
# command line and run stride on it. (does not call stride2pymol, nor require
# importing cmd from pymol
#
if __name__ == '__main__':
  pdb_file = sys.argv[1]
  sec_dict,res_dict = stride_extract(pdb_file)
  chain_list = sec_dict.keys()
  chain_list.sort()
  for chain in chain_list:
    res_list = sec_dict[chain].keys()
    res_list.sort()
    for resnum in res_list:
      print "stride %s %s %s %s stride=%s  pymol_ss='%s'" % (pdb_file,chain,res_dict[chain][resnum],resnum,sec_dict[chain][resnum],stride_2_pymol[sec_dict[chain][resnum]])

# if not calling directly, then extend the def for use in PyMOL
else:
  cmd.extend("stride2pymol",stride2pymol)
