from pymol import cmd
import glob

def load_models(files,obj,discrete=0):
  """
  load_models <files>, <object>, <discrete=0>

  loads multiple files (using filename globbing)
  into a single object (e.g. from modelling or NMR).

  use discrete=1 if you want to color individual states separately

  e.g. load_models prot_*.pdb, prot
  """
  file_list = glob.glob(files)
  if file_list:
    file_list.sort()
    for name in file_list:
      cmd.load(name,obj,discrete=discrete)
  else:
    print "No files found for pattern %s" % files

cmd.extend('load_models',load_models)
