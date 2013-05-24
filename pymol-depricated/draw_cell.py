#! /usr/bin/python

from cctbx import uctbx, sgtbx
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

#from Numeric import *

def set_to_zero(a):
  if abs(a) < 1e-10:
    a=0
  return a

def draw_cell(obj,radius=0.2):
  """
  From pymol issue the "run draw_cell.py" command to load the script,
  then issue the "draw_cell(object,<optional radius>)" command 
  to actually run it and create the cgo object showing the unit cell
  border for the space group specified by molecular object 'object'.

  e.g. load 1avv.pdb
       run draw_cell.py
       draw_cell 1avv 0.5   (or draw_cell('1avv',.5))

  see also help(draw_cell_param) to draw the cell border for 
  user-defined cell dimensions (i.e. not loaded from a pdb file)

  See also "help(draw_cell_param) to draw the cell border by
  specifying the unit cell parameters directly (i.e. not loaded from
  a pdb file).
  """
  radius=float(radius)
  cell_info=cmd.get_symmetry(obj)
  draw_cell_param(cell_info[0:6],radius)

def draw_cell_param(cell_param_list,radius=0.2):
  """
  If you wish to draw the unit cell border for any cell without the
  need to load a pdb file, then do this:

  e.g. run draw_cell.py
       draw_cell_param((45.2,45.2,70.8,90.,90.,120.),0.5)

  to generate the cell border for this trigonal space group "p 31 2 1"
  with a radius of 0.5A.  Labels for the origin, and A, B and C axes
  will appear as well.  The perimeter of the cell is colored with the
  RGB components corresponding to the A,B,C components.
  """
  
  U=uctbx.unit_cell((cell_param_list))

  vert_000 = map(set_to_zero,U.orthogonalize((0.,0.,0)))
  vert_100 = map(set_to_zero,U.orthogonalize((1.,0.,0)))
  vert_010 = map(set_to_zero,U.orthogonalize((0.,1.,0)))
  vert_001 = map(set_to_zero,U.orthogonalize((0.,0.,1)))
  vert_110 = map(set_to_zero,U.orthogonalize((1.,1.,0)))
  vert_011 = map(set_to_zero,U.orthogonalize((0.,1.,1)))
  vert_101 = map(set_to_zero,U.orthogonalize((1.,0.,1)))
  vert_111 = map(set_to_zero,U.orthogonalize((1.,1.,1)))

#  vert_000 = map(None,U.orthogonalize((0.,0.,0)))
#  vert_100 = map(None,U.orthogonalize((1.,0.,0)))
#  vert_010 = map(None,U.orthogonalize((0.,1.,0)))
#  vert_001 = map(None,U.orthogonalize((0.,0.,1)))
#  vert_110 = map(None,U.orthogonalize((1.,1.,0)))
#  vert_011 = map(None,U.orthogonalize((0.,1.,1)))
#  vert_101 = map(None,U.orthogonalize((1.,0.,1)))
#  vert_111 = map(None,U.orthogonalize((1.,1.,1)))

  #print vert_000

  #CYLINDER = ['CYLINDER']
  #radius = [0.2]
  #print radius
  cell = [] 
  cell.append(CYLINDER)
  cell = cell + vert_000 + vert_100 + [radius] + [0,0,0] + [1,0,0]
  cell.append(CYLINDER)
  cell = cell + vert_000 + vert_010 + [radius] + [0,0,0] + [0,1,0]
  cell.append(CYLINDER)
  cell = cell + vert_000 + vert_001 + [radius] + [0,0,0] + [0,0,1]
  cell.append(CYLINDER)
  cell = cell + vert_100 + vert_110 + [radius] + [1,0,0] + [1,1,0]
  cell.append(CYLINDER)
  cell = cell + vert_100 + vert_101 + [radius] + [1,0,0] + [1,0,1]
  cell.append(CYLINDER)
  cell = cell + vert_010 + vert_110 + [radius] + [0,1,0] + [1,1,0]
  cell.append(CYLINDER)
  cell = cell + vert_010 + vert_011 + [radius] + [0,1,0] + [0,1,1]
  cell.append(CYLINDER)
  cell = cell + vert_001 + vert_101 + [radius] + [0,0,1] + [1,0,1]
  cell.append(CYLINDER)
  cell = cell + vert_001 + vert_011 + [radius] + [0,0,1] + [0,1,1]
  cell.append(CYLINDER)
  cell = cell + vert_110 + vert_111 + [radius] + [1,1,0] + [1,1,1]
  cell.append(CYLINDER)
  cell = cell + vert_101 + vert_111 + [radius] + [1,0,1] + [1,1,1]
  cell.append(CYLINDER)
  cell = cell + vert_011 + vert_111 + [radius] + [0,1,1] + [1,1,1]

  cmd.load_cgo(cell,"cell")
  #return cell

  text = [COLOR, 1.0, 0.0, 1.0,]

  #wire_text(text,plain,[-5.,-5.,-1],'Origin',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((1.05,0.0,0.0))),'A',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((0.0,1.05,0.0))),'B',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((0.0,0.0,1.05))),'C',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])

  cyl_text(text,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[1.0,0.0,1.0])
  cyl_text(text,plain,map(None,U.orthogonalize((1.05,0.0,0.0))),'A',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[1.0,0.0,0.0])
  cyl_text(text,plain,map(None,U.orthogonalize((0.0,1.05,0.0))),'B',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[0.0,1.0,0.0])
  cyl_text(text,plain,map(None,U.orthogonalize((0.0,0.0,1.05))),'C',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[0.0,0.0,1.0])

  cmd.load_cgo(text,'text')

cmd.extend("draw_cell",draw_cell)
