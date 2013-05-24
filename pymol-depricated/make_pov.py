# make_pov.py                                                                                       

from pymol import cmd

def make_pov(file):
  """
  Do "run make_pov.py" from within pymol and then execute the script                                
  with "make_pov('povray.inp')" to create the povray.inp file.                                      
 
  make note of the dimensions written out by this script and use the
  same (or the same ratio) within povray:
  e.g. output (the "vol" lines only appear if spheres are present):
                                                                                                    
   PyMOL>make_pov('povray.inp')
    RayRenderPOV: w 1100 h 900 f   63.349 b  102.310
    RayRenderPOV: vol  -17.316   17.316  -14.168
    RayRenderPOV: vol   14.168   63.349  102.310
    RayRenderPovRay: processed 714 graphics primitives.
    Ray: total time: 0.03 sec. = 123885.9 frames/hour. (0.03 sec. accum.)
 
  So for example, you would run povray with:
    povray +Ipovray.inp +W1100 +H900
  or
    povray +Ipovray.inp +W550 +H450
  """
  (header,data) = cmd.get_povray()
  povfile=open(file,'w')
  povfile.write(header)
  povfile.write(data)
  povfile.close()

cmd.extend('make_pov',make_pov)
