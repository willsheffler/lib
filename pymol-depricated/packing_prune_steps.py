import pymol
from pymol import cmd
import sys,os,random
from math import floor

cmd.delete('all')

#loadPackingPDB("/Users/sheffler/project/newtest/prune0.pdb","prune")
#for ii in range(1,9):
#    cmd.load("/Users/sheffler/project/newtest/prune"+`ii`+".pdb","prune")

for ii in range(8):
    loadPackingPDB("/Users/sheffler/project/newtest/prune"+`ii`+".pdb","p"+`ii`)

cmd.hide('all')
cmd.show('spheres')
cmd.select("proteins","resi 0-499")
cmd.select("cavities","resi 500-999")
cmd.color('red','cavities')
cmd.color('white','proteins')
