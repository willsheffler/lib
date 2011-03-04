import pymol
from math import *


for rms in range(10):
    d = '/users/sheffler/project/visualize_scores/1ahn/rms'+`rms`+'/'
    files = filter(lambda x: x.endswith('.pdb'),os.listdir(d))
    print files
    for fname in files[:1]:
        print rms,fname
        cmd.load(d+fname)
