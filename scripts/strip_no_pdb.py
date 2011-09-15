#!/usr/bin/python
import sys,os

for i in sys.argv[1:]:
    os.system("grep pdb "+i+" > tmppythonfileDONTMOVEME")
    os.rename("tmppythonfileDONTMOVEME",i)
