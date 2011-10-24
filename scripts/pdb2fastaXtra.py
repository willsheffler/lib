#!/usr/bin/python

import string,os
from sys import argv,stdout
from os import popen,system
from os.path import exists
from amino_acids import longer_names


pdbnames = argv[1:]

#chainid = ' '
#if len(argv)>2:
#    chainid = argv[2]

for pdbname in pdbnames:
	if pdbname.startswith("-"): continue # skip args
	
	if not pdbname.endswith('.pdb') and not pdbname.endswith(".gz"):
		pdbname += '.pdb'

	outfile = pdbname

	removechain = 0
	if argv.count('-nochain'):
		removechain = 1

	netpdbname = pdbname
	assert( exists(netpdbname))
	#print 'Reading ... '+netpdbname

	lines = []
	if netpdbname.endswith(".gz"):
		lines = popen("zcat "+netpdbname).readlines()
	else:
		lines = open(netpdbname,'r').readlines()

	oldresnum = '   '
	count = 0;
	extra = dict()
	orig = ""
	mod = ""
	lastnum = -9999
	for line in lines:
		if line.startswith("REMARK 465"):
			aa = line[15:18]
			try: resi = int(line[22:27])
			except ValueError: continue
			extra[resi] = longer_names[aa.replace("MSE","MET")]
		if (len(line)>20): # and (chainid == line[21]):
			line_edit = line
			if line.startswith('TER'):
				if argv.count('-allchains'): fastaid.write('_')
				else: break
			elif line.startswith("ENDMDL"):
				break
			elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
				line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
				if (line_edit[12:14] == 'SE'):
					line_edit = line_edit[0:12]+' S'+line_edit[14:]
				if len(line_edit)>75:
					if (line_edit[76:78] == 'SE'):
						line_edit = line_edit[0:76]+' S'+line_edit[78:]
			if line_edit[0:4] == 'ATOM':
				resnum = line_edit[23:26]
				if not resnum == oldresnum:
					for k in sorted(extra.keys()):
						if lastnum < k < int(resnum):
							mod  += extra[k]
							orig += " "							
					count = count + 1
					longname = line_edit[17:20]
					if longer_names.has_key(longname):
						mod  += longer_names[longname]
						orig += longer_names[longname]
					else:
						mod  += 'X'
						orig += 'X'
					oldresnum = resnum
					lastnum = int(resnum)
				newnum = '%3d' % count
				line_edit = line_edit[0:23] + newnum + line_edit[26:]
				if removechain:
					line_edit = line_edit[0:21]+' '+line_edit[22:]
	for k in sorted(extra.keys()):
		if count < k:
			mod += extra[k]
	if True:#mod != orig:
		if argv.count('-list'): print os.path.basename(netpdbname).replace(".gz","").replace(".pdb","")+": ",
		else: print '>'+pdbname+'\n'
		print mod
		print orig

	#outid.close()
	#fastaid.close()
