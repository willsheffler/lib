#!/usr/bin/python
## make mammoth structure alignments


import string
from glob import glob
from sys import argv,stderr,exit
from os import popen,system
from os.path import exists
from operator import add
from math import sqrt

#############################
def Help():
    print '\n'
    print '-'*75
    print 'USAGE: %s <pdb1> <pdb2> {... <pdbN>} > <superposition-pdb>'%argv[0]
    print '\n will superimpose pdbs 2-N onto pdb1 using maxsub, so seqs should agree'
    print '-'*75
    print '\n\n'
    exit()

if len(argv) <=2:
    Help()


RENUMBER_ATOMS = 0
SHOW_MODEL_0 = 1

model_count = 0
atom_count = 0

args = argv[1:]
R_DEFINED = 0
if args.count('-R'):
    pos = args.index('-R')
    rmsd_threshold = float(args[pos+1])
    del args[pos]
    del args[pos]
    distance_threshold = rmsd_threshold # + 3
    R_DEFINED = 1
else:
#    rmsd_threshold = 0.0
    rmsd_threshold = 4.0

if args.count('-D'):
    pos = args.index('-D')
    distance_threshold = float(args[pos+1])
    assert(R_DEFINED)

CALC_PER_RESIDUE_DEVIATIONS = 0
if args.count('-per_res'):
    pos = args.index('-per_res')
    del( args[ pos] )
    CALC_PER_RESIDUE_DEVIATIONS = 1

if R_DEFINED:
    stderr.write( 'using distance threshold %5.1f A; reporting number of residues within %5.1f A \n' % (distance_threshold, rmsd_threshold))

COPY_RESNUM = 0
if args.count('-copy_resnum'):
    del args[args.index('-copy_resnum')]
    COPY_RESNUM = 1


RENUMBER_ATOMS = 0
if args.count('-renumber_atoms'):
    del args[args.index('-renumber_atoms')]
    RENUMBER_ATOMS = 1

COPY_HETATM = 0
if args.count('-copy_hetatm'):
    del args[args.index('-copy_hetatm')]
    COPY_HETATM = 1

if args.count('-N'):
    pos = args.index('-N')
    slicenum = int( args[pos+1])
    del args[pos]
    del args[pos]
    stderr.write( 'using first %d residues\n' % slicenum )
    slice  = 1
else:
    rmsd_threshold = 4.0
    slice = 0


subset_residues = []
use_subset = 0

if args.count('-subset'):
    pos = args.index('-subset')
    del args[pos]
    use_subset = 1
    stderr.write( 'using a subset of residues: '  )
    goodint = 1
    while goodint:

        try:
            subset_start = int(args[pos].split('t')[0])
            subset_end = int(args[pos].split('t')[1])
            del args[pos]

            print subset_start, subset_end, range(subset_end-subset_start)

            for s in range(subset_end-subset_start+1):
                subset_residues.append( subset_start+s )
                #stderr.write('%d ' % subset_start+s )

            #subset_residue = int(args[pos])
            #subset_residues.append( subset_residue )
            #subset_residues.append( subset_residue - 1 )
            #subset_residues.append( subset_residue + 1)
            #stderr.write('%d ' % subset_residue )
        except:
            goodint = 0

    stderr.write( '\n'  )

if args.count('-1'):
    del args[args.index('-1')]
    SHOW_MODEL_0 = 0

pdb_list = args
pdb1 = pdb_list[0]

for pdb in pdb_list:
    if not exists( pdb ):
        stderr.write( 'Could not find '+pdb+'\n' )
        Help()



for pdb in pdb_list[1:]:

    pdb1_to_superimpose = pdb1
    pdb_to_superimpose = pdb

    if slice:
        command = '~rvernon/scripts/pdbslice.py %s 1 %d blah_' % (pdb1, slicenum)
        system(command)
        pdb_to_superimpose = 'blah_'+pdb1

    if use_subset:

        command = '~rvernon/scripts/pdbslice.py '+pdb1
        command += ' -subset '
        for i in subset_residues:
            command += ' %d ' % i
        command += ' blah_ '
        system(command)
        pdb1_to_superimpose = 'blah_'+pdb1

        command = '~rvernon/scripts/pdbslice.py '+pdb
        command += ' -subset '
        for i in subset_residues:
            command += ' %d ' % i
        command += ' blah_ '
        system(command)
        pdb_to_superimpose = 'blah_'+pdb


    if R_DEFINED:
        command = '~rvernon/mammoth2/mammoth_rna -R %f -D %f -p %s -e %s 2> /dev/null | grep PSI.end'\
                          %(rmsd_threshold,distance_threshold,pdb1_to_superimpose,pdb_to_superimpose)

    else:
        command = '~rvernon/mammoth2/mammoth_rna -p %s -e %s 2> /dev/null | grep PSI.end'\
            %(pdb1_to_superimpose,pdb_to_superimpose)
#        command = '/work/pbradley/maxsub/maxsub -p %s -e %s 2> /dev/null | grep PSI.end'\
#                  %(pdb_to_superimpose,pdb)

#    if rmsd_threshold:
#        command = '/work/pbradley/mammoth2/test.out -R %f -p %s -e %s 2> /dev/null | grep PSI.end'\
#              %(rmsd_threshold,pdb1,pdb)
#    else:
#        command = '/work/pbradley/mammoth/mastodon -p %s -e %s 2> /dev/null | grep PSI.end'\
#                  %(pdb1,pdb)

#    stderr.write(command+'\n')
    lines = popen(command).readlines()

    if slice:
        command = 'rm blah_'+pdb1
        system(command)

    if use_subset:
        command = 'rm blah_'+pdb1
        system(command)

    if not lines:
        stderr.write('empty file? %s\n'%pdb)
        continue

    l = string.split(lines[0])

    stderr.write('%s -vs- %s: %s over %s residues\n'%(pdb_list[0],pdb,l[7],l[3]))


    file = 'maxsub_sup.pdb'
    matrix = map(lambda x:map(float,string.split(x)[1:]),
                 popen('grep -A3 "Transformation Matrix" %s'%file).readlines()[1:4])


    P_translation = map(float,
                        string.split(popen('grep -A1 "Translation vector (Pred" %s'\
                                           %file).readlines()[1])[1:])

    E_translation = map(float,
             string.split(popen('grep -A1 "Translation vector (Exp" %s'\
                                %file).readlines()[1])[1:])


    def E_transform(v,matrix,tP,tE):
        ans = [0.0]*3
        for i in range(3):
            for j in range(3):
                ans[i] = ans[i] + matrix[i][j]*(v[j]-tE[j])
            ans[i] = ans[i] + tP[i]

        return ans

    if model_count == 0:
        model0_resnum = []
        hetatm_lines = []

        model0_xyzs = []

        if SHOW_MODEL_0:
            print 'MODEL     %4d'%model_count
            model_count = model_count+1

            data = open(pdb1,'r')
            line = data.readline()

            prev_resnum = ''
            while line:
                if line[:4] in ['ATOM','HETA']:
                    atom_count = atom_count + 1
                    print '%s%5d%s'%(line[:6],atom_count,line[11:-1])

                    if line[12:16]==' CA ' and CALC_PER_RESIDUE_DEVIATIONS: model0_xyzs.append( [float(line[30:38]), float(line[38:46]), float(line[46:54])] )

                    if COPY_RESNUM:
                        resnum = line[22:26]
                        if resnum != prev_resnum:
                            model0_resnum.append( resnum )
                        prev_resnum = resnum

                    if line[:4] == 'HETA':
                        hetatm_lines.append(line[:66])

                elif line[:6] == 'ENDMDL':break
                line = data.readline()
            data.close()
            print 'ENDMDL'
        else:
            model_count = model_count + 1

    print 'MODEL     %4d'%model_count
    model_count = model_count+1
    data = open(pdb,'r')
    line = data.readline()

    prev_resnum = ''
    rescount = -1

    atom_count = 0
    while line:
        if line[:4] in ['ATOM','HETA']:
                atom_count = atom_count + 1
                pos = E_transform(map(float,[line[30:38],line[38:46],line[46:54]]),
                                  matrix,
                                  P_translation,
                                  E_translation)

                current_resnum = line[22:26]
                if prev_resnum != current_resnum:
                    rescount += 1
                prev_resnum = current_resnum
                if COPY_RESNUM:
                    new_resnum = model0_resnum[rescount]
                    line = line[0:22]+new_resnum+line[26:]

                if line[12:16]==' CA ' and CALC_PER_RESIDUE_DEVIATIONS:
                    stderr.write( '%4d %8.4f\n' % (rescount+1, \
                                                       sqrt( ( model0_xyzs[rescount][0] - pos[0] )*( model0_xyzs[rescount][0] - pos[0] ) +
                                                             ( model0_xyzs[rescount][1] - pos[1] )*( model0_xyzs[rescount][1] - pos[1] ) +
                                                             ( model0_xyzs[rescount][2] - pos[2] )*( model0_xyzs[rescount][2] - pos[2] ) ) ) )


                if RENUMBER_ATOMS:
                    print '%s%5d%s%8.3f%8.3f%8.3f%s'\
                          %(line[:6],atom_count,line[11:30],pos[0],pos[1],pos[2],line[54:-1])
                else:
                    print '%s%s%8.3f%8.3f%8.3f%s'\
                          %(line[:6],line[6:30],pos[0],pos[1],pos[2],line[54:-1])

        elif line[:6] == 'ENDMDL':break
        line = data.readline()

    if COPY_HETATM:
        for line in hetatm_lines:
            atom_count += 1
            if RENUMBER_ATOMS:
                print '%s%5d%s' % (line[:6],atom_count,line[11:])
            else:
                print line

    data.close()

    print 'ENDMDL'


system('rm maxsub_sup.pdb maxsub_sup2.pdb rasmol.tcl')
