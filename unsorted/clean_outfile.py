#!/usr/bin/python

import string
from sys import argv,stderr,exit
from os import popen,system
from os.path import basename,exists
from glob import glob

def Help():
    print
    print 'Usage: clean_outfile <file1> <file2> ... [-all] [-fixtag]'
    print '  -all     Get decoys marked with F_ (failed filters)'
    print '  -fixtag  Make decoy tag match its scoreline.'
    print
    print
    exit()

if len(argv)<2:
    Help()


fix_tag = 0
if argv.count('-fixtag'):
    pos = argv.index('-fixtag')
    del(argv[pos])
    fix_tag = 1


split_files = 0
if argv.count('-split'):
    pos = argv.index('-split')
    del(argv[pos])
    split_files = 1
num_file = 1

leave_in_dir = 0
if argv.count('-leave_in_dir'):
    pos = argv.index('-leave_in_dir')
    del(argv[pos])
    leave_in_dir = 1


score_cut_defined = 0
if argv.count('-score'):
    pos = argv.index('-score')
    del(argv[pos])
    score_cut = float( argv[pos])
    del(argv[pos])
    score_cut_defined = 1

alldecoys = 0
if argv.count('-all'):
    pos = argv.index('-all')
    del(argv[pos])
    alldecoys = 1

max_decoys = 99999999999999
if argv.count('-max_decoys'):
    pos = argv.index('-max_decoys')
    max_decoys = int(argv[pos+1])
    del(argv[pos])

figure_out_index_from_file = 0
try:
    S_index = int( argv[-1] ) # Is it an integer?
    del(argv[-1])
except:
    figure_out_index_from_file = 1


infiles_input = argv[1:]

infiles = []
for infile in infiles_input:
    if exists(infile):
        infiles.append(infile)
    else:
        glob_infiles = glob('/net/boinc/results/'+infile[:8]+'/'+infile)
        glob_infiles.sort()
        if len(glob_infiles) == 0:
            glob_infiles = glob('/net/boinc/results_ralph/'+infile[:8]+'/'+infile)
            glob_infiles.sort()

        infiles += glob_infiles


for infile in infiles:
    if not exists(infile):
        Help()


    if infile[-4:] == '.bz2' :
        inp = popen('bzcat '+infile)
        infile = infile[:-4]
    else:
        inp = open(infile,'r')


    outfile =  infile.replace( '.out', '.clean.out' )
    outfile =  outfile.replace( '.sc', '.clean.sc' )
    if alldecoys:
        outfile = infile.replace( '.out', '.clean.all.out' )
    if score_cut_defined:
        outfile = infile.replace( '.out', '.scorecut.clean.out' )

    if not leave_in_dir: outfile = basename(outfile)

    if split_files:
        outfile = outfile.replace('.out','.1.out');

    print 'Cleaning silent file ==> ',outfile
    out = open(outfile,'w')

    writeout = 0
    headerlines = []
    line = inp.readline()
    headerlines += line
    line = inp.readline()
    headerlines += line
    for headerline in headerlines:
        out.write(headerline)


    count = 0
    line = inp.readline()
    if figure_out_index_from_file:
        S_index = line.find('S_')
        if (S_index<0):
            S_index = line.find('F_')
        print "Looked in SCORE tags, and found index of description: ",S_index
        if S_index == -1:
            exit()
    while line:
        if line[0:6] == 'SCORE:':
            if score_cut_defined:
                try:
                    made_the_score_cut = 0
                    score = float( string.split(line)[1] )
                    if (score > score_cut):
                        made_the_score_cut = 1
                except:
                    made_the_score_cut = 0

            if len(line) > S_index and \
                   (not score_cut_defined or made_the_score_cut) and \
                   (line[S_index] == 'S' or \
                    (line[S_index] == 'F' and alldecoys)):

                writeout = 1
                count += 1

                if count > max_decoys:
                    out.close()
                    break

                if (split_files and count % 5000 == 0): #Time for a new file.
                    out.close();
                    command = 'gzip -rf '+outfile
                    print(command)
                    system(command)

                    num_file += 1;
                    outfile = outfile.replace( '.%d.out' % (num_file-1),  '.%d.out' % num_file)
                    print 'Cleaning silent file ==> ',outfile
                    out = open( outfile, 'w')
                    for headerline in headerlines:
                        out.write(headerline)

                if fix_tag:
                    tag = string.split(line)[-1]
            else:
                writeout = 0

        if writeout:
            if fix_tag:
                tag_index = line.find('S_')
                if tag_index > 0:
                    line = line[:tag_index]+tag+'\n'
            out.write(line)

        line = inp.readline()

    if split_files:
        command = 'gzip -rf '+outfile
        print(command)
        system(command)

    print 'Found ', count, ' decoys.'
