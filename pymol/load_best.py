#! /usr/bin/python

import os,sys,stats
from pymol import cmd

def load_obj(filename,objname,energy,state,discrete_flag):
    print 'load %s, %s  # %g into state: %d' % (filename, objname, energy,state)
    cmd.load(filename,objname,state,'pdb',discrete=discrete_flag)
    return
  
def load_best(files,obj,number=0,level=0,discrete=0):
  """
  load_best(files,obj,number=0,level=0,discrete=0)

  Loads the best (lowest energy) Modeller files into a single pymol
  molecular object.  You can specify the number of structures to load,
  or an upper energy cutoff (level) below which structures will be
  loaded.

  Set discrete=1 to have the structures loaded with the discrete flag 
  set (allows colouring on individual state properties, for example).

  """
  #input = sys.argv[1]
  grep_cmd = 'grep OBJECTIVE ' + files
  grepped = os.popen(grep_cmd).readlines()

  level=float(level)
  number=int(number)

  entries = {}
  for f in grepped:
    line = f.split(':')
    entries[line[0]] = line[-1][0:-1]

  data = []
  for my_key in entries:
    data.append((my_key,entries[my_key]))

  data.sort(lambda x,y: cmp(x[1],y[1]))

# create array of data to send to statistics calculation
  x=[]
  for y in data:
    x.append(map(float,[y[1]]))

  mean,stdev,median,maximum,minimum = stats.stats(x)


  numcol = len(mean)
  numpts = len(x)

  print "Num pts: ", numpts, "Num cols: ", numcol
  print "col:    Mean:       Stdev:      Median:      Max:        Min:"
  for j in range(numcol):
    print "%4d  %10.8g  %10.8g  %10.8g  %10.8g  %10.8g" % \
    (j, mean[j], stdev[j], median[j],maximum[j], minimum[j])



  count = 0

  if number != 0:
    print "Loading %d lowest energy structures: " % number
    for i in range(number):
      load_obj(data[i][0], obj, float(data[i][1]),i+1,discrete)
      count += 1
  else:
    if level != 0.0:
      print "Loading structures with energy below specified level: ",level
      limit = level
    else:
      print "Loading structures with energy below the median: ",median
      limit = median

    for i in range(len(data)):
      if float(data[i][1]) < limit:
        load_obj(data[i][0], obj, float(data[i][1]),i+1,discrete)
        count += 1

  print 'Loaded %d structures into %s' % (count, obj)


cmd.extend('load_best',load_best)
