#!/usr/bin/env python

import sys,os,time
from optparse import OptionParser
from threading import Thread

parser = OptionParser()
parser.add_option("-t","--threads",dest="num_thread",help="number of threads to spawn",default=1,metavar='INT')
parser.add_option
opts,args = parser.parse_args()

MY_MAX_THREADS = int(opts.num_thread)
MY_NUM_THREADS = int(0)

class cmdThread(Thread):
  def __init__(self,cmd):
    self.cmd = cmd
    Thread.__init__(self)
  
  def run(self):
    global MY_NUM_THREADS
    os.system(self.cmd)
    MY_NUM_THREADS -= 1
  

def main():
  cmds = open(sys.argv[-1]).readlines()
  pool = map(cmdThread,cmds)
  pool.reverse()
  while pool:
    global MY_NUM_THREADS
    if MY_NUM_THREADS >= MY_MAX_THREADS:
      time.sleep(1)
    else: 
      pool.pop().start()
      MY_NUM_THREADS += 1
  
if __name__ == '__main__':
  main()


