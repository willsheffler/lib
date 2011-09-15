#!/usr/bin/python
import random,sys

list = [x.strip() for x in open(sys.argv[1]).readlines()]
random.shuffle(list)
for l in list:
    print l
