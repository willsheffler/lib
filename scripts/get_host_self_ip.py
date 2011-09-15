#!/usr/bin/python

import os,sys,re

addr = sys.argv[1]

cmd = """ ssh -f %(addr)s "ping -c1 %(addr)s | grep 'bytes from' " """%vars()

pout = os.popen( cmd )
s = pout.read().strip()
pout.close()

#print '"'+s+'"'
ip = re.match(".*\((.+)\).*",s,re.M|re.DOTALL).group(1)
print ip

#ip = {}
#ip['bes' ] = "192.168.78.250"
#ip['niau'] = "192.168.127.71"
#ip['hapy'] = "192.168.133.250"
#ip['yah' ] = "192.168.77.250"
#ip['apep'] = "192.168.132.250"
#ip['gebb'] = "192.168.73.250"
#ip['ptah'] = "192.168.127.66"
#ip['isis'] = "192.168.127.72"
#ip['maat'] = "192.168.69.250"
#
#import sys
#print ip[sys.argv[1]]
#