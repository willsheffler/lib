import numpy as np
import itertools as it
import functools as ft
import pandas as pd
import glob,re,os.path,time
from IPython.parallel import Client
from IPython.core.display import HTML
import matplotlib.pyplot as plt
import cgi
import json
import urllib2



def encode_for_url(s):
	# s = s.replace(r'%',r"%25")
	# s = s.replace(r"$",r"%24")
	# s = s.replace(r"&",r"%26")
	# s = s.replace(r"+",r"%2B")
	# s = s.replace(r",",r"%2C")
	# s = s.replace(r"/",r"%2F")
	# s = s.replace(r":",r"%3A")
	# s = s.replace(r";",r"%3B")
	# s = s.replace(r"=",r"%3D")
	# s = s.replace(r"?",r"%3F")
	# s = s.replace(r"@",r"%40")
	s = s.replace(r' ',r"%20")
	s = s.replace(r'"',r"%22")
	# s = s.replace(r'<',r"%3C")
	# s = s.replace(r'>',r"%3E")
	# s = s.replace(r'#',r"%23")
	s = s.replace(r"{",r"%7B")
	s = s.replace(r"}",r"%7D")
	# s = s.replace(r"|",r"%7C")
	# s = s.replace("\\",r"%5C")
	# s = s.replace(r"^",r"%5E")
	# s = s.replace(r"~",r"%7E")
	# s = s.replace(r"[",r"%5B")
	# s = s.replace(r"]",r"%5D")
	# s = s.replace(r"`",r"%60")
	return s


def json_map(t):
	d = dict((kv for kv in it.izip(t.file,t.json)))
	del t['json']
	return d

def read_sicdock_impl(fnames):
	if type(fnames) is type(""): fnames = (fnames,)
	df = None
	pressed = dict()
	for fn in fnames:
		if not os.path.exists(fn):
			print "skipping bad file",fn
			continue
		with open(fn) as f:  buf = f.read()
		buf = re.sub(r"^.*\| ","",buf,flags=re.M)
		compress = None
		if fn.endswith('.gz'): compress = 'gzip'
		if fn.endswith('.bz'): compress = 'bz2'
		t = pd.read_table(fn,delimiter="\s+",verbose=True,compression=compress)
		t['pdb1'] = [x[ 3: 7] for x in t.file]
		t['pdb2'] = [x[13:17] for x in t.file]
		t['order'] = [ int(x.split('_')[-1]) for x in t.file]
		if df: df.append(t)
		else:  df = t
	return df

def read_sicdock(fnames,cache=True):
	if type(fnames) is type(""): fnames = (fnames,)
	df = None
	for fn in fnames:
		if not os.path.exists(fn):
			print "skipping bad file",fn
			continue
		t = None
		storefile = fn+".df"
		if cache and os.path.exists(storefile):
			print "read df"
			#return pd.read_hdf(storefile,'table')
			t = pd.DataFrame.load(storefile)
		else:
			print "read txt"
			t = read_sicdock_impl(fn)
			t.save(storefile)
		if df: df.append(t)
		else:  df = t
	return df

def linkrescore(s,jmap):
	try: j = encode_for_url( jmap[s] )
	except: j = "ERROR"
	url = 'http://localhost:50001/rescore/bouquet='+s+'%20'+j
	return r'<a href="'+url+'">'+s+'</a>'
	# url = 'javascript:void((function(){(new Image()).src = +url+' ;})());'
	s = s.replace("C2_","").replace("C3_","").replace("C4_","").replace("C5_","").replace("C6_","")
	s = s.replace("_1_","_").replace("_2_","_").replace("_3_","_").replace("_4_","_").replace("_5_","_")
	s = s.replace(" "," ").upper()
	# return '<form method="get" action="'+url+'"><button type="submit">'+s+'</button>'
	# </form><FORM><INPUT Type="BUTTON" VALUE="Get '+s+'" ONCLICK="window.open('+url+')"></FORM>'


def showlinks(df,jmap):
	fmtr = ft.partial(linkrescore,jmap=jmap)
	h = df.to_html(formatters=dict(file=fmtr),escape=False)
	return HTML(h)

if __name__ == '__main__':
	df = read_sicdock("/work/sheffler/production/I53/1woz_3l8r.tcdock")
	rescore_dat = json_map(df)
	print df.head(1)
	#showlinks(df.head(2),rescore_dat)
	# jstr = '{"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[78.11381972729157,48.27699558254881,0.0,211.7174744114613,148.2825255885387,174.0],"source":"/gpfs/DS3524-1/WORK/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_1ejb_1.pdb.gz","zero":[0.0,0.0,0.0,121.7174744114613,238.2825255885387,89.99999743867906]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[89.93994039385353,3.552713678800501e-15,34.35400028431229,174.3486645570009,174.3486645570009,28.91321390332912],"source":"/gpfs/DS3524-1/WORK/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1woz_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}'
	# print encode_for_url(jstr)






