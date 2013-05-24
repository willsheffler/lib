import tornado.ioloop
import tornado.web
import json,os,time,threading



WD = "/work/sheffler/tmp/sicdock_service/"
os.chdir(WD)

BIN = "/work/sheffler/bin/sicdock_service"
serv_in,serv_out = os.popen2(BIN)

class RescoreThread(threading.Thread):
	"""docstring for RescoreThread"""
	def __init__(self, cmd):
		super(RescoreThread, self).__init__()
		self.cmd = cmd
	def run(self):
		jstr = """    totsc CLSH __ang1 __ang2  __dis1  __dis2 _closed _intra1 _intra2 _mCOV _mEE _mEXT2   _mRES dN0 dN1 dN2 mDUMP mEXT1 mEXT3  mLOG  mNUM  mRAW  mSQRT nM1 nM2 nT1 nT2 n_A n_E n_H n_L   rmi1  rmi2  rmi3  rmi4 scrgen2  tdis tr0 tr1 tr2  y_dori y_ofst  y_ori
C5_2obx_1_C3_1f23_1_I53R_1            93.815    0   65.0   13.0  -71.20  -62.41 0.98403    0.03    0.00 0.552    0  3.295   94.90 139  14   6  1050   9.6 1.255  23.9    60 112.5   59.8  21  27  41  46  87   0  69  18   73.1  12.5  -9.4  16.4  93.815  30.2 222 150  72  -43.61   3.70  101.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.56244984463271,-37.42965244594426,0.0,211.7174744114613,148.2825255885387,6.999999999999961],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-58.30230741451545,0.0,-22.2694998097998,39.37348783066911,39.37348783066911,117.8943691678887],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
C5_2obx_1_C3_1f23_1_I53R_6            90.808    3   62.0   20.0  -71.33  -62.53 0.98403    0.00    0.00 0.522    0  3.204   90.66 138   8   4   476  10.0 1.416  23.4    56  99.3   44.9  19  28  43  47  90   1  75  14   73.2  12.5  -9.4  16.4  90.808  30.8 222 150  72  -43.70   3.50  101.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.68094666110041,-37.50288750607995,3.552713678800501e-15,211.7174744114613,148.2825255885387,10.00000000000005],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-58.4163820240972,0.0,-22.31307243340665,44.4318126435404,44.4318126435404,122.7645564043373],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
C5_2obx_1_C3_1f23_1_I53R_30           86.288    0   63.0   35.0  -71.11  -71.11 0.98403    0.03    0.00 0.506    0  3.551   85.52 139  10   4   199   9.4 1.332  16.4    58  93.9   52.3  20  25  41  48  89   0  71  18   75.0  11.5  -6.2  14.6  86.288  29.5 222 150  72  -45.57   3.70   90.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.49012670302983,-37.38495428626051,0.0,211.7174744114613,148.2825255885387,9.000000000000036],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-66.42937830811695,0.0,-25.37376466217557,58.14718852103299,58.14718852103299,131.5744637854795],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
"""
		i,o = os.popen(self.cmd)
		i.write(jstr)

class MainHandler(tornado.web.RequestHandler):
	def get(self):
		self.write("top level site")

class RescoreHandler(tornado.web.RequestHandler):
	def get(self,bouquet):
		try:
			bouquet_json = json.loads(bouquet)
			print "write to sicdock_service"
			serv_in.write(bouquet_json)
			print "done write"
			self.write(str(bouquet_json))
		except:
			self.write("bad json:<br>"+bouquet)


application = tornado.web.Application([
	(r"/", MainHandler),
	(r"/rescore/bouquet=(.+)", RescoreHandler),
])

if __name__ == "__main__":
	jstr = """    totsc CLSH __ang1 __ang2  __dis1  __dis2 _closed _intra1 _intra2 _mCOV _mEE _mEXT2   _mRES dN0 dN1 dN2 mDUMP mEXT1 mEXT3  mLOG  mNUM  mRAW  mSQRT nM1 nM2 nT1 nT2 n_A n_E n_H n_L   rmi1  rmi2  rmi3  rmi4 scrgen2  tdis tr0 tr1 tr2  y_dori y_ofst  y_ori
C5_2obx_1_C3_1f23_1_I53R_1            93.815    0   65.0   13.0  -71.20  -62.41 0.98403    0.03    0.00 0.552    0  3.295   94.90 139  14   6  1050   9.6 1.255  23.9    60 112.5   59.8  21  27  41  46  87   0  69  18   73.1  12.5  -9.4  16.4  93.815  30.2 222 150  72  -43.61   3.70  101.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.56244984463271,-37.42965244594426,0.0,211.7174744114613,148.2825255885387,6.999999999999961],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-58.30230741451545,0.0,-22.2694998097998,39.37348783066911,39.37348783066911,117.8943691678887],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
C5_2obx_1_C3_1f23_1_I53R_6            90.808    3   62.0   20.0  -71.33  -62.53 0.98403    0.00    0.00 0.522    0  3.204   90.66 138   8   4   476  10.0 1.416  23.4    56  99.3   44.9  19  28  43  47  90   1  75  14   73.2  12.5  -9.4  16.4  90.808  30.8 222 150  72  -43.70   3.50  101.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.68094666110041,-37.50288750607995,3.552713678800501e-15,211.7174744114613,148.2825255885387,10.00000000000005],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-58.4163820240972,0.0,-22.31307243340665,44.4318126435404,44.4318126435404,122.7645564043373],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
C5_2obx_1_C3_1f23_1_I53R_30           86.288    0   63.0   35.0  -71.11  -71.11 0.98403    0.03    0.00 0.506    0  3.551   85.52 139  10   4   199   9.4 1.332  16.4    58  93.9   52.3  20  25  41  48  89   0  71  18   75.0  11.5  -6.2  14.6  86.288  29.5 222 150  72  -45.57   3.70   90.0 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.49012670302983,-37.38495428626051,0.0,211.7174744114613,148.2825255885387,9.000000000000036],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-66.42937830811695,0.0,-25.37376466217557,58.14718852103299,58.14718852103299,131.5744637854795],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}
"""
	# application.listen(50001)
	# tornado.ioloop.IOLoop.instance().start()

	output = serv_out.readline()
	while not output.startswith("READY"):
		output = serv_out.readline()
		print "server says:",output,

	serv_in.write(jstr)
	serv_in.flush()

	while not os.path.exists("C5_2obx_1_C3_1f23_1_I53R_1.pdb"):
		print "wait for output"
		time.sleep(1)

	print WD+"C5_2obx_1_C3_1f23_1_I53R_1.pdb"


	output = serv_out.readline()
	while output != "RESULT":
		output = serv_out.readline()
		print "server says:",output,

	result = output.split()[1]

	print "got file location",result

	output = serv_out.readline()
	while not output.startswith("READY"):
		output = serv_out.readline()
		print "server says:",output,

	print "DONE"


