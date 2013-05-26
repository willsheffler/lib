#!/usr/bin/env python

import tornado.ioloop
import tornado.web
import json,os,time,threading,subprocess



WD = "/work/sheffler/tmp/sicdock_server/"
os.chdir(WD)

BIN = "/work/sheffler/bin/sicdock_stream"
subproc = subprocess.Popen(BIN, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)

class MainHandler(tornado.web.RequestHandler):
	def get(self):
		self.write("top level site")

class RescoreHandler(tornado.web.RequestHandler):
	def get(self,bouquet_input):
		splt = bouquet_input.split()
		tag = splt[0]
		jbouq = "".join(splt[1:])
		try:
			bouquet_json = json.loads(jbouq)
		except Exception as excpt:
			return self.write("bad bouquet specification:<br>"+jbouq)

		# self.write("tag: "+tag+"<br>")
		try:
			output = subproc.stdout.readline()
			while not output.startswith("READY"):
				output = subproc.stdout.readline()
				print "sicdock_stream says:",output,
		except Exception as excpt:
			return self.write("problem reading from sicdock_stream before write<br>")
		# self.write("about to write to sicdock_stream<br>")
		try:
			subproc.stdin.write(bouquet_input+"\n")
			subproc.stdin.flush()
		except Exception as excpt:
			return self.write("error writing to sicdock_stream:<br>"+str(excpt))

		msg = ""
		try:
			result = subproc.stdout.readline().split()
			assert result[0]=="RESULT"
			pdbfile = result[1]
			pdbfile = os.path.basename(pdbfile)
		except Exception as excpt:
			msg = "problem reading from sicdock_stream subprocess after write<br>"+str(result)

		# newurl = "http://www.google.com/"
		# newurl = "file:///Users/sheffler/tmp/sicdock_stream/C5_2qw7_1_C3_1woz_1_I53F_11.pdb"
		while not os.path.exists(WD+tag+".pdb"): time.sleep(0.1)
		tail = os.popen("tail -n1 "+WD+tag+".pdb").read().strip()
		while not tail=="ENDMDL":
			print tail
			time.sleep(0.1)
			tail = os.popen("tail -n1 "+WD+tag+".pdb").read().strip()

		newurl = "http://localhost:50001/static/"+tag+".pdb"
		self.write('<html><head><meta http-equiv="refresh" content="0; url='+newurl+'"></head><body>'+msg+'<br>forwarding to '+newurl+'</body></html>')

application = tornado.web.Application([
	(r"/", MainHandler),
	(r"/rescore/bouquet=(.+)", RescoreHandler),
	(r'/static/(.+)', tornado.web.StaticFileHandler, {'path':WD}),

])

if __name__ == "__main__":
	try:
		application.listen(50001)
		tornado.ioloop.IOLoop.instance().start()
	finally:
		print "killing sicdock_stream subprocess"
		subproc.kill()

	# jstr = """C5_2obx_1_C3_1f23_1_I53R_1 {"Bouquet":{"stems":{"Stems":[{"Stem":{"Rose":{"position":[-60.56244984463271,-37.42965244594426,0.0,211.7174744114613,148.2825255885387,6.999999999999961],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C5/C5_2obx_1.pdb.gz","zero":[0.0,0.0,0.0,301.7174744114613,121.7174744114613,90.00000256132094]},"SymOper":{"FixedAxisOper":{"axis":[0.8506508083520371,0.5257311121191383,0.0],"nfold":5}}}},{"Stem":{"Rose":{"position":[-58.30230741451545,0.0,-22.2694998097998,39.37348783066911,39.37348783066911,117.8943691678887],"source":"/work/sheffler/data/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/C3_1f23_1.pdb.gz","zero":[0.0,0.0,0.0,90.0,270.0,69.09484062937167]},"SymOper":{"FixedAxisOper":{"axis":[0.9341723589627162,0.0,0.3568220897730885],"nfold":3}}}}]},"sym_gen_depth":4,"sym_gen_dist":9000000000.0}}"""

	# output = subproc.stdout.readline()
	# while not output.startswith("READY"):
	# 	output = subproc.stdout.readline()
	# 	print "server says:",output,

	# subproc.stdin.write(jstr+"\n")
	# subproc.stdin.flush()

	# while not os.path.exists("C5_2obx_1_C3_1f23_1_I53R_1.pdb"):
	# 	print "wait for output"
	# 	time.sleep(1)

	# print WD+"C5_2obx_1_C3_1f23_1_I53R_1.pdb"


	# output = subproc.stdout.readline()
	# while output != "RESULT":
	# 	output = subproc.stdout.readline()
	# 	print "server says:",output,

	# result = output.split()[1]

	# print "got file location",result

	# output = subproc.stdout.readline()
	# while not output.startswith("READY"):
	# 	output = subproc.stdout.readline()
	# 	print "server says:",output,

	# print "DONE"


