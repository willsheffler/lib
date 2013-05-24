import pdb
import sys,urllib2,zlib,os,threading,glob,math
from pymol import cmd
from Tkinter import *
sys.path.append("/Users/sheffler/lib/pymol")
from util import Mat,Vec,getres,aa_3_1,aa_1_3,redopent,rot,trans,com,getaa

savedir = "mutalyze_PYMOL_WORKING/"

class MutalyzeGui(Frame):
	def __init__(self,manager,master=None):
		Frame.__init__(self,master)		
		self.pack()
		self.manager = manager
		self.createWidgets(self.manager.d.muts)
	
	def createWidgets(self,d):
		self.bf = Frame(self)
		self.quit = Button(self.bf,text="DONE"       ,command=self.quit       ).pack({"side": "top"})
		self.rsta = Button(self.bf,text="Reset Align",command=self.manager.reset_align).pack({"side": "top"})
		self.zoom = Button(self.bf,text="Zoom"       ,command=self.manager.zoom       ).pack({"side": "top"})
		self.cena = Button(self.bf,text="Center"     ,command=self.manager.cenall     ).pack({"side": "top"})
		self.topr = Button(self.bf,text="To Primary" ,command=self.manager.toprimary  ).pack({"side": "top"})
		self.charge = StringVar()
		self.charge.set("foo")
		self.crgu = Button(self.bf,textvariable=self.charge,command=self.manager.chargeupdate).pack({"side": "top"})
		self.bf.pack({'side':'top'})
		self.mutframe = Frame(master=self)
		for m in d.muts:
			m.gui = self
			f = Frame(self.mutframe)
			e = Entry(f,width=8,textvariable=t)
			e.insert(0,m.resfileline())
			e.pack({'side':'right'})
			lab = m.rnat.resn+`m.rnat.resi`+" to "+m.rdes.resn+`m.rdes.resi`
			b3 = Button(f,text="S"         ,command=m.showpack).pack({'side':'right'})
			b  = Button(f,text=lab,width=15,command=m.focus   ).pack({'side':'right'})
			b2 = Button(f,text="R"         ,command=m.revert  ).pack({'side':'right'})
			f.pack({'side':'top'})
		self.mutframe.pack({"side":"left"})
	


class MutalyzeManager(object):
	_m,_prevm = None,None
	def setm(s,v): s._m,s._prevm=v,s._m
	m     = property(fget=lambda s: s._m,fset=setm)
	prevm = property(fget=lambda s: s._prevm)
	_d,_prevd = None,None
	def setd(s,v): s._d,s._prevd=v,s._d
	d     = property(fget=lambda s: s._d,fset=setd)
	prevd = property(fget=lambda s: s._prevd)
	def __init__(self,master=None):
		self.dir = "./output/"
		files = [f for f in os.listdir(self.dir) if f.endswith('.pdb') or f.endswith('.pdb.gz')]
		while not files:
			self.dir = tkSimpleDialog.askstring("Input","where is your mutalyze dir?")
			files = [f for f in os.listdir(self.dir) if f.endswith('.pdb') or f.endswith('.pdb.gz')]
		self.des = [Design(self.dir+f,self) for f in files]
		self.d = self.des[0]
		self.d.load()
		self.gui = None
		self.toggle = { 'nat':1,'lines':1,'cartoon':1,'all':1,'view':0 }
		cmd.extend("a"      ,lambda a:         self.addaa(a)       )	
		cmd.extend("c"      ,lambda:           self.centeractive() )
		cmd.extend("car"    ,lambda:           self.togglerep("cartoon"))
		cmd.extend("charge" ,lambda:           self.chargeupdate() )
		cmd.extend("g"      ,lambda:           self.rungui()       )
		cmd.extend("h"      ,lambda:           self.help()         )
		cmd.extend("i"      ,lambda:           self.showinfo()     )
		cmd.extend("iall"   ,lambda:           self.showallinfo()  )
		cmd.extend("j"      ,lambda:           self.showallinfo()  )
		cmd.extend("l"      ,lambda:           self.togglerep("lines"))		
		cmd.extend("n"      ,lambda resi=None: self.next(resi)     )
		cmd.extend("note"   ,lambda txt=None:  self.note(txt)      )
		cmd.extend("cnote"  ,lambda:           self.clearnote()    )
		cmd.extend("nd"     ,lambda pid=None:  self.nextdesign(pid))
		cmd.extend("o"      ,lambda:           self.orient()       )
		cmd.extend("p"      ,lambda:           self.prev()         )
		cmd.extend("pd"     ,lambda:           self.prevdesign()   )
		cmd.extend("ps"     ,lambda:           self.toprimary()    )
		cmd.extend("r"      ,lambda:           self.revert()       )
		cmd.extend("ra"     ,lambda:           self.reset_align()  )
		cmd.extend("s"      ,lambda:           self.showpack()     )
		cmd.extend("sav"    ,lambda:           self.save()         )
		cmd.extend("u"      ,lambda:           self.showbup()      )
		cmd.extend("v"      ,lambda:           self.toggleall()    )
		cmd.extend("vn"     ,lambda:           self.togglenat()    )
		cmd.extend("resfile",lambda:           self.printresfile() )
		self.helpstr = """
cmd: a <AA>          add aa to resfile for active res
cmd: c               center active res (toggle)
cmd: car             toggle cartoon
cmd: charge          print total charge (is it correct???)
cmd: g               run gui       
cmd: h               help      
cmd: i               show info for active res (toggle)
cmd: j               show all info for active res (toggle)
cmd: l               toggle lines	
cmd: n <resi or 0>   next design pos, or resi, or first (if 0)
cmd: note <note>     record note for active res
cmd: cnote           clear notes for active res
cmd: nd <ofst=+1>    next design, saves current in mutalyze_PYMOL_WORKING
cmd: o               toggle global view
cmd: p               prev design pos
cmd: pd              prev design   
cmd: ps              change selection to primary subunit
cmd: r               revert active res (toggle)
cmd: ra              reset to global nat alignment
cmd: s               show spheres around active res (toggle)
cmd: u               show bur uns polar (toggle)
cmd: v               view all (toggle); sometimes useful to see info clearly
cmd: vn              view native (toggle)
cmd: resfile         print resfile to screen
		""".strip()
		print 'MutalyzeManager initialized managing %i designs'%len(self.des)
	
	def orient(self):
		if self.toggle['view']:
			cmd.set_view(self.toggle['view'])
			self.toggle['view'] = None
		else:
			self.toggle['view'] = cmd.get_view()
			cmd.orient(self.d.obj+" or pnt*")
	
	def centeractive(self):
		if self.toggle['view']:
			cmd.set_view(self.toggle['view'])
			self.toggle['view'] = None
		else:
			self.toggle['view'] = cmd.get_view()			
			cmd.center(self.m.rdes.sel())
	
	def help(self):
		print self.helpstr
	
	def rungui(self):
		root = Tk()
		self.gui = MutalyzeGui(self,root)
		self.gui.mainloop()
		root.destroy()
	
	def printresfile(self):
		print "####################### BEGIN RESFILE #######################"
		print self.d.resfile()
		print "######################## END RESFILE ########################"
	
	def reset_align(self):
		cmd.align(self.d.muts[0].rnat.obj,self.d.muts[0].rdes.obj)
	
	def showinfo(self):
		self.m.showinfo()
	
	def showallinfo(self):
		self.m.showallinfo()
	
	def toggleall(self):
		if self.toggle['all']: cmd.disable("(not pseudo*)")
		else:                  cmd.enable ("(not pseudo*)")
		self.toggle['all'] = not self.toggle['all']
	
	def togglenat(self):
		if self.toggle['nat']: cmd.disable(self.d.pid)
		else:                  cmd.enable (self.d.pid)
		self.toggle['nat'] = not self.toggle['nat']
	
	def togglerep(self,rep):
		if self.toggle[rep]: cmd.hide(rep)
		else:                cmd.show(rep)
		self.toggle[rep] = not self.toggle[rep]		
	
	def addaa(self,aa):
		self.m.addaa(aa)
	
	def note(self,txt=None):
		if txt is None: self.printresfile()
		else: self.m.notes.append(txt)
	
	def clearnote(self):
		self.m.notes = []
	
	def zoom(self):
		cmd.zoom(self.d.muts[0].rdes.obj+" or pnt*")
	
	def cenall(self):
		cmd.center(self.d.muts[0].rdes.obj+" or pnt*")
	
	def toprimary(self):
		cs = [x[0] for x in getres('sele')]
		if not cs or any([x!=cs[0] for x in cs]): return
		c = cs[0]
		if c == "B": rot("all",Vec(1,0,0),  72)
		if c == "C": rot("all",Vec(1,0,0), 144)
		if c == "D": rot("all",Vec(1,0,0),-144)
		if c == "E": rot("all",Vec(1,0,0), -72)
		r = getres('sele')[0][1]
		cmd.select("chain A and resi %i"%r)
	
	def chargeupdate(self):
		pos = cmd.select("tmpcharge",self.d.obj+" and resn ARG+LYS and name CA")
		neg = cmd.select("tmpcharge",self.d.obj+" and resn ASP+GLU and name CA")
		cmd.delete("tmpcharge")
		c = "Total Charge: %i"%(pos-neg)
		print c
		if self.gui: self.gui.charge.set(c)
	
	def save(self):
		self.d.save(manual=True)
	
	def nextdesign(self,pid=None):
		newd = None
		try:
			self.d.save()
			if pid:
				try: 
					print "SWITCH TO DES",pid				
					newd = self.desmap[pid]
				except KeyError as e:
					print "no design for pdb",pid
					print e
					return
			else:
				try: 
					print "SWITCH TO DES #",self.des.index(self.d)+1+1,"of",len(self.des)
					newd = self.des[self.des.index(self.d)+1]
				except IndexError as e:
					print "FAIL: SWITCH TO DES #",self.des.index(self.d)+1+1,"of",len(self.des)
					print e
					print "DONE??!?!?!"
					return
			if not newd:
				print "NOT SWITCHING DESIGNS!"
				return
			newd.load()
			return
			#print "now active:",self.d		
		except Exception as e:
			print "FAIL ON",newd.fname if newd else "after "+self.d.fname
			print e
			raise
			return
	
	def prevdesign(self):
		self.d.save()
		newd = None
		try: 
			print "SWITCH TO DES #",self.des.index(self.d)-1+1,"of",len(self.des)
			newd = self.des[self.des.index(self.d)-1]
		except IndexError as e:
			print "FAIL: SWITCH TO DES #",self.des.index(self.d)-1+1,"of",len(self.des)
			print e
			return
		if not newd:
			print "NOT SWITCHING DESIGNS!"
			return
		newd.load()
		#print "now active:",self.d		
	
	def next(self,resi=None):
		if resi is not None:
			if resi == '0': 
				self.d.muts[0].focus()
			else:
				try: self.d.mutmap[int(resi)].focus()
				except KeyError as e:
					print "no design pos resi",resi
					print e
					return
		else:
			self.shiftworkingmut(+1)
	
	def prev(self):
		self.shiftworkingmut(-1)
	
	def shiftworkingmut(self,delta):
		try: self.d.muts[self.d.muts.index(self.m)+delta].focus()
		except ValueError: self.d.muts[0].focus()
		except IndexError as e:
			print "DONE!!!"
			return
		print "now active:",self.m
	
	def revert(self):
		self.m.revert()
	
	def showpack(self):
		self.m.showpack()
	
	def showbup(self):
		self.d.showbup()
	


class Residue(object):
	def __init__(self,resi,chain,obj):
		self.resi = resi
		self.chain = chain
		self.obj = obj
		self.resn = cmd.get_model(self.sel()).atom[0].resn
	
	def sel(self):
		return "(resi %s and chain %s and %s)"%(self.resi,self.chain,self.obj)
	
	def __str__(self):
		return "%s/%s/%s/%i"%(self.obj,self.chain,self.resn,self.resi)
	
	def coords(self):
		return (Vec(a.coord) for a in cmd.get_model(self.sel()).atom)
	
	def isdiffenent(a,b):
		return a.resn != b.resn
	


class DesignPos(object):
	def __init__(self,rdes,rnat):
		self.rdes = rdes
		self.rnat = rnat
		self.oldcrd = None
		self.notes = []
		self.aas = [self.aa()]
		self.scores = {}
		self.res_scores = {}
		self.toggle = {'sph':False,'info':False,'infokind':None}
	
	def __str__(self):
		return str(self.rdes)+" (native %s / %i)"%(self.rnat.resn,self.rnat.resi)
	
	def focus(self):
		self.manager.m = self
		if self.manager.prevm: cmd.hide('sticks',self.manager.prevm.rdes.sel())
		if self.manager.prevm: cmd.color('green',self.manager.prevm.rdes.sel()+" and elem C")
		fr = self.rnat.obj+" and name n+ca+c and resi %i-%i"%(self.rnat.resi-10,self.rnat.resi+10)		
		to = self.rdes.obj+" and name n+ca+c and resi %i-%i"%(self.rdes.resi-10,self.rdes.resi+10)
		cmd.align(fr,to)
		cmd.show('sticks',self.rdes.sel())
		cmd.color('white',self.rdes.sel()+" and elem C")
		cmd.center(self.rdes.sel())
		v = list(cmd.get_view())
		v[11] = -80.0
		cmd.set_view(tuple(v))
		cmd.clip('slab',50,self.rdes.sel())
	
	def showpack(self):
		print self.toggle['sph']
		#cmd.hide('spheres','not '+self.manager.d.bupsel)
		cmd.hide('spheres')
		if self.toggle['sph']:
			self.toggle['sph'] = False
			return
		fr = self.rnat.obj+" and name n+ca+c and resi %i-%i"%(self.rnat.resi-10,self.rnat.resi+10)		
		to = self.rdes.obj+" and name n+ca+c and resi %i-%i"%(self.rdes.resi-10,self.rdes.resi+10)
		cmd.align(fr,to)
		cmd.center(self.rdes.sel()+" or "+self.rnat.sel())
		#cmd.show('sticks',self.rnat.sel())
		cmd.show('spheres',"byres ((not "+self.rnat.obj+") within 6 of "+self.rdes.sel()+")")
		self.toggle['sph'] = True
	
	def showinfo(self):
		m = dict(self.scores)
		m['hbond'     ] = sum((self.res_scores[k] for k in self.res_scores if k.startswith("hbond_")))
		m['rep'       ] = self.res_scores      ['fa_rep'     ]		
		m['tot_sc'    ] = self.manager.d.scores['sc'         ]
		m['tot_rholes'] = self.manager.d.scores['packing'    ]
		m['tot_area'  ] = self.manager.d.scores['sc_int_area']
		self.showinfomap(m,"basic")
	
	def showallinfo(self):
		m = dict(self.scores)
		m.update(self.res_scores)
		m["tot_score"        ] = self.manager.d.scores["score"        ]
		m["tot_fa_atr"       ] = self.manager.d.scores["fa_atr"       ]
		m["tot_fa_rep"       ] = self.manager.d.scores["fa_rep"       ]
		m["tot_fa_sol"       ] = self.manager.d.scores["fa_sol"       ]
		m["tot_fa_intra_rep" ] = self.manager.d.scores["fa_intra_rep" ]
		m["tot_pro_close"    ] = self.manager.d.scores["pro_close"    ]
		m["tot_fa_pair"      ] = self.manager.d.scores["fa_pair"      ]
		m["tot_hbond_sr_bb"  ] = self.manager.d.scores["hbond_sr_bb"  ]
		m["tot_hbond_lr_bb"  ] = self.manager.d.scores["hbond_lr_bb"  ]
		m["tot_hbond_bb_sc"  ] = self.manager.d.scores["hbond_bb_sc"  ]
		m["tot_hbond_sc"     ] = self.manager.d.scores["hbond_sc"     ]
		m["tot_dslf_ss_dst"  ] = self.manager.d.scores["dslf_ss_dst"  ]
		m["tot_dslf_cs_ang"  ] = self.manager.d.scores["dslf_cs_ang"  ]
		m["tot_dslf_ss_dih"  ] = self.manager.d.scores["dslf_ss_dih"  ]
		m["tot_dslf_ca_dih"  ] = self.manager.d.scores["dslf_ca_dih"  ]
		m["tot_rama"         ] = self.manager.d.scores["rama"         ]
		m["tot_omega"        ] = self.manager.d.scores["omega"        ]
		m["tot_fa_dun"       ] = self.manager.d.scores["fa_dun"       ]
		m["tot_p_aa_pp"      ] = self.manager.d.scores["p_aa_pp"      ]
		m["tot_ref"          ] = self.manager.d.scores["ref"          ]
		m["tot_ddG"          ] = self.manager.d.scores["ddG"          ]
		m["tot_air_energy"   ] = self.manager.d.scores["air_energy"   ]
		m["tot_air_fa_atr"   ] = self.manager.d.scores["air_fa_atr"   ]
		m["tot_air_fa_rep"   ] = self.manager.d.scores["air_fa_rep"   ]
		m["tot_air_fa_dun"   ] = self.manager.d.scores["air_fa_dun"   ]
		m["tot_unsat_hbs"    ] = self.manager.d.scores["unsat_hbs"    ]
		m["tot_des_pos"      ] = self.manager.d.scores["des_pos"      ]
		m["tot_packing"      ] = self.manager.d.scores["packing"      ]
		m["tot_avg_deg"      ] = self.manager.d.scores["avg_deg"      ]
		m["tot_sasa_int_area"] = self.manager.d.scores["sasa_int_area"]
		m["tot_sc_int_area"  ] = self.manager.d.scores["sc_int_area"  ]
		m["tot_sc"           ] = self.manager.d.scores["sc"           ]
		self.showinfomap(m,"all")
	
	def showinfomap(self,scmap,infokind):
		cmd.delete("pseudo*")
		cmd.remove("labelpseudo*")
		self.toggle['info'] = self.toggle['info'] and self.toggle['infokind']==infokind #and cmd.select("pseudo*") 
		if not self.toggle['info']:
			cen = com(self.rdes.sel())
			v = cmd.get_view()
			sp = math.sqrt(-v[11]/200.0)
			delta = Vec(v[1],v[4],v[7])*sp
			cen = cen - 0.5*len(scmap)*sp*delta
			kys = scmap.keys()
			kys.sort()
			for i,s in enumerate(reversed(kys)):
				p = cen + delta*(float(i)*sp)
				txt = s+": %2.3f"%scmap[s]
				txt.strip("0")
				cmd.pseudoatom(name="labelpseudo",pos=[p.x,p.y,p.z],label=txt)		
		self.toggle['info'] = not self.toggle['info']
		self.toggle['infokind'] = infokind
	
	def aa(self): return aa_3_1[self.rdes.resn]	
	
	def addaa(self,aa):	self.aas.extend(aa)
	
	def revert(self):
		v = cmd.get_view()
		cmd.remove(self.rdes.obj+" and not chain A")
		m = cmd.get_model(self.rdes.obj)
		n = self.oldcrd
		if not n:
			cmd.create("tmp12345",self.rnat.sel())
			cmd.align("tmp12345",self.rdes.sel())
			n = cmd.get_model("tmp12345").atom
			cmd.delete("tmp12345")
		di,ni = 0,0
		while m.atom[di].resi != str(self.rdes.resi): di += 1
		dj = di
		while m.atom[dj].resi == str(self.rdes.resi): dj += 1
		if self.oldcrd: self.oldcrd = None
		else:           self.oldcrd = m.atom[di:dj]
		m.atom = m.atom[:di] + n + m.atom[dj:]
		for i in range(di,di+len(n)):
			m.atom[i].resi  = str(self.rdes.resi)
			m.atom[i].chain = str(self.rdes.chain)
		cmd.load_model(m,self.rdes.obj,1)
		cmd.save("tmp.pdb",self.rdes.obj)
		cmd.delete(self.rdes.obj)
		cmd.load("tmp.pdb",self.rdes.obj,1)
		cmd.show('car',self.rdes.obj)
		cmd.show('lines')		
		redopent(self.rdes.obj)
		cmd.set_view(v)
		print "revert removing "+"".join(self.aas)+" from aas!"
		self.aas = [getaa('A',self.rdes.resi,self.manager.d.obj)]
	
	def resfileline(self,sp=1):
		return "%4i A PIKAA "%self.rdes.resi + "".join(self.aas).ljust(sp) + " # "+" # ".join(self.notes)
	


class Design(object):
	def __init__(self,fname,manager):
		self.fname = fname
		self.manager = manager
		self.muts = None
		self.remembermm = None
		self.remembermv = None
		self.mutmap = {}
		self.toggle = {'bup':False}
		self.scores = {}
		self.obj = os.path.basename(fname)
		if self.obj.endswith('gz' ): self.obj = self.obj[:-3]
		if self.obj.endswith('pdb'): self.obj = self.obj[:-4]
		self.pid = self.obj[:4]
		self.resscoreslbl = os.popen("grep '^label fa_atr' "+self.fname).read().split()
		self.resscoreslbl[-1] = "score" # instead of total
	
	def getloadname(self):
		i = 0
		while os.path.exists(sd+self.obj+"_%i"%i): i+=1
		mansave = sd+self.obj+"_%i.pdb"%(i-1)
		if i==0: mansave = None
		if mansave:
			print "LOAD FROM SAVED",mansave
			return mansave
		elif os.path.exists("mutalyze_PYMOL_WORKING/"+self.obj+".pdb"):
			print "LOAD FROM AUTOSAVED","mutalyze_PYMOL_WORKING/"+self.obj+".pdb"
			return "mutalyze_PYMOL_WORKING/"+self.obj+".pdb"
		else:
			return self.fname
	
	def getsavename(self,sd):
		i = 0
		while os.path.exists(sd+self.obj+"_%i"%i): i+=1
		return sd+self.obj+"_%i"%i
	
	def load(self):
		self.manager.d = self
		self.manager.m = None
		self.manager.m = None # to also set prevm = None		
		# cmd.delete(self.manager.prevd.obj)
		# cmd.delete(self.manager.prevd.pid)
		# print self.manager.d.obj
		# print self.manager.prevd.obj		
		# sys.exit()
		cmd.delete("all")
		cmd.load(self.getloadname(),self.obj,1)
		cmd.remove(self.obj+" and not chain A")
		cmd.fetch(self.pid)
		cmd.remove("het or hydro")
		cmd.color('gray',self.pid+' and elem C')
		self.scores['rms'] = cmd.align(self.pid,self.obj+" and chain A")[0]
		cmd.hide('ev')
		cmd.show('lines')
		cmd.show("car")
		#cmd.disable(self.pid)
		redopent(self.obj)
		self.recalc_design_pos()
		self.read_data_dir('avg_deg')
		self.read_data_dir('ddG')
		self.read_data_dir('rot_boltz')
		self.read_bup()
		self.read_tot_scores()
		self.read_res_scores()
		cmd.orient(self.obj+" or pnt*")
		self.manager.m = self.remembermm
		if self.manager.m: self.manager.m.focus()
		if self.remembermv: cmd.set_view(self.remembermv)
		
	
	def save(self,manual=False,savepse=False):
		sd = "mutalyze_MANUAL_SAVE/" if manual else savedir
		if not os.path.exists(sd): os.mkdir(sd)
		redopent(self.obj)
		sname = self.getsavename(sd)
		if savepse: cmd.save(sname+".pse")
		cmd.save (sname+".pdb",self.obj+" or pnt*")
		with open(sname+".resfile","w") as o: o.write(self.resfile()+"\n")
		#print "saved pse, pdb, and resfile with notes"
		self.remembermm = self.manager.m
		self.remembermv = cmd.get_view()
	
	def resfile(self):
		self.resupdate()
		sp = max(len(dp.aas) for dp in self.muts)
		rf = "AUTO\nNATRO\n\nstart\n"
		for m in self.muts: rf += m.resfileline(sp)+'\n'
		return rf
	
	def resupdate(self):
		for m in self.muts:
			aa = getaa('A',m.rdes.resi,self.obj)
			aa3 = aa_1_3[aa]
			if aa3 != m.rdes.resn and aa3 != m.rnat.resn:
				m.aas = [aa]
	
	def get_design_pos(self):
		if not self.muts: recalc_design_pos()
		return self.muts
	
	def recalc_design_pos(self):
		res = [Residue(x[1],x[0],self.obj) for x in getres(self.obj+" and chain A")]
		re2 = []
		#first,last = 0,-1
		c0,r0,cN,rN = None,None,None,None
		while r0 is None:
			try: 
				c0,r0 = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[0].sel()))[0]
			except IndexError:
				#first += 1
				res = res[1:]
		while rN is None:
			try: 
				cN,rN = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[-1].sel()))[0]
			except IndexError:
				#last -= 1
				res = res[:-1]
		assert c0 is not None and r0 is not None and cN is not None and rN is not None 
		assert c0==cN # and first < 10 and last > -10
		c = c0
		try:
			for j,i in enumerate(range(r0,rN+1)):
				r2 = Residue(i,c,self.pid)
				re2.append(r2)
		except IndexError:
			raise Exception("BAD NATIVE")
		self.muts = [DesignPos(r,r2) for r,r2 in zip(res,re2) if r.isdiffenent(r2)]
		for m in self.muts: m.manager = self.manager
		self.mutmap = {}
		for m in self.muts: self.mutmap[m.rdes.resi] = m
		for m in self.mutmap.values(): assert m in self.muts			
	
	def read_data_dir(self,dname,mytype=float):
		tmpf = glob.glob("%s/*%s*.%s"%(dname,self.obj.replace(".mutalyze",""),dname))
		assert len(tmpf) == 1
		#print "reading",dname,'from',tmpf[0]
		fres,pres = set(),set()
		with open(tmpf[0]) as f:
			ref = None
			for line in f.readlines():
				lab,val = line.split()
				val = mytype(val)-ref if ref else mytype(val)
				if lab == dname: ref = mytype(val)
				if lab[-1] not in list("0123456789"): lab = lab[:-3]
				try: resi = int(lab[3:])
				except Exception: continue
				fres.add(resi)
				if resi in self.mutmap:
					self.mutmap[resi].scores[dname] = val
					if self.mutmap[resi].rdes.resn != lab[:3]:
						print "WARNING: resn's disagree: file: %s pymol: %s"%(lab[:3],self.mutmap[resi].rdes.resn)
		for k in self.mutmap.keys(): pres.add(k)
		for i in fres|pres: 
			#if i not in pres and i in fres: print dname+":","resi in file not in pymol:%4i"%i
			#if i in pres and i not in fres: print dname+":","resi in pymol not in file:%4i"%i
			#if i in pres and i in fres:     print dname+":","resi in pymol and in file:%4i GOOD!"%i
			for i in pres-fres: self.mutmap[i].scores[dname] = float('nan')
	
	def read_tot_scores(self):
		tmpf = glob.glob("output/*%s*.sc"%(self.obj.replace(".mutalyze","")))
		if len(tmpf)==0:
			print "WARNING: no mutalyze scorefile found:","output/*%s*.sc"%(self.obj.replace(".mutalyze",""))
			return
		assert len(tmpf) == 1
		#print "reading total scores from",tmpf[0]
		with open(tmpf[0]) as f:
			lbl = f.readline().split()
			val = f.readline().split()
			for l,v in zip(lbl[1:-1],val[1:-1]):
				try: self.scores[l] = float(v)
				except Exception as e: pass#print e
	
	def read_res_scores(self):
		lines = os.popen("egrep '^[A-Z]{3}(_p:[^ ]+)?_[0-9]+ -?[0-9]' %s"%self.fname).readlines()
		#print "reading residue scores from",self.fname
		for line in lines:
			line = line.split()
			resi = int(line[0].split("_")[-1])
			for k,v in zip(self.resscoreslbl[1:],line[1:]):					
				if resi in self.mutmap:
					try: 
						self.mutmap[resi].res_scores[k] = float(v)
					except Exception as e:
						pass#print e
		# for m in self.muts:
		# 	for k in self.resscoreslbl:
		# 		if k not in m.res_scores: 
		# 			m.res_scores[k] = float('nan')
	
	def read_bup(self):
		dname = "hbonds"
		tmpf = glob.glob("%s/*%s*.%s"%(dname,self.obj.replace(".mutalyze",""),dname))
		assert len(tmpf) == 1
		#print "reading",dname,'from',tmpf[0]
		self.bup = []
		with open(tmpf[0]) as f:
			for line in f.readlines():
				lab,val = line.split()
				resi = int(lab[3:])
				for aname in val.split(','):
					self.bup.append((resi,aname))
		self.bupsel = self.obj+" and ("+" or ".join("(resi %i and name %s)"%(r,a) for r,a in self.bup)+")"
	
	def showbup(self):
		if not self.toggle['bup']: cmd.show("spheres",self.bupsel)
		else:                        cmd.hide("spheres",self.bupsel)
		self.toggle['bup'] = not self.toggle['bup']
	


def mutalyze():
	manager = MutalyzeManager()


def __init__(self):
	self.menuBar.addmenuitem('Plugin','command','Mutalyze',label='m...ma...mutalyze!',command=mutalyze)
	cmd.extend('m',mutalyze)




