import sys,urllib2,zlib,os,threading
from pymol import cmd
from Tkinter import *
sys.path.append("/Users/sheffler/lib/pymol")
from util import Mat,Vec,getres,aa_3_1,aa_1_3,redopent,rot,trans

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
		for m in mutations:
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
	def __init__(self,master=None):
		self.dir = "./output/"
		files = [f for f in os.listdir(self.dir) if f.endswith('.pdb') or f.endswith('.pdb.gz')]
		if not files: raise IOException("ast")
		self.des = [Design(self.dir+f,self) for f in files]
		self.d = self.des[0]
		self.m = None
		self.d.load()
		cmd.extend("n",lambda: self.next())
		cmd.extend("p",lambda: self.prev())
		cmd.extend("s",lambda: self.showpack())
		cmd.extend("i",lambda: self.info())
		cmd.extend("r",lambda: self.revert())
		cmd.extend("c",lambda: self.chargeupdate())
		cmd.extend("g"  ,lambda: self.gui())
		cmd.extend("gui",lambda: self.gui())		
		cmd.extend("addaa",lambda a: self.addaa(a))
		cmd.extend("a"    ,lambda a: self.addaa(a))		
		cmd.extend("note",lambda txt=None: self.note(txt))
		cmd.extend("resfile",lambda: self.resfile())
		print 'MutalyzeManager initialized managing %i designs'%len(self.des)
	
	def gui(self):
		root = Tk()
		app = MutalyzeGui(self,root)
		app.mainloop()
		root.destroy()    	
	
	def resfile(self):
		sp = max(len(m.aas) for m in self.d.muts)
		print "####################### BEGIN RESFILE #######################"
		print "AUTO\nNATRO\n\nstart"
		for m in self.d.muts: print m.resfileline(sp)
		print "######################## END RESFILE ########################"
	
	def reset_align(self):
		cmd.align(self.d.muts[0].rnat.obj,self.d.muts[0].rdes.obj)
		
	def info(self):
		hide('label','all')
	
	def addaa(self,aa):
		self.m.addaa(aa)
	
	def note(self,txt=None):
		if txt is None: print "\n".join(self.m.notes)
		else: self.m.notes.append(txt)
	
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
		pos = cmd.select("tmpcharge",self.des.obj+" and resn ARG+LYS and name CA")
		neg = cmd.select("tmpcharge",self.des.obj+" and resn ASP+GLU and name CA")
		cmd.delete("tmpcharge")
		c = "Total Charge: %i"%(pos-neg)
		print c
		self.charge.set(c)
	
	def next(self):
		try: self.m = self.d.muts[self.d.muts.index(self.m)+1]
		except ValueError: self.m = self.d.muts[0]			
		except IndexError as e:
			print "DONE!!!"
			return
		self.m.focus()
	
	def revert(self):
		self.m.revert()
	
	def showpack(self):
		self.m.showpack()
	

class Residue(object):
	def __init__(self,resi,chain,obj):
		self.resi = resi
		self.chain = chain
		self.obj = obj
		self.resn = cmd.get_model(self.sel()).atom[0].resn
	
	def sel(self):
		return "(resi %s and chain %s and obj %s)"%(self.resi,self.chain,self.obj)
	
	def __str__(self):
		return "%s/%s/%s/%i"%(self.obj,self.chain,self.resn,self.resi)
	
	def coords(self):
		return (Vec(a.coord) for a in cmd.get_model(self.sel()).atom)
	
	def isdiffenent(a,b):
		print a,b
		return a.resn != b.resn
	

class Mutation(object):
	def __init__(self,rdes,rnat):
		self.rdes = rdes
		self.rnat = rnat
		self.sph = False
		self.oldcrd = None
		self.notes = []
		self.aas = [self.aa()]
	
	def focus(self):
		self.manager.m = self
		cmd.hide('sticks')
		fr = self.rnat.obj+" and name n+ca+c and resi %i-%i"%(self.rnat.resi-10,self.rnat.resi+10)		
		to = self.rdes.obj+" and name n+ca+c and resi %i-%i"%(self.rdes.resi-10,self.rdes.resi+10)
		cmd.align(fr,to)
		cmd.center(self.rdes.sel()+" or "+self.rnat.sel())
		cmd.show('sticks',self.rdes.sel())
		#cmd.show('sticks',self.rnat.sel())
	
	def showpack(self):
		cmd.hide('spheres')
		if self.sph:
			self.sph = False
			return
		fr = self.rnat.obj+" and name n+ca+c and resi %i-%i"%(self.rnat.resi-10,self.rnat.resi+10)		
		to = self.rdes.obj+" and name n+ca+c and resi %i-%i"%(self.rdes.resi-10,self.rdes.resi+10)
		cmd.align(fr,to)
		cmd.center(self.rdes.sel()+" or "+self.rnat.sel())
		#cmd.show('sticks',self.rnat.sel())
		cmd.show('spheres',"byres ((not "+self.rnat.obj+") within 6 of "+self.rdes.sel()+")")
		self.sph = True
	
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
	
	def resfileline(self,sp=1):
		return "%4i A PIKAA "%self.rdes.resi + "".join(self.aas).ljust(sp) + " # "+" # ".join(self.notes)
	

class Design(object):
	def __init__(self,fname,manager):
		self.fname = fname
		self.manager = manager
		self.muts = None
		self.obj = os.path.basename(fname)
		if self.obj.endswith('gz' ): self.obj = self.obj[:-3]
		if self.obj.endswith('pdb'): self.obj = self.obj[:-4]
	
	def load(self):
		cmd.delete(self.obj)
		cmd.load(self.fname)
		cmd.remove(self.obj+" and not chain A")
		self.pid = self.obj[:4]
		cmd.fetch(self.pid)
		cmd.remove("het or hydro")
		cmd.color('gray',self.pid+' and elem C')
		assert cmd.align(self.pid,self.obj+" and chain A")[0] < 0.6
		cmd.hide('ev')
		cmd.show('lines')
		cmd.show("car")
		redopent(self.obj)
		cmd.zoom(self.obj)
		self.recalc_mutations()
	
	def get_mutations(self):
		if not self.muts: recalc_mutations()
		return self.muts
	
	def recalc_mutations(self):
		res = [Residue(x[1],x[0],self.obj) for x in getres(self.obj+" and chain A")]
		re2 = []
		c,r0 = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[ 0].sel()))[0]
		c,rN = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[-1].sel()))[0]
		for i in range(r0,rN+1):
			r2 = Residue(i,c,self.pid)
			re2.append(r2)
		self.muts = [Mutation(r,r2) for r,r2 in zip(res,re2) if r.isdiffenent(r2)]
		for m in self.muts: m.manager = self.manager
	

def mutalyze():
	manager = MutalyzeManager()

def __init__(self):
	self.menuBar.addmenuitem('Plugin','command','Mutalyze',label='m...ma...mutalyze!',command=mutalyze)
	cmd.extend('m',mutalyze)

cmd.extend('m',mutalyze)

