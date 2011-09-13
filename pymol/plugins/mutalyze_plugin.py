from pymol import cmd
import sys, urllib2, zlib, os
#from easygui import *
from Tkinter import *

sys.path.append("/Users/sheffler/lib/pymol")
from util import Mat,Vec,getres,aa_3_1,aa_1_3,redopent,rot,trans

def __init__(self):
   self.menuBar.addmenuitem('Plugin', 'command',
                            'Mutalyze',
                            label = 'm...ma...mutalyze!',
                            command = lambda s=self : mutalyze(s))

class Application(Frame):
	def say_hi(self):
	    print "hi there, everyone!"
	
	def reset_align(self):
		cmd.align(self.mut[0].rnat.obj,self.mut[0].rdes.obj)
	
	def createWidgets(self,mutations):
		self.mut = mutations
		self.bf = Frame(self)
		self.quit = Button(self.bf,text="DONE",command=self.quit).pack({"side": "top"})
		self.rsta = Button(self.bf,text="Reset Align",command=self.reset_align).pack({"side": "top"})
		self.zoom = Button(self.bf,text="Zoom",command=self.zoom).pack({"side": "top"})
		self.cena = Button(self.bf,text="Center",command=self.cenall).pack({"side": "top"})
		self.topri = Button(self.bf,text="To Primary",command=self.toprimary).pack({"side": "top"})
		self.charge = StringVar()
		self.charge.set("foo")
		self.crgu = Button(self.bf,textvariable=self.charge,command=self.chargeupdate).pack({"side": "top"})
		self.bf.pack({'side':'top'})
		self.mutframe = Frame(master=self)
		for m in mutations:
			f = Frame(self.mutframe)
			b = Button(f,width=15)
			b['text'] = m.rnat.resn+`m.rnat.resi`+" to "+m.rdes.resn+`m.rdes.resi`
			b['command'] = m.focus
			b2 = Button(f)
			b2['text'] = "R"
			b2['command'] = m.revert
			b3 = Button(f)
			b3['text'] = "S"
			b3['command'] = m.focussph
			e = Entry(f,width=8)
			e.insert(0,"%4d A PIKAA %s # "%(m.rdes.resi,aa_3_1[m.rdes.resn]))
			e.pack({'side':'right'})
			b3.pack({'side':'right'})
			b.pack({'side':'right'})
			b2.pack({'side':'right'})
			f.pack({'side':'top'})
		self.mutframe.pack({"side":"left"})
		
	
	def __init__(self,des,master=None):
		Frame.__init__(self, master)		
		self.pack()
		self.des = des
		self.mut = des.get_mutations()
		self.createWidgets(self.mut)
	
	def zoom(self):
		cmd.zoom(self.mut[0].rdes.obj+" or pnt*")
	
	def cenall(self):
		cmd.center(self.mut[0].rdes.obj+" or pnt*")
	
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
	
	def focus(self):
		cmd.hide('sticks')
		fr = self.rnat.obj+" and name n+ca+c and resi %i-%i"%(self.rnat.resi-10,self.rnat.resi+10)		
		to = self.rdes.obj+" and name n+ca+c and resi %i-%i"%(self.rdes.resi-10,self.rdes.resi+10)
		cmd.align(fr,to)
		cmd.center(self.rdes.sel()+" or "+self.rnat.sel())
		cmd.show('sticks',self.rdes.sel())
		#cmd.show('sticks',self.rnat.sel())
	
	def focussph(self):
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
		redopent(self.rdes.obj)
		cmd.set_view(v)
		

class Design(object):
	def __init__(self,obj):
		v = cmd.get_view()
		self.obj = obj
		cmd.remove(self.obj+" and not chain A")
		self.pid = obj[:4]
		cmd.fetch(self.pid)
		cmd.remove("het or hydro")
		cmd.color('gray',self.pid+' and elem C')
		assert cmd.align(self.pid,self.obj+" and chain A")[0] < 0.6
		cmd.hide('ev')
		cmd.show('lines')
		cmd.show("car")
		redopent(self.obj)
		cmd.set_view(v)
	
	def get_mutations(self):
		res = [Residue(x[1],x[0],self.obj) for x in getres(self.obj+" and chain A")]
		re2 = []
		c,r0 = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[ 0].sel()))[0]
		c,rN = getres('%s within 1.0 of (%s and name ca)'%(self.pid,res[-1].sel()))[0]
		for i in range(r0,rN+1):
			r2 = Residue(i,c,self.pid)
			re2.append(r2)
		return [Mutation(r,r2) for r,r2 in zip(res,re2) if r.isdiffenent(r2)]
	

def mutalyze(app):
	obj = [x for x in cmd.get_object_list() if x[4:9]==".pdb."]
	assert len(obj) == 1
	obj = obj[0]
	
	root = Tk()
	app = Application(Design(obj),root)
	app.mainloop()
	root.destroy()
	
    
