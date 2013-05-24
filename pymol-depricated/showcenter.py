from pymol import cmd

def showcenters(res,crad='2',rad='5'):
    cmd.hide('all')
    sel = "centers within %(crad)s of %(res)s/sd05"%vars()
    cmd.show('spheres',sel)
    cmd.color('purple', sel)
    cmd.show('spheres','(not centers within %(rad)s of (%(sel)s)) or (%(sel)s)'%vars())
    cmd.show('lines','(name CA,N,C,O)')
    cmd.center("%(res)s/sd05"%vars())

cmd.extend('showcenters',showcenters)
