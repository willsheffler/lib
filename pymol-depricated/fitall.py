import pymol,re


def sup( reference, whichobjs='.*' ):
    sel_template = '^'+whichobjs+'$'
    r = ''
    #selections = cmd.get_names('selections')
    #for s in selections:
    #    if re.match(sel_template,s) is not None:
    #        cmd.fit()
    models = cmd.get_names()
    for m in models:
        print m
        if re.match(whichobjs,m) is not None:
            cmd.fit(m,reference)

cmd.extend('sup',sup)
