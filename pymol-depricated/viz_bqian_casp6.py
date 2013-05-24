import pymol

def viz_bqian_casp6(prot):
    cmd.delete('all')
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/nat_3way.pdb'  ,'n3way_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/good_3way.pdb' ,'g3way_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/bad_3way.pdb'  ,'b3way_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/nat_bsasa.pdb' ,'nbsasa_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/good_bsasa.pdb','gbsasa_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/bad_bsasa.pdb' ,'bbsasa_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/good_bsasa.pdb','grms_'+prot)
    cmd.load('~/project/modelpicking/bqian_casp6/'+prot+'/out/bad_bsasa.pdb' ,'brms_'+prot)
    cmd.hide('lines')
    for ii in "b3way g3way nbsasa gbsasa bbsasa grms brms".split():
        cmd.align(ii,'n3')
    cmd.show('spheres')
    cmd.center('visible')
    cmd.do('resrms %s, %s'%('grms','n3'))
    cmd.do('resrms %s, %s'%('brms','n3'))


def viz_bqian_casp7(prot):
    cmd.delete('all')
    md = "/users/bqian/submit/"+prot+'/submit/'
    cpdbd = "/users/bqian/submit/"+prot+'/analysis/colored_pdbs/'
    native = ""
    model  = ""
    robetta = ""
    for f in os.listdir(md):
        if f.endswith('.pdb'):
            print 'asd',f
            if f.startswith('model') and f.endswith('1.pdb'):
                model = f
                print 'asdf',f,model
            elif len(f) == 8:
                native = f
            elif f.startswith('robetta'):
                robetta = f
    print model, native, robetta
    cmd.load(cpdbd+native[:-4]+"_bsasa"+'.pdb','native')
    cmd.load(cpdbd+model[:-4]+"_bsasa"+'.pdb','mb_model_bsasa')
    cmd.load(cpdbd+"robetta_bsasa"+'.pdb',    'rb_robetta_3way')
    cmd.load(cpdbd+model[:-4]+"_3way"+'.pdb','m3_model_3way')
    cmd.load(cpdbd+"robetta_3way"+'.pdb',     'r3_robetta_bsasa')
    cmd.load(cpdbd+model[:-4]+"_bsasa"+'.pdb','mr_model_rms')
    cmd.load(cpdbd+"robetta_bsasa"+'.pdb',    'rr_robetta_rms')        

    for ii in "mb rb mr rr".split():
        cmd.align(ii,'native')
    cmd.hide('all')
    cmd.show('cartoon')
    cmd.do('select mb')
    cmd.do('center selected')
    cmd.do('delete sele')
    cmd.do('resrms %s, %s'%('mr','native'))
    cmd.do('resrms %s, %s'%('rr','native'))
    cmd.color('white','native')

def v(sel,rep='spheres'):
    cmd.hide('all')
    cmd.show(rep,sel)
    cmd.center('visible')

cmd.extend('viz',viz_bqian_casp7)
cmd.extend('v',v)

