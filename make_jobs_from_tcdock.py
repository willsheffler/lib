#!/usr/bin/env python

import numpy as np
import pandas as pd
import glob
import re
import StringIO
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option('-o',dest='wd', default='.', help='dump setup scripts here')
op,args = parser.parse_args()



tplt = """./rosetta_scripts.static.linuxiccrelease @input/hyak_design.flags \
-s                input/scaffolds/%(pdb1)s_%(pdb2)s.pdb.gz \
-in::file::native input/scaffolds/%(pdb1)s_%(pdb2)s.pdb.gz \
-parser:script_vars a1a2=%(a1)f,%(a2)f r1r2=%(r1)f,%(r2)f flip_axis=%(flip)s,0 path=%(outdir)s sym=I32 symdof1=JCT00 symdof2=JCD00 \
-jd2:checkpoint_file %(outdir)s/%(tag)s.ckpt \
-out:path:all %(outdir)s \
-mute all \
> %(outdir)s/%(tag)s.log\n"""

def readtcdock(fn,**kwargs):
    with open(fn) as f:  buf = f.read()
    buf = re.sub(r"^.*\| ","",buf,flags=re.M)
    colname="""  tag                                    score   sc/nc sc/nwtd     rsc   nr   nc    nwtd    rsc0  nr0  nc0   nwtd0    cb1  nr1  nc1   nwtd1    cb2  nr2  nc2   nwtd2   diam
    cover   tdis   n1  a1       r1   n2  a2       r2 ori   mom1    mom2    mom3    mom4  v0.4  v0.8  v1.2  v1.6  v2.0  v2.4  v2.8  v3.2  v3.6  v4.0    CBScore      MH_RAW  MH_SPREAD  MH_SSPAIR
        M_EE       M_EL       M_HE       M_HH       M_HL       M_LL  M_SS_PAIR      M_num     M_xpsc  M_xpsc_ar  M_xpsc_av  M_xpsc_rt   no""".replace("/","_").replace(".","").split()
    t = pd.read_table(StringIO.StringIO(buf),delimiter="\s+",names=colname,verbose=True,**kwargs)
    t['pdb1'] = ["_".join(x.split("_")[0:3]) for x in t.tag]
    t['pdb2'] = ["_".join(x.split("_")[3:6]) for x in t.tag]
    t['arch'] = ["_".join(x.split("_")[6:-1]) for x in t.tag]
    t["forward"]   = [x[-1]=='F' for x in t.arch]
    t['arch'] = [x[:-1] for x in t.arch]
    t["hitnum"] = [int(x.split("_")[-1]) for x in t.tag]
    return t



seenit = set()
try: os.makedirs(op.wd)
except OSError: pass
with open(op.wd+"/setup.sh",'w') as setup:
    setup.write("""shuffle input/alljobs.sh > input/randjobs.sh
split -a3 -d -l 80 input/randjobs.sh input/genjobs/split_jobs_
for i in $(/bin/ls input/genjobs/split_jobs_*); do
    sed -e "s=JOB_FILE_REPLACEME=`basename $i`=g" input/template.pbs > input/genjobs/submit_`basename $i`;
done

""")
    with open(op.wd+"/alljobs.sh",'w') as jobs:
        for tcfile in args:
            df = readtcdock(tcfile)
            for index,row in df.iterrows():
            	a1,a2,r1,r2 = row['a1'],row['a2'],row['r1'],row['r2']
            	flip = '0' if row['forward'] else 'y'
            	outdir = "output/%(tag)s"%row
                tag = row['tag']
                pdb1 = row['pdb1']
                pdb2 = row['pdb2']
                p12 = pdb1+pdb2
                if p12 not in seenit:
                    seenit.add(p12)
                    setup.write("awk '{print substr($_,1,21)\"A\"substr($_,23,1000)}' /ssd8/pdb/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C3/%(pdb1)s.pdb | grep ATOM  > input/scaffolds/%(pdb1)s_%(pdb2)s.pdb\n"%row )
                    setup.write("awk '{print substr($_,1,21)\"B\"substr($_,23,1000)}' /ssd8/pdb/PISA_PDB_SYMMETRIC_SCAFFOLDS_BALE/C2/%(pdb2)s.pdb | grep ATOM >> input/scaffolds/%(pdb1)s_%(pdb2)s.pdb\n"%row )
                    setup.write("gzip input/scaffolds/%(pdb1)s_%(pdb2)s.pdb\n"%row)
                setup.write("mkdir -p %(outdir)s\n"%vars())
                jobs.write(tplt%vars())




# for i in $(/bin/ls input/genjobs/split_jobs_*); do sed -e "s=JOB_FILE_REPLACEME=`basename $i`=g" input/template.pbs > input/genjobs/submit_`basename $i`; done


