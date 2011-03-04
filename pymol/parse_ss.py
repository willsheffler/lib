import pymol

def parsess(model_name,label='',size_thresh=3,base_sel='',sstypes=('H','S')):
    if label=='':
        label = model_name
    if base_sel=='':
        base_sel = "name CA or name O or name N or name C"
    m = cmd.get_model(model_name)
    ss = ['']*1000
    for a in m.atom:
        if a.ss!='':
            ss[int(a.resi)] = a.ss
    while ss[-1]=='':
        ss.pop()
    p = {}
    for x in ss:
        p[x] = []
    p[ss[0]].append(0)
    for ii in range(1,len(ss)):
        if ss[ii-1]!=ss[ii]:
            p[ss[ii-1]].append(ii-1)
            p[ss[ii]].append(ii)
    p[ss[len(ss)-1]].append(len(ss)-1)
    selstr = []
    count = {'H':0,'L':0,'S':0}
    for k in sstypes:
        R = [(p[k][ii],p[k][ii+1]) for ii in range(len(p[k])) if not ii%2]
        for ii in range(len(R)):
            s = R[ii][0]
            e = R[ii][1]
            if e-s+1 >= size_thresh:
                count[k] += 1
                selstr.append((label+k.lower()+str(count[k]),
                        "( ( resi %i-%i ) and ( %s ) and ( %s ) )"%(s,e,model_name,base_sel)))
    for name,sel in selstr:
        cmd.select(name,sel)
    return selstr

def test():
    print 'test'

cmd.extend("parsess",parsess)

