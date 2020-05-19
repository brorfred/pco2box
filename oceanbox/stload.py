
import os 
import numpy as np

def save(S,name):


    savedict ={}

    for at in S.__dict__:
        if type(S.__dict__[at]) == str:
            savedict[at] = S.__dict__[at]
        elif type(S.__dict__[at]) == list:
            savedict[at] = S.__dict__[at]
        elif type(S.__dict__[at]) == int:
            savedict[at] = S.__dict__[at]
        elif type(S.__dict__[at]) == float:
            savedict[at] = S.__dict__[at]
        elif type(S.__dict__[at]) == np.ndarray:
            savedict[at] = S.__dict__[at]
        elif type(S.__dict__[at]) == np.ma.core.MaskedArray:
            savedict[at] = S.__dict__[at].data
            savedict[at + '_mask'] = S.__dict__[at].mask
    
    np.savez(name,**savedict)

def load(name,S=-999, t1=None, t2=None, k1=None, k2=None,
         i1=None,i2=None,j1=None,j2=None):

    if S == -999:
        class S:pass
    tl = np.load(name)

    pref = os.path.split(name)[1][:-4]
    try:
        if not '__tm__' in tl.files:
            f = np.load(pref + "_ijk.npz")
            __tm__ = f['tm'].item()
            __km__ = f['km'].item()
            __im__ = f['im'].item()
            __jm__ = f['jm'].item()   
            
        if not '__shape__' in tl.files:
            __shape__ = np.load(pref + "_shape.npz")["shp"].item()
        load = sliceload
    except IOError:
        def load(var): return tl[var]

    def sliceload(var):
        if __shape__[var] == (__tm__,__km__,__im__,__jm__):
            return tl[f][t1:t2,k1:k2,i1:i2,j1:j2]
        elif __shape__[var] == (__tm__,__im__,__jm__):
            return tl[f][t1:t2,i1:i2,j1:j2]
        elif __shape__[var] == (__km__,__im__,__jm__):
            print(var + " kij")
            return tl[f][k1:k2,i1:i2,j1:j2]
        elif __shape__[var] == (__im__,__jm__):
            return tl[f][i1:i2,j1:j2]
        elif  len(__shape__[var]) == 1:
            if __shape__[var] == (__jm__):
                return tl[f][j1:j2]
            elif __shape__[var] == (__im__):
                return tl[f][i1:i2]
        return tl[f]
        
    for f in tl.files:
        if f + '_mask' in tl.files:
            field = load(f)
            mask =  load(f + '_mask')
            setattr(S, f, np.ma.masked_array(field,mask))
        if '_mask' in f:
            pass
        else:
           setattr(S, f, load(f))
    return S
