#!/usr/bin/env python
#
import pickle
import Numeric, cephes, RandomArray
import sys

class Test:

    """ 
    There are two reasons why we don't rely on test.regrtest:
    first, putting the expected results inside the test_
    script would lead to very small coverage, or VERY HUGE test_ 
    files; second, I liked the idea of trying to evenly cover the
    configuration space, avoiding deterministic lattices; third, I never
    pickled variables, and wanted to try!  """

    def __init__(self,fn,fnname,**args):
        self.name=fnname
        self.reffile='ref_'+self.name+'.pkl'
        self.call=fn
        self.vars_read=(0==1)
        if args.has_key('ref'):
            self.ref = args['ref']
        else:
            self.ref = (0 == 1)
        self.in_vars=args['in_vars']
        self.out_vars=args['out_vars']
        if args.has_key('tries'):
            self.tries=args['tries']
        else:
            self.tries=100

    def _readref(self):
        if not self.ref:
            f=open(self.reffile,'r')
            p=Numeric.Unpickler(f)
            for t in self.in_vars.keys():
                self.in_vars[t]=p.load()
            for t in self.out_vars.keys():
                self.out_vars[t]=p.load()
            f.close()
            self.vars_read=(0==0)

    def _genref(self):
        if self.ref:
            f=open(self.reffile,'w')
            p=Numeric.Pickler(f)
            for t in self.in_vars.keys():
                self.in_vars[t]=self._gen_array(self.in_vars[t])
                p.dump(self.in_vars[t])
            self._compute()
            if type (self.result) != type (()): self.result=self.result,
            for t in self.result:
                p.dump(t)
            f.close
            
    def _gen_array(self,limits):
        seed=RandomArray.seed
        random=RandomArray.uniform
        for t in limits:
            if type(t)==type(0.+0.j): _complex=(0==0)
            else: _complex=(0==1)
        if _complex:
            seed()
            minr=min(limits[0].real,limits[1].real)
            maxr=max(limits[0].real,limits[1].real)
            mini=min(limits[0].imag,limits[1].imag)
            maxi=max(limits[0].imag,limits[1].imag)
            a=random(minr,maxr,(self.tries,))+0.j
            a.imag=random(mini,maxi,(self.tries,))
        else:
            minr=min(limits[0],limits[1])
            maxr=max(limits[0],limits[1])
            a=random(minr,maxr,(self.tries,))
        return a

    def _compute(self):
        self.result=apply(self.call,tuple(self.in_vars.values()))

    def test(self):
        self.max_rel_dev=[]
        if self.ref:
            self._genref()
        else:
            if not self.vars_read:
                self._readref()
            self._compute()
            if type (self.result) != type(()): self.result=self.result,
            for t in range(len(self.out_vars.keys())):
                dev=abs(self.result[t]-self.out_vars[self.out_vars.keys()[t]])
                ref=abs(self.result[t]+self.out_vars[self.out_vars.keys()[t]])/2
                mx_dev_idx=Numeric.argmax(dev)
                if dev[mx_dev_idx] > 0.:
                    if ref[mx_dev_idx] > 0.:
                        self.max_rel_dev.append(dev[mx_dev_idx]/ref[mx_dev_idx])
                    else:
                        self.max_rel_dev.append(1.e+38)
        if len(self.max_rel_dev)>0:
            return max(self.max_rel_dev)
        return 0
