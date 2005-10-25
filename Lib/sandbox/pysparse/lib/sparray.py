# -*- coding: iso-8859-1 -*-
from types import IntType,SliceType
import operator
import spmatrix
import Numeric

class sparray:
    """
    d-dimensionnal sparse array emulation by a long sparse vector.
    supports syntax like:
    
    a[2,1,6,5]=2
    a[:,5,6,6]=list of appropriate length (one slice at a time only)
    b=a[:,5,6,6] (b is a numeric array)
    b=a[n], a[n]=6  if n in range
    
    a=sparray((2,6,9,8),dicto,shifts)
        where, optionnally, dicto is a dictionnary whose keys are tuple
        in range of (2,6,9,8), and shifts is a tuple
        to shift origins in case the smallest coordinate in dicto
        is not (0,...,0)
    
    """
    def __init__(self,dim,dicto=None,shifts=None):
        """
        attributes: shifts, dims, data, is1D, length
        methods : dump
        """

        self.shifts=shifts
        
        if type(dim)==type(()) or type(dim)==type([]):
            self.data = spmatrix.ll_mat(reduce(operator.mul, dim), 1)
            self.dims = dim        
            if dicto:
                for k, v in dicto.iteritems():
                    shk = map(operator.__sub__, k, shifts)
                    self.data[self.comp(shk), 0]=v
        
        elif type(dim)==IntType:
            self.data = spmatrix.ll_mat(dim,1)
            self.dims = dim        
            if dicto:
                for k, v in dicto.iteritems():
                    shk = k - shifts
                    self.data[shk,0] = v
        
        self.is1D = type(self.dims)==IntType
        
        

    def __get_shape0(self):return self.data.shape[0]
    length = property(__get_shape0, doc="sparray linear length")
    
    def decomp(self,ind):
        "from linear to multi indice"
        a = ind
        l = len(self.dims)
        res = [0]*l
        for i in range(l - 1, -1, -1):
            a, b = divmod(a,self.dims[i])
            res[i] = b
        return tuple(res)
    
    def comp(self,indice):
        "from multi indice to linear"
        l = len(self.dims)
        a = 0
        for i in range(l-1):
            a += reduce(operator.mul, self.dims[i+1:]) * indice[i]
        a += indice[l-1]
        return a
        
    def __setitem__(self, k, value):
        if type(k) is IntType:
            self.data[k, 0] = value
            return
        
        vec = map(lambda x: type(x) is SliceType, k)
        
        if True in vec: # suppose only one slice
            ii = vec.index(True)
            indices = []
            k = list(k)
            comp = self.comp
            for i in range(self.dims[ii]):
                k[ii] = i
                self.data[comp(k), 0] = value[i]
        else:
            self.data[self.comp(k),0]=value
    
    
    def __getitem__(self,k):
        """
        output a Numeric vector if slice in coordinates
        (one slice present only)
        """
        if type(k) is IntType: return self.data[k, 0]
        
        vec = map(lambda x: type(x) is SliceType, k)
        
        if True in vec: #suppose only one slice
            ii=vec.index(True)
            indices=[]
            k = list(k)
            rep = Numeric.zeros((self.dims[ii],), 'd')
            for i in range(self.dims[ii]):
                k[ii] = i
                rep[i] = self.data[self.comp(k), 0]
            return rep
        else:
            return self.data[self.comp(k), 0]
            
    def dump(self):
        print self.data


