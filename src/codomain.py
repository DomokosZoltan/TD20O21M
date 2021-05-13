# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 19:53:15 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep

import mstruct as mst
#import siglib as sgl

 
class Codomain:

    def __init__(self, codmstr):
       
        self.codmstr = deep(codmstr)
            
        return


    def prod(self, rhs):
        
        res = None
        
        codmstr = self.codmstr.prod(rhs.codmstr)
        
        res = Codomain(codmstr)
        
        return res
    

    def replace(self, rhs, dims, src_offs):
        
        self.codmstr.replace(rhs.codmstr, dims, src_offs)
        
        return
    

    def restrict(self, dims):
        
        res = None
        
        codmstr = self.restrict(dims)
        
        res = Codomain(codmstr)
        
        return res


    def fcval_at(self, domain, inds):
        
        sample_pt = domain.sample_at(inds)
        
        res = 1.0 + 0.0j
        
        for obj in self.codmstr:
            
            dims = obj[0]
            
            fc = obj[1]
            
            sl = mst.dims_to_sl(dims)
            
            res *= fc(sample_pt[sl] )
        
        return res

