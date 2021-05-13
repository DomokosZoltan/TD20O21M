# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 19:53:13 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep
import numpy as np


class SignalIterator:
    
    
    def __init__(self, signal, domain, order = 'F', which = 'bs', mode = 'ms'):
        
        self.signal = signal
        
        self.domain = domain
        
        self.__it_dm = domain.__iter__(order, which)
        
        self.mode = mode
        
        return
    
    
    def __next__(self):
        
        res = None
        
        try:
        
            inds = self.__it_dm.nd_index()
            
            self.__it_dm.__next__()
            
            sigval = self.signal.fcval_at(self.domain, inds)
            
            meas = self.domain.measure_at(inds)
            
            if 'ms' == self.mode:
                
                res = (meas, sigval)
                
            else:
                
                res = sigval
        
        except StopIteration:
        
            raise StopIteration
        
        return res
    

#TODO: add flat iterator and ND iterator to this class
class Signal:

    def __init__(self, labels, coeffs, bsgens):
        
        self.labels = deep(labels)
        
        #coefficients of base functions
        self.coeffs = deep(coeffs) #Can be simple array, as product/replace etc.
        #is not supported for Signals. 
        
        self.bsgens = deep(bsgens) #Generator object that maps a base
        #function to a given multidim. label.

        return


    #iterate the signal w. respect to the inputted domain
    def __iter__(self, dm, order = 'F', which = 'bs', mode = 'ms'):
        
        res = None
        
        res = SignalIterator(self, dm, 'F', 'bs', 'ms')
        
        return res


    def fcval_at(self, dm_src, dm_dst, inds):
        
        res = None

        with np.nditer(self.coeffs, flags = ['multi_index'], 
                       op_flags=['readonly'], order='F') as it:
                    
            for coeff in it:
                    
                bsinds = list(it.multi_index)
                    
                bsfunc = self.bsgens.gen(dm_src, bsinds)
                bsfcval = bsfunc.fcval_at(dm_dst, inds)
                    
                if res is None:
                        
                    res = coeff * bsfcval
                        
                else:
                        
                    res += coeff * bsfcval
        
        return res
    
    
    def scmult(self, sc):
        
        res = None
        
        labels = deep(self.labels)
        
        coeffs = deep(self.coeffs)
        
        bsgens = deep(self.bsgens)
        
        for bsinds in labels:
            
            coeffs[bsinds] *= sc
        
        res = Signal(labels, coeffs, bsgens)
        
        return res 
    
    
    def add(self, rhs):
        
        res = None

        labels = deep(self.labels)
        
        coeffs = deep(self.coeffs)
        
        bsgens = deep(self.bsgens)
        
        for bsinds in labels:
            
            coeffs[bsinds] += rhs.coeffs[bsinds]
        
        res = Signal(labels, coeffs, bsgens)
        
        return res


    def tenzor_prod(self, rhs):
        
        res = None
        
        #TODO
        
        return res 
    

    def scalef(self):
        
        #TODO
        
        pass

#HELPERS
def sigvec(sigfc = None, src_dm = None, dst_dm = None, lod = None, lims = None, seldims = None):
        
    dmsp = dst_dm.ns
        
    dmls = []
            
    for obj in dmsp:
            
        dmls.append( obj[1] )
        
    #Create the numpy array that fits dims in limits
    dmvec = np.zeros(shape = dmls, dtype = np.float32)
    codmvec = np.zeros(shape = dmls, dtype = np.complex64)
        
    dmls = []
       
    #Fill the array with domain point coordinates
    with np.nditer(codmvec, flags = ['multi_index'], 
                       op_flags=['writeonly'], order='F') as it:
            
        for elm in it:
                
            inds = list(it.multi_index)
                
            dmls.append(dst_dm.sample_at(inds) )
                
            elm[...] = sigfc(src_dm, dst_dm, inds)
            
    dmvec = np.asarray(dmls)
        
    res = [dmvec, codmvec] 
        
    #Return result
    return res


#TODO: move this function to codomain.py
def cdmvec(sigfc = None, domain = None, lod = None, lims = None, seldims = None):
        
    dmsp = domain.ns
        
    dmls = []
            
    for obj in dmsp:
            
        dmls.append( obj[1] )
        
    #Create the numpy array that fits dims in limits
    dmvec = np.zeros(shape = dmls, dtype = np.float32)
    codmvec = np.zeros(shape = dmls, dtype = np.complex64)
        
    dmls = []
       
    #Fill the array with domain point coordinates
    with np.nditer(codmvec, flags = ['multi_index'], 
                       op_flags=['writeonly'], order='F') as it:
            
        for elm in it:
                
            inds = list(it.multi_index)
                
            dmls.append(domain.sample_at(inds) )
                
            elm[...] = sigfc(domain, inds)
            
    dmvec = np.asarray(dmls)
        
    res = [dmvec, codmvec] 
        
    #Return result
    return res


#TODO: review this workaround and find a final solution, consider
#using a 'dict'.
def scalev(sigvec, domain, scmd):
    
    dmsp = []
    
    for obj in domain.ns:
        
        dmsp.append( obj[1] )
    
    res = np.zeros(shape = dmsp, dtype = np.complex64)

    if scmd == 'sqrt(dx)':
            
        sc = lambda x: np.sqrt(domain.measure_at(list(x) ) )
            
    elif  scmd == 'sqrt(1/dx)':
            
        sc = lambda x: 1.0 / np.sqrt(domain.measure_at(list(x) ) )
        
    elif scmd == 'sqrt(1/n)' :
        
        full_n = 1
        for n in domain.ns:
            full_n *= n
        
        sc = lambda x: 1.0 / np.sqrt(full_n)
    
    elif scmd == 'sqrt(avg_dx)':
    
        avg_dx = 1.0
        for i in range(0, domain.ns.last_did() + 1):
            avg_dx *= (  (domain.bandlims.get_obj((i,i))[1]\
                          - domain.bandlims.get_obj((i,i))[0]) / domain.n(i) )
        
        sc = lambda x: np.sqrt(avg_dx)
        
    elif scmd == 'sqrt(1/avg_dx)':
    
        avg_dx = 1.0
        for i in range(0, domain.ns.last_did() + 1):
            avg_dx *= ( (domain.bandlims.get_obj((i,i))[1]\
                         - domain.bandlims.get_obj((i,i))[0]) / domain.n(i) )
        
        sc = lambda x: 1.0 / np.sqrt(avg_dx)
    
    elif scmd == 'sqrt(avg_df)':
    
        sc = lambda x: 0.5 #TODO: use bounds and cardinalities
        
    elif scmd == 'sqrt(1/avg_df)':
    
        sc = lambda x: 1.0/0.5 #TODO: use bounds and cardinalities
    
    else:
        #TODO: raise exception
        pass

    
    with np.nditer(res, flags = ['multi_index'], op_flags=['writeonly'], order='C') as it:
            
        for elm in it:
            
            elm[...] = sigvec[it.multi_index] * sc(it.multi_index)
    
    return res