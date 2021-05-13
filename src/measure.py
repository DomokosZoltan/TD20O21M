# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:55:11 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep
import numpy as np


import mstruct as mst


def cube(volume):
    
    res = 0.0
    
    if not volume is None:
        res = 1.0
        
    for ind in range(0, len(volume) ):
        res *= ( volume[ind][1][0] - volume[ind][0][0] ) 

    return res


def sector(volume):
    
    res = 1.0
    
    #edge 'a'
    asx = volume[0][0][0]
    asy = volume[0][0][1]
    
    aex = volume[0][1][0]
    aey = volume[0][1][1]
    
    #edge 'b'
    bsx = volume[1][0][0]
    bsy = volume[1][0][1]
    
    bex = volume[1][1][0]
    bey = volume[1][1][1]
    
    #length of 'a' = delta_r
    delta_r = np.sqrt( (aex - asx) ** 2.0 + (aey - asy) ** 2.0)
    
    #start of 'a' is r0
    r0 = np.sqrt(asx ** 2.0 + asy ** 2.0)
    
    #r = r0 + delta_r
    r = r0 + delta_r
    
    #length of 'b'
    hypot = np.sqrt( (bex - bsx) ** 2.0 + (bey - bsy) ** 2.0)
    
    #dot
    dot = (aex - asx) * (bex - bsx) + (aey - asy) * (bey - bsy)
    
    if 0.0 < delta_r:
    
        beta = np.arccos( dot / (delta_r * hypot) )
        
        delta_phi = 2.0 * beta - np.pi 
        
        res = ( (r ** 2.0) - (r0 ** 2.0) ) * delta_phi / 2.0
            
    else:
        
        res = 0.0
    
    return res
    

class Measure():
    
    def __init__(self, measstr = None):
        
        if measstr is None:
            
            measfcs = [cube]
        
            fcsdims = [ (0,0) ]
            
        else:
            
            measfcs = deep(measstr._MdimStruct__objs)
        
            fcsdims = deep(measstr._MdimStruct__dims)
        
        self.__measstr = mst.MdimStruct(measfcs, fcsdims)
        
        return
    
    
    def measure_at(self, volume, bounds, inds):
        
        res = 1.0
        
        vols = volume.volume_at(bounds, inds)
        
        for obj in self.__measstr:
            
            dims = obj[0]
            
            fc = obj[1]
    
            sl = mst.dims_to_sl(dims)
            
            res *= fc(vols[sl] )
            
        return res            


    def prod(self, rhs):
        
        res = None
        
        res = Measure( self.__measstr.prod(rhs.__measstr) )
        
        return res
    
    
    def replace(self, rhs, dims, src_offs = 0):
        
        self.__measstr.replace(rhs.__measstr, dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
        
        res = Measure(self.__measstr.restrict(dims) )
        
        return res
    
