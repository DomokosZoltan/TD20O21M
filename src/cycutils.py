# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 22:27:56 2021

@author: Domokos Zoltán
"""

import numpy as np


def get_lim_lod(lod, did, base = 32, fix_cyc = None):
    
    if fix_cyc is None:
        
        lim = int(base - 1)
    
    else:
    
        lim = int(base / pow(2.0, lod[did]) - 1)
    
    return lim


def get_limit(test_lim = None, per = None):
        
    if test_lim is None:
            
        intlim = np.iinfo(np.uintc).max
        
    else:
            
        intlim = test_lim
        
    if per is None:
            
        lim = intlim
        
    else:
            
        num_full_per = np.floor_divide( (intlim + 1), per)
            
        lim = (num_full_per * per) - 1
    
    return lim    
    

def nat_cyclic_label(offset, pos, per = None, test_lim = None):
        
    res = 0
        
    lim = get_limit(test_lim, per)
        
    offset = np.mod(offset, lim + 1)
    pos = np.mod(pos, lim + 1)
        
    if (offset + pos) < lim + 1:
        res = offset + pos
    else:
        res = (offset - (lim + 1) ) + pos
    
    return res
    

def int_cyclic_label(offset, pos, per = None, test_lim = None):
        
    res = 0
        
    lim = get_limit(test_lim, per)
        
    offset = np.mod(offset, lim + 1)
    pos = np.mod(pos, lim + 1)
        
    if (offset + pos) < lim + 1:
        res = offset + pos
    else:
        res = (offset - (lim + 1) ) + pos
    
    if res >= (lim + 1) / 2:
        
        res -= (lim + 1)
    
    return res

