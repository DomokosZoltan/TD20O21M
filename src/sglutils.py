# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:44:05 2021

@author: Domokos Zoltán
"""

import numpy as np
import cmath as cm

from scipy.special import j0, j1
from scipy.special import jn_zeros


def expfc(X, f, phi, ampl, symb):
    
    s = 1.0j
    if 'f' == symb:
        s = -1.0j
    
    w = 2.0 * np.pi / X
    c = phi * ampl
    
    res = lambda x : c * np.sqrt(1.0 / X) * cm.exp( s * w * f * x[0])
    
    return res


def stepfc(dx, shf, phi, ampl, symb):
    
    c = phi * ampl
    
    res = lambda x: c * np.sqrt(1.0 / dx)\
                    if (shf * dx) <= x[0] and x[0] < ( (shf + 1) * dx)\
                    else 0.0 + 0.0j
    
    return res


def cyc_stepfc(dx, X, shf, phi, ampl, symb):
    
    c = phi * ampl
    
    res = lambda x: c * np.sqrt(1.0 / dx)\
                    if (shf * dx) <= np.remainder(x[0], X) and np.remainder(x[0], X) < ( (shf + 1) * dx)\
                    else 0.0 + 0.0j
                    
    return res


def impfc(dx, shf, phi, ampl, symb):
    
    eps = 0.00001
    
    c= phi * ampl    
    
    res = lambda x: c * np.sqrt(1.0 / dx)\
                    if np.abs(1.0 - ( (shf * dx) / x) ) < eps\
                    else 0.0 + 0.0j
    
    return res
    

def spher_bessfc(s, fa, R, A, phi, ampl, symb, alph = 0.0):
    
    zer = jn_zeros(0, 16)[s]
    
    def fc(pt):
        
        res = None
        
        #pt is given in cartesian coords. so we have to aquire radius and angle
        r = np.sqrt(pt[0] ** 2.0 + pt[1] ** 2.0) 
        
        if r == 0.0:
         
            a = 0.0
        
        else:
             
            if 0 <= pt[1]:
            
                a = np.arccos( pt[0] / r )
                
            else:
                
                a = np.arccos( - pt[0] / r ) + np.pi
    
        c = phi * ampl
    
        wx = 1.0
        
        norm = 1.0 # / 32.0
        
        def wave_rad(x): 
            
            u = x * zer
            
            res = norm * (2.0 * x) / ( j1(zer) ** 2.0 ) * j0(u)
            
            return res
            
        wave_ang = lambda x : c * cm.exp(wx * 1.0j * fa * x)
        
        res =  wave_rad( r ) * wave_ang( a )
        
        return res
    
    return fc    


def rev_bessfc(y, S, phi, ampl, symb):
    
    zer = jn_zeros(0, 16) 
    
    res = lambda x : j0( y * zer[int(x[0])] )
    
    return res


def spher_expfc(fr, fa, R, A, phi, ampl, symb):
    
    s = 1.0j
    if 'f' == symb:
        s = -1.0j
    
    def fc(pt):
        
        res = None
        
        #pt is given in cartesian coords. so we have to aquire radius and angle
        r = np.sqrt(pt[0] ** 2.0 + pt[1] ** 2.0) 
        
        if r == 0.0:
         
            a = 0.0
        
        else:
             
            if 0 <= pt[1]:
            
                a = np.arccos( pt[0] / r )
                
            else:
                
                a = np.arccos( -pt[0] / r ) + np.pi
    
        wr = 2.0 * np.pi / R
        
        c = phi * ampl
    
        if 0.0 == r:
            
            wave_rad = lambda x : c * np.sqrt(1.0 / R) * cm.exp( s * wr * fr * x)
            
            wave_ang = lambda x : 1.0 + 0.0j
            
        else:
        
            wx = 1.0
            
            wave_rad = lambda x : c * np.sqrt(1.0 / R) * cm.exp( s * wr * fr * x)
            
            wave_ang = lambda x : c * cm.exp( s * wx * fa * x)
        
        res =  wave_rad( r ) * wave_ang( a )
        
        return res
    
    return fc


def spher_stepfc(dr, da, R, A, shf_r, shf_a, phi, ampl, symb):
    
    c = phi * ampl
    
    def fc(pt):
        
        #pt is given in cartesian coords. so we have to aquire radius and angle
        r = np.sqrt(pt[0] ** 2.0 + pt[1] ** 2.0) 
        
        if r == 0.0:
         
            a = 0.0
        
        else:
             
            if 0 <= pt[1]:
            
                a = np.arccos( pt[0] / r )
                
            else:
                
                a = np.arccos( -pt[0] / r ) + np.pi
        
        eps = 0.000001 #TODO: remove this wk.around
        
        step_rad = lambda x: c if (shf_r * dr) <= x + eps and x < ( (shf_r + 1) * dr) else 0.0 + 0.0j
                    
        step_ang = lambda x: c if (shf_a * da) <= x + eps and x < ( (shf_a + 1) * da) else 0.0 + 0.0j
        
        dx = dr * da
        
        res = np.sqrt(1.0 / dx) * step_rad( r ) * step_ang( a )
        
        return res
    
    return fc


def rectfc():
    #TODO
    pass


def sincfc():
    #TODO
    pass