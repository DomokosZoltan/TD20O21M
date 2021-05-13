# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:21:44 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep

import mstruct as mst
import codomain as codm
import sglutils as sgl


class BsFcGen:
    
    
    def __init__(self, genstr):
        
        self.genstr = deep(genstr)
        
        return
    
    
    def gen(self, domain, inds):

        res = None
        
        for obj in self.genstr:
            
            dims = obj[0]
            
            genfc = obj[1]

            sl = mst.dims_to_sl(dims)

            dims = (0, dims[1] - dims[0] ) 

            codmstr = mst.MdimStruct(genfc(domain, inds[sl] ), [dims] )

            if res is None:
                
                res = codm.Codomain(codmstr)
                
            else:
                
                res = res.prod( codm.Codomain(codmstr) )
                
        return res
    
    
    def prod(self, rhs):
        
        res = None
        
        res = self.genstr.prod(rhs.genstr)
        
        return res
    
    
    def replace(self, rhs, dims, offs):
        
        self.genstr.replace(rhs, dims, offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
        
        res = self.genstr.restrict(dims)
        
        return res


def get_sphero_gen(src_smb = None, dst_smb = None):
    
    res = None
    
    if src_smb == 'ra' and dst_smb == 'ra':
    
        def generator(src_dm, inds):

            res = None
            
            R = src_dm.bandlims.get_obj( (0,0) )[1]
            A = src_dm.bandlims.get_obj( (1,1) )[1]
            
            nr = src_dm.ns.get_obj( (0,0) )
            na = src_dm.ns.get_obj( (1,1) ) 
            
            dr = R / nr
            da = A / na
            
            shf_r = inds[0]
            shf_a = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            res = [ sgl.spher_stepfc(dr, da, R, A, shf_r, shf_a, phi, ampl, 'ra') ]
            
            return res
    
    elif src_smb == 'sp' and dst_smb == 'sp':
        
        def generator(src_dm, inds):
            
            res = None
            
            S = src_dm.bandlims.get_obj( (0,0) )[1]
            P = src_dm.bandlims.get_obj( (1,1) )[1]
            
            ns = src_dm.ns.get_obj( (0,0) )
            np = src_dm.ns.get_obj( (1,1) ) 
            
            ds = S / ns
            dp = P / np
            
            shf_s = inds[0]
            shf_p = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'f'
            
            sfc = sgl.cyc_stepfc(ds, S, shf_s, phi, ampl, symb)
            
            pfc = sgl.cyc_stepfc(dp, P, shf_p, phi, ampl, symb)
            
            res = [ lambda x : sfc( [ x[0] ] ) * pfc( [ x[1] ] ) ]
            
            return res
    
    elif src_smb == 'sp' and dst_smb == 'ra':
        
        def generator(src_dm, inds):
            
            res = None
            
            S = src_dm.bandlims.get_obj( (0,0) )[1]
            P = src_dm.bandlims.get_obj( (1,1) )[1]
            
            ns = src_dm.ns.get_obj( (0,0) )
            np = src_dm.ns.get_obj( (1,1) )
            
            ds = S / ns
            dp = P / np
            
            R = 1.0 / ds
            A = 1.0 / dp
            
            s = inds[0]
            p = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'ra'
            
            res = [ sgl.spher_expfc(s, p, R, A, phi, ampl, symb) ]
            
            return res
    
    else:
        
        def generator(src_dm, inds):
            
            res = None
            
            R = src_dm.bandlims.get_obj( (0,0) )[1]
            A = src_dm.bandlims.get_obj( (1,1) )[1]
            
            nr = src_dm.ns.get_obj( (0,0) ) 
            na = src_dm.ns.get_obj( (1,1) ) 
            
            dr = R / nr
            da = A / na
            
            S = 1.0 / dr
            P = 1.0 / da
            
            r = inds[0]
            a = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'f'
            
            sfc = sgl.expfc(S, r, phi, ampl, symb)
            
            pfc = sgl.expfc(P, a, phi, ampl, symb)
            
            res = [ lambda x : sfc( [ x[0] ] ) * pfc( [ x[1] ] ) ]
            
            return res
        
    res = generator
    
    return res


def get_codm_gen(src_smb = None, dst_smb = None):
    
    res = None
    
    if src_smb == 't' and dst_smb == 't':
        
        def generator(src_dm, ind):
            
            res = None
            
            T = src_dm.bandlims.get_obj( (0,0) )[1]
            n = src_dm.ns.get_obj( (0,0) ) 
            dt = T/n
            shf = ind[0]
            phi = 1.0 + 0.0j
            ampl = 1.0
            symb = 't'
            
            res = [ sgl.cyc_stepfc(dt, T, shf, phi, ampl, symb) ]
            
            return res
        
    elif src_smb == 'f' and dst_smb == 'f':
        
        def generator(src_dm, ind):
            
            res = None
            
            F = src_dm.bandlims.get_obj( (0,0) )[1]
            n = src_dm.ns.get_obj( (0,0) )
            df = F/n
            shf = ind[0]
            phi = 1.0 + 0.0j
            ampl = 1.0
            symb = 'f'
            
            res = [ sgl.cyc_stepfc(df, F, shf, phi, ampl, symb) ]
            
            return res
    
    elif src_smb == 't' and dst_smb == 'f':
    
        def generator(src_dm, ind):
            
            res = None
            
            T = src_dm.bandlims.get_obj( (0,0) )[1]
            n = src_dm.ns.get_obj( (0,0) ) 
            dt = T/n
            F = 1.0 / dt
            t = ind[0]
            phi = 1.0 + 0.0j
            ampl = 1.0
            symb = 'f'
            
            res = [ sgl.expfc(F, t, phi, ampl, symb) ]
            
            return res
    
    else:
    
        def generator(src_dm, ind):
            
            res = None
            
            F = src_dm.bandlims.get_obj( (0,0) )[1]
            n = src_dm.ns.get_obj( (0,0) )
            df = F/n
            T = 1.0 / df
            f = ind[0]
            phi = 1.0 + 0.0j
            ampl = 1.0
            symb = 't'
            
            res = [ sgl.expfc(T, f, phi, ampl, symb) ]
            
            return res
        
    res = generator
    
    return res


def get_bess_gen(src_smb = None, dst_smb = None):
    
    res = None
    
    if src_smb == 'ra' and dst_smb == 'ra':
    
        def generator(src_dm, inds):

            res = None
            
            R = src_dm.bandlims.get_obj( (0,0) )[1]
            A = src_dm.bandlims.get_obj( (1,1) )[1]
            
            nr = src_dm.ns.get_obj( (0,0) )
            na = src_dm.ns.get_obj( (1,1) ) 
            
            dr = R / nr
            da = A / na
            
            shf_r = inds[0]
            shf_a = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            res = [ sgl.spher_stepfc(dr, da, R, A, shf_r, shf_a, phi, ampl, 'ra') ]
            
            return res
    
    elif src_smb == 'sp' and dst_smb == 'sp':
        
        def generator(src_dm, inds):
            
            res = None
            
            S = src_dm.bandlims.get_obj( (0,0) )[1]
            P = src_dm.bandlims.get_obj( (1,1) )[1]
            
            ns = src_dm.ns.get_obj( (0,0) )
            np = src_dm.ns.get_obj( (1,1) ) 
            
            ds = S / ns
            dp = P / np
            
            shf_s = inds[0]
            shf_p = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'f'
            
            sfc = sgl.stepfc(ds, shf_s, phi, ampl, symb)
            
            pfc = sgl.stepfc(dp, shf_p, phi, ampl, symb)
            
            res = [ lambda x : sfc( [ x[0] ] ) * pfc( [ x[1] ] ) ]
            
            return res
    
    elif src_smb == 'sp' and dst_smb == 'ra':
        
        def generator(src_dm, inds):
            
            res = None
            
            S = src_dm.bandlims.get_obj( (0,0) )[1]
            P = src_dm.bandlims.get_obj( (1,1) )[1]
            
            ns = src_dm.ns.get_obj( (0,0) )
            np = src_dm.ns.get_obj( (1,1) )
            
            ds = S / ns
            dp = P / np
            
            R = 1.0 / ds
            A = 1.0 / dp
            
            s = inds[0]
            p = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'ra'
            
            res = [ sgl.spher_bessfc(s, p, R, A, phi, ampl, symb) ]
            
            return res

    else:
        
        def generator(src_dm, inds):
            
            res = None
            
            R = src_dm.bandlims.get_obj( (0,0) )[1]
            A = src_dm.bandlims.get_obj( (1,1) )[1]
            
            nr = src_dm.ns.get_obj( (0,0) ) 
            na = src_dm.ns.get_obj( (1,1) ) 
            
            dr = R / nr
            da = A / na
            
            S = 1.0 / dr
            P = 1.0 / da
            
            r = inds[0]
            a = inds[1]
            
            phi = 1.0 + 0.0j
            ampl = 1.0
            
            symb = 'f'
            
            sfc = sgl.rev_bessfc(r * dr, S, phi, ampl, symb)
            
            pfc = sgl.expfc(P, a, phi, ampl, symb)
            
            res = [ lambda x : sfc( [ x[0] ] ) * pfc( [ x[1] ] ) ]
            
            return res
        
    res = generator
    
    return res