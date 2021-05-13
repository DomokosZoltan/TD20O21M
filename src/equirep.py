# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 19:54:57 2021

@author: Domokos Zoltán
"""

import copy as cp
import numpy as np

from math import isnan
from math import isinf

import mstruct as mst
import labels as ls
import sglspace as sp
import domain as dm
import sglgen as sgg

import view as vw #TODO: this is for debugging, remove later

#These checks are always the same, regardless of the rules.
def icheck_X_n(X,n):
    eps=0.0000001
    if isnan(X) or isnan(float(n)):
        raise ValueError('NaN input')
    elif 10000 < X or isinf(float(n)):
        raise ValueError('inf input')
    elif np.isclose(X,0.0,eps) or\
            0==int(n):
        raise ValueError('input close to zero')
    elif 0.0>X or 0>int(n):
        raise ValueError('negative input')
    elif not float(n).is_integer():
        raise ValueError('n is not integer')
    elif 5000<int(n):
        raise ValueError('n too high')
    elif np.mod(int(n),2) != 0:
        #n is not even
        raise ValueError('n is odd')
    else:
        pass
    return
    
    
def icheck_X_Y(X,Y):
    eps=0.0000001
    if isnan(X) or isnan(Y):
        raise ValueError('NaN input')
    elif 10000 < X or 10000 < Y:
        raise ValueError('Inf. product')
    elif np.isclose(X,0.0,eps) or\
        np.isclose(Y,0.0,eps):
        raise ValueError('input close to zero')
    elif 0.0>X or 0.0>Y:
        raise ValueError('input negative')
    else:
        pass
    return
    
    
def ocheck_n(n):
    if isnan(float(n)):
        raise ValueError('NaN n')
    elif isinf(float(n)):
        raise ValueError('inf n')
    elif 0==int(n):
        raise ZeroDivisionError #TODO: fix this
    elif 0>int(n):
        raise ValueError('zero or negative n')
    elif 5000<int(n):
        raise ValueError('n too high') #TODO: define a config file for settings
    elif not float(n).is_integer():
        raise ValueError('n is not integer')
    elif np.mod(int(n),2) != 0:
        #n is not even
        raise ValueError('n is odd')
    else:
        pass
    return n
    
    
def ocheck_Y(Y):
    eps=0.0000001
    if isnan(Y):
        raise ValueError('NaN Y')
    elif 10000 < Y:
        raise ValueError('inf Y')
    elif 0.0>Y:
        raise ValueError('zero or negative Y')
    elif np.isclose(Y,0.0,eps):
        raise ValueError('almost zero Y')
    else:
        pass
    return Y
    

def get_XX_rule(src_smb = None, fix_smb = None, did = 0, core_symbs = ['t'] ):
    
    def rule_XX(src_val, fix_val):
            
        icheck_X_Y(src_val, fix_val)
            
        n = src_val * fix_val
            
        if src_smb == core_symbs[did]:
                
            bandlims = [0.0, src_val]
                
        elif fix_smb in core_symbs:
                
            bandlims = [0.0, fix_val]
            
        ocheck_n(n)
            
        return [n, bandlims]
        
    res = rule_XX

    return res


def get_nX_rule(fix_smb = None, did = 0, core_symbs = ['t'] ):
    
    def rule_nX(src_val, fix_val):
            
        icheck_X_n(fix_val, src_val)
            
        n = src_val
            
        if fix_smb == core_symbs[did]:
                
            bandlims = [0.0, fix_val]
                
        else:
                
            bandlims = [0.0, (n / fix_val) ]
            
        ocheck_Y(bandlims[1])
            
        return [n, bandlims]
        
    res = rule_nX
        
    return res


def get_Xn_rule(src_smb = None, did = 0, core_symbs = ['t'] ):
    
    def rule_Xn(src_val, fix_val):
            
        icheck_X_n(src_val, fix_val)
            
        n = fix_val
            
        if src_smb == core_symbs[did]:
                
            bandlims = [0.0, src_val]
                
        else:
                
            bandlims = [0.0, (n / src_val) ]
            
        ocheck_Y(bandlims[1])
            
        return [n, bandlims]
        
    res = rule_Xn
        
    return res


#TODO: move these functions to sglgen
def get_cartes_gen(src_smb = None, dst_smb = None):
    
    def get_dst_dm(src_dm):

        X = src_dm.bandlims.get_obj( (0,0) )[1]
       
        n = src_dm.ns.get_obj( (0,0) )
            
        dx = X / n      
        
        Y = 1.0 / dx
        
        dst_bandlims = [ [0, Y] ]
        dst_ns = [n]
        
        dst_dm = dm.factory(dst_bandlims, dst_ns)
        
        return dst_dm
    
    return get_dst_dm


def get_polar_gen(src_smb = None, dst_smb = None):
    
    fc = None
    
    #TODO: merge cases, if possible
    if 'ra' == src_smb and 'sp' == dst_smb:
    
        def get_dst_dm(src_dm):
            
            R = src_dm.bandlims.get_obj( (0,0) )[1]
            A = src_dm.bandlims.get_obj( (1,1) )[1]
       
            nr = src_dm.ns.get_obj( (0,0) )
            na = src_dm.ns.get_obj( (1,1) )
            
            dr = R / nr 
            da = A / na 
        
            S = 1.0 / dr
            P = 1.0 / da
        
            dst_dm = dm.factory( [ [0, S] ], [nr], [0.0], 'cartes')
            
            dst_dm = dst_dm.prod( dm.factory( [ [0, P] ], [na], [0.0], 'cartes') )
            
            return dst_dm
        
        fc = get_dst_dm
        
    elif 'sp' == src_smb and 'ra' == dst_smb:
        
        def get_dst_dm(src_dm):
            
            S = src_dm.bandlims.get_obj( (0,0) )[1]
            P = src_dm.bandlims.get_obj( (1,1) )[1]
       
            ns = src_dm.ns.get_obj( (0,0) )
            np = src_dm.ns.get_obj( (1,1) )
            
            ds = S / ns 
            dp = P / np 
        
            R = 1.0 / ds
            A = 1.0 / dp
        
            dst_bandlims = [ [0, R], [0, A] ]
            dst_ns = [ns, np]
            
            dst_dm = dm.factory(dst_bandlims, dst_ns, [0.0000001, 0.0], 'polar')
            
            return dst_dm
        
        fc = get_dst_dm
        
    else:
        #TODO: raise exception
        pass
    
    return fc


class Table:

    def __init__(self, domain = None, symb = None, words = None,
                         rls = None, codmgens = None, dmgens = None, depr = False):
        
        if not symb is None:
        
            self.core_symb = symb
            
        else:
            
            self.core_symb = mst.MdimStruct( ['t'], [ (0,0) ] )
            
        if not words is None:
            
            self.core_words = words
            
        else:
            
            self.core_words = cp.deepcopy(self.core_symb)
        
        if not domain is None:
            
            self.core_dm = domain
        
        else:
            
            self.core_dm = dm.Domain()
            
        if not rls is None:
            
                self.rls = cp.deepcopy(rls)
       
        else:

            self.rls = mst.MdimStruct( [ {'t' : {'t' : None,
                                                 'f' : get_XX_rule('t', 'f'),
                                                 'n': get_Xn_rule('t') },
                                            
                                          'f' : {'t' : get_XX_rule('f', 't'),
                                                 'f' : None, 
                                                 'n': get_Xn_rule('f') },
                                            
                                          'n' : {'t' : get_nX_rule('t'),
                                                 'f' : get_nX_rule('f'), 
                                                 'n': None } } ],
    
                                        [ (0,0) ] )

        #Base function codomain generators
        if not ( codmgens is None ):
            
                self.gens = cp.deepcopy(codmgens)
                
        else:
    
            self.gens = mst.MdimStruct( [ {'t' : {'t' : sgg.get_codm_gen('t', 't'),
                                                  'f' : sgg.get_codm_gen('f','t') },
    
                                           'f' : {'t' : sgg.get_codm_gen('t', 'f'),
                                                  'f' : sgg.get_codm_gen('f', 'f') } } ],
    
                                        [ (0,0) ] )
        
        #Base function domain generators
        if not dmgens is None:
        
            self.dmgens = dmgens
            
        else:

            self.dmgens = mst.MdimStruct( [ {'t' : {'t' : None,
                                                    'f' : get_cartes_gen('t', 'f') },
    
                                             'f' : {'t' : get_cartes_gen('f', 't'),
                                                    'f' : None } } ],
    
                                          [ (0,0) ] )

            
        cards = cp.deepcopy(self.core_dm.ns)        
        self.bsubsp_labels = None
        for card in cards:
            
            new_labels = ls.LabelSet(mst.MdimStruct( [ card[1] ], [ (0,0) ] ) )
            
            if self.bsubsp_labels is None:
                
                self.bsubsp_labels = new_labels
                
            else:
                
                self.bsubsp_labels = self.bsubsp_labels.prod(new_labels)
    
        #TODO: remove this deprecated member, as soon as tests will be adapted
        #to new class structure.
        if depr:
        
            dmgen = self.dmgens.get_obj((0,0))['t']['f']
            
            self.domains = {'t' : self.core_dm, 'f' : dmgen(self.core_dm) }
    
        return
    
    
    #gets values and symbols
    def update(self, src_val, did, cmd, depr = False):
        
        dims = (did, did)
        
        #TODO: implement these as function arguments 
        src = cmd['src']
        fix = cmd['fix']
        
        #1. Get fix value
        if 'n' == fix[did]:
        
            fix_val = self.core_dm.n(did)
        
        elif self.core_symb.get_obj(dims) == fix[did]:
        
            fix_val = self.core_dm.bandlims.get_obj(dims)[1]
        
        else:
            
            wdims = self.core_words.get_dims(did)           
            crword = self.core_words.get_obj(wdims)
                
            dmgen = self.dmgens.get_obj(wdims)[crword][fix]
            fix_dm = dmgen(self.core_dm)
            
            fix_val = fix_dm.bandlims.get_obj(dims)[1]
        
        #2. Apply rule
        rule = self.rls.get_obj(dims)[src[did] ][fix[did] ]
        dst_val = rule(src_val, fix_val)
        
        #3. Prepare results
        new_n = dst_val[0]
        new_bandlims = cp.deepcopy(dst_val[1] )
        
        #4. Core update
        self.core_dm.update(new_n, new_bandlims, did)
                
        #TODO: find final sol. this is a workaround
        rlabels = ls.LabelSet(mst.MdimStruct( [new_n], [ (0, 0) ] ) )
        self.bsubsp_labels.replace(rlabels, dims)
        
        #TODO: remove these depr. members and adapt test cases in test_equirep
        if depr:
            
            dmgen = self.dmgens.get_obj((0,0))['t']['f']
            self.domains = {'t' : self.core_dm, 'f' : dmgen(self.core_dm) }
        
        return


    def gen_domain(self, symbstr):
        
        res = None
        
        for obj in symbstr:
            
            dims = obj[0]
            symb = obj[1]
            
            rcore_dm = self.core_dm.restrict(dims)
            
            bsymb = self.core_words.get_obj(dims)
            
            if bsymb != symb:
                
                dmgen = self.dmgens.get_obj(dims)[bsymb][symb]
                next_dm = dmgen(rcore_dm)
            
            else:
                
                next_dm = rcore_dm
            
            if res is None:
                
                res = next_dm
            
            else:
                
                res = res.prod(next_dm)

        return res
            
        
    def get_bsgen(self, srcstr, dststr):
        
        res = None
        
        res_genstr = None
        
        for itsrc in srcstr:
            
            dims = itsrc[0]
            
            src = itsrc[1]
            dst = dststr.get_obj(dims)
            
            loc_dims = (0, dims[1]-dims[0] ) 
            
            codmgen = self.gens.get_obj(dims)[src][dst]
            next_genstr = mst.MdimStruct( [codmgen], [loc_dims] )
            
            if res_genstr is None:
                
                res_genstr = next_genstr
                
            else:
                
                res_genstr = res_genstr.prod(next_genstr)
                
        res = sgg.BsFcGen(res_genstr)

        return res


    def base_change(self, srcstr, dststr, curr_src_dm, curr_dst_dm, signal, md = 'cartes'):
        
        res = None

        sigvec = sp.sigvec(signal.fcval_at, curr_src_dm, curr_dst_dm)[1]
        
        sigvec = sp.scalev(sigvec, curr_dst_dm, 'sqrt(avg_dx)')
    
        #flatten the vectorized signal
        sigvec = sigvec.flatten('F')
        
        #get the base function generator
        bsgen = self.get_bsgen(srcstr, dststr)
        
        dst_dm = self.gen_domain(dststr)
        
        #iterator for the base function set labels
        it_subsp = self.bsubsp_labels.__iter__('F')
        
        #prepare inverse trf. mtx., this will be filled with base functions
        mtxshp = (sigvec.shape[0], sigvec.shape[0] )
        itrf_mtx = np.zeros(shape = mtxshp, dtype = np.complex64)
        
        done = False
        while not done:
            
            try:
                
                inds = it_subsp.nd_index()
                flind = it_subsp.flat_index()
                
                next(it_subsp)
                
                bsfc = bsgen.gen(dst_dm, inds)
                
                bsigvec = sp.cdmvec(bsfc.fcval_at, curr_dst_dm)[1]
                
                #vw.plot_heatvec(curr_dst_dm, bsigvec, md)
                
                bsigvec = sp.scalev(bsigvec, curr_dst_dm, 'sqrt(avg_dx)')
       
                #flatten base
                bsigvec = bsigvec.flatten('F')
                
                #copy base vector in matrix
                itrf_mtx[:, flind] = np.transpose(bsigvec)
                
            except StopIteration:
                
                done = True
        
        #Solve equation by inverting itrf_mtx, it must be an invertible transform.
        trf_mtx = np.linalg.inv(itrf_mtx)    
        
        flat_coeffs = np.dot(trf_mtx, sigvec)
        
        #Assemble the reconstructed Signal using the flat_coeffs
        
        #TODO: this is a w.around, replace it with a final sol. (from here)
        ndim_shape = []
        for n in dst_dm.ns:
            ndim_shape.append( n[1] )
        #(to here)
        
        labels = self.bsubsp_labels
        
        coeffs = np.reshape(flat_coeffs, ndim_shape, 'F')
        
        bsgen = self.get_bsgen(dststr, dststr)
        
        res = sp.Signal(labels, coeffs, bsgen)
        
        return res


    def product(self, rhs):
        
        res = None
        
        crdm = self.core_dm.prod(rhs.core_dm)
        
        crsymbs = self.core_symb.prod(rhs.core_symb)
        
        crwords = self.core_words.prod(rhs.core_words)
        
        rls = self.rls.prod(rhs.rls)
        
        codmgens = self.gens.prod(rhs.gens)
        
        dmgens = self.dmgens.prod(rhs.dmgens)
        
        res = Table(crdm, crsymbs, crwords, rls, codmgens, dmgens)
        
        return res
    
    
    def replace(self, rhs, dims, src_offs):
        
        self.core_dm.replace(rhs.core_dm, dims, src_offs)
        
        self.core_symb.replace(rhs.core_symb, dims, src_offs)
        
        self.core_words.replace(rhs.core_words, dims, src_offs)
        
        self.rls.replace(rhs.rls, dims, src_offs)
        
        self.gens.replace(rhs.gens, dims, src_offs)
        
        self.dmgens.replace(rhs.dmgens, dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
         
        res = None
        
        crdm = self.core_dm.restrict(dims)
        
        crsymbs = self.core_symb.restrict(dims)
        
        crwords = self.core_words.restrict(dims)
        
        rls = self.rls.restrict(dims)
        
        codmgens = self.gens.restrict(dims)
        
        dmgens = self.dmgens.restrict(dims)
        
        res = Table(crdm, crsymbs, crwords, rls, codmgens, dmgens)
        
        return res
    
    