# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 12:51:02 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep

import mstruct as mst
import labels as lb
import cycutils as cyu


class LbFixParams:

    def __init__(self, card, lim):
        
        self.lim = lim[0]
        
        self.card = card[0]
        
        return
    
    
class LbModParams:
    
    def __init__(self, mcard, offs, per):
        
        self.mcard = mcard[0]
        
        self.offs = offs[0]
        
        self.per = per[0] 
        
        return
    

class PtFixParams:

    def __init__(self, t, shf, band):
        
        self.T = band[0]
        
        self.t = deep(t)
        
        self.per = len(t)
        
        self.shf_t = shf[0]
        
        return
    
    
class PtModParams:
    
    def __init__(self, sc):
        
        self.sc_t = sc[0]
        
        return


#BOUNDARY
def bound_check(self, par, lod, new_bounds, restr_base, loc_did):
    #TODO
    return 'ok'


#DEPTH
def depth_check(self, par, lod, incr, bounds, restr_base, loc_did):
    #TODO
    return 'ok'


class LbGen:
    
    def __init__(self, fix_params):
        
        self.__fix_params = fix_params
        
        return


    def par(self):
        
        res = deep(self.__fix_params)
        
        return res


    def transform(self, trf_par):
        
        res = deep(self)
        
        res.__fix_params.lim = trf_par.lim
        
        res.__fix_params.card = trf_par.card
        
        return res


    def generate(self, fix_par, mod_par = None):
        
        res = None
    
        lim = fix_par.lim
        
        if mod_par is None:
        
            offs = 0
            card = fix_par.card
            per = card
        
        else:
            
            offs = mod_par.offs
            card = mod_par.mcard
            per = mod_par.per
            
        def labfc(x):    
            
            clab = cyu.int_cyclic_label(offs, x[0], per, lim)
            
            res = [ clab ]
            
            return res
         
        fcdims = [ (0, 0) ]
        
        cards = mst.MdimStruct( [card], fcdims )
        labels = mst.MdimStruct( [labfc], fcdims ) 
        
        res = lb.LabelSet(cards, labels)
        
        return res


    def label_trf(self, par, incr, new_bounds, rlod, rlabels, loc_did):
        
        res = None
        
        if incr is None:
            
            res = self.bound_labels(par, new_bounds, rlod, rlabels, loc_did)
            
        elif new_bounds is None:
            
            res = self.depth_labels(par, incr, rlod, rlabels, loc_did)
            
        return res
    

    def bound_labels(self, par, new_bounds, rlod, rlabels, loc_did):
    
        labels = None
        
        offs = [ new_bounds[0] ]
            
        per = [ rlabels.part_card(loc_did) ]
            
        card = [ 1 + abs(new_bounds[1] - new_bounds[0] ) ]
            
        mod_par = LbModParams(card, offs, per)
        
        labels = self.generate(par, mod_par)
    
        return labels
    
    
    def project_labels(self, fix_par, incr, sel_bnds, rlod, rlabels, loc_did):
        
        res = None
        
        fac = pow(2, incr)
        
        start = sel_bnds[0]
        
        end = sel_bnds[1]
        
        if 1 == end % 2:
            
            end += 1
            
        else:
            #TODO: raise
            pass
        
        start = fac * start
        
        end = fac * end
        
        end -= 1
        
        res = [start, end]
        
        return res
        
    
    def depth_labels(self, fix_par, incr, rlod, rlabels, loc_did): 
          
        labels = None
        
        lod = rlod.get_obj( (0, 0) )
        
        fac = pow(2, lod + incr)
        
        card = [ fac * fix_par.card]
            
        per = [ fac * fix_par.card]
        
        offs = [0]
        
        mod_par = LbModParams(card, offs, per)
        
        labels = self.generate(fix_par, mod_par)
        
        return labels


class PtGen:
    
    
    def __init__(self, fix_params):
        
        self.__fix_params = fix_params
        
        return  
    
    
    def par(self):
        
        res = deep(self.__fix_params)
        
        return res
    
    
    def transform(self, trf_par):
        
        res = deep(self)
        
        res.__fix_params.T = trf_par.T
        
        res.__fix_params.dt = trf_par.dt
        
        res.__fix_params.shf = trf_par.shf_t
        
        return res
    
    
    def generate(self, fix_par, mod_par = None):
         
        res = None
        
        shf_t = fix_par.shf_t
        
        t = fix_par.t
        p = fix_par.per
        T = fix_par.T
        
        if mod_par != None:
            
            sc_t = mod_par.sc_t
       
        else:
        
            sc_t = 1.0
            
        def linearfc(inds):
            
            x = inds[0] % p
            s = inds[0] // p
            
            rt = sc_t * ( t[x] + (s * T + shf_t) )
                        
            return [rt]
           
        res = mst.MdimStruct( [linearfc], [ (0,0) ] )   
        
        return res


    def point_trf(self, par, incr, new_bounds, rlod, rlabels, loc_did):
        
        res = None
        
        if incr is None:
            
            res = self.bound_points(par, new_bounds, rlod, rlabels, loc_did)
            
        elif new_bounds is None:
            
            res = self.depth_points(par, incr, rlod, rlabels, loc_did)
            
        return res
    
    
    def bound_points(self, par, new_bounds, lod, rpoints, loc_did):
        
        res_points = None
    
        res_points = deep(rpoints)
        
        return res_points
    
    
    def depth_points(self, fix_par, incr, rlod, rpoints, loc_did):
        
        res_points = None
    
        sc = [ pow(2.0, -lod[1] ) for lod in rlod]
        
        sc[loc_did] *= pow(2.0, -incr)
    
        mod_par = PtModParams(sc)
        
        res_points = self.generate(fix_par, mod_par)
        
        
        return res_points