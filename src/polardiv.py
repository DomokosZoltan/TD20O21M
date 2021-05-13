# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 19:30:10 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep
import numpy as np
   
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

    def __init__(self, diff, shf, band):
        
        self.A = 2.0 * np.pi
        self.R = band[0]
        
        self.dr = diff[0]
        self.da = diff[1]
        
        self.shf_r = shf[0]
        self.shf_a = shf[1]
        
        return
    
    
class PtModParams:
    
    def __init__(self, sc):
        
        self.sc_r = sc[0]
        self.sc_a = sc[1]
        
        return


#BOUNDARY
def bound_check(self, par, lod, new_bounds, restr_base, loc_did):
    #TODO
    return 'ok'


#DEPTH
def depth_check(self, par, lod, incr, bounds, restr_base, loc_did):
    #TODO
    return 'ok'


class PolarLbGen:
    
    def __init__(self, fix_params):
        
        self.__fix_params = fix_params
        
        return


    def par(self):
        
        res = deep(self.__fix_params)
        
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
    
    
    def project_labels(self, par, incr, lod, restr_base, loc_did):
        #TODO
        labels = None
        
        return labels
        
    
    def depth_labels(self, fix_par, incr, rlod, rlabels, loc_did): 
          
        labels = None
        
        lod = rlod.get_obj( (0, 0) )
        
        curr_card = rlabels.part_card(0)
        
        card = [ pow(2, lod + incr) * curr_card ]
            
        per = [ pow(2, lod + incr) * curr_card ]
        
        offs = [0]
        
        mod_par = LbModParams(card, offs, per)
        
        labels = self.generate(fix_par, mod_par)
        
        return labels


class PolarPtGen:
    
    
    def __init__(self, fix_params):
        
        self.__fix_params = fix_params
        
        return  
    
    
    def par(self):
        
        res = deep(self.__fix_params)
        
        return res
    
    
    def generate(self, fix_par, mod_par = None):
         
        res = None
        
        shf_r = fix_par.shf_r
        shf_a = fix_par.shf_a
        
        dr = fix_par.dr
        da = fix_par.da
          
        if mod_par != None:
            
            sc_r = mod_par.sc_r
            sc_a = mod_par.sc_a
       
        else:
        
            sc_r = 1.0
            sc_a = 1.0
        
        def polarfc(inds):
                        
            r = inds[0] * (sc_r * dr) + (sc_r * shf_r)
            a = inds[1] * (sc_a * da) + (sc_a * shf_a)
                        
            x = r * np.cos(a)
            y = r * np.sin(a)
                        
            return [x, y]
           
        res = mst.MdimStruct( [polarfc], [ (0,1) ] )   
        
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
    
        sc = [ pow(-2.0, lod[1] ) for lod in rlod]
        
        sc[loc_did] *= pow(-2.0, incr)
    
        mod_par = PtModParams(sc)
        
        res_points = self.generate(fix_par, mod_par)
        
        
        return res_points

