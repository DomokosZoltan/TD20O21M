# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 23:03:06 2021

@author: Domokos Zoltán
"""

import mstruct as mst
from copy import deepcopy as deep


class Generator:
    
    def __init__(self, lbgens, ptgens):
        
        self.__lbgens = lbgens        
        
        self.__ptgens = ptgens 
        
        return


    def transform_lb(self, trf_par, dims):
        
        lbgens = deep(self.__lbgens)
        
        for gen in lbgens:
                 
            obj = gen[1]

            curr_dims = gen[0]
            
            if curr_dims == dims:
                
                new_obj = obj.transform(trf_par)
        
        rlbgen = mst.MdimStruct( [new_obj], [dims] )
        
        lbgens.replace(rlbgen, dims, dims[0] )
        
        return lbgens
    
    
    def transform_pt(self, trf_par, dims):
        
        ptgens = deep(self.__ptgens)
        
        for gen in ptgens:
                 
            obj = gen[1]

            curr_dims = gen[0]
            
            if curr_dims == dims:
                
                new_obj = obj.transform(trf_par)
        
        rptgen = mst.MdimStruct( [new_obj], [dims] )
        
        ptgens.replace(rptgen, dims, dims[0] )
        
        return ptgens
    
    
    def get_lod_zero_lb(self):
        
        labels = None
        
        for gen in self.__lbgens:
                 
            obj = gen[1]

            if labels is None:
                     
                labels = obj.generate(obj.par() )
                     
            else:
                     
                labels = labels.prod( obj.generate(obj.par() ) )
                
        return labels
    
    
    def get_lod_zero_pt(self):
        
        points = None
        
        for gen in self.__ptgens:
                 
            obj = gen[1]

            if points is None:
                     
                points = obj.generate(obj.par() )
                     
            else:
                     
                points = points.prod( obj.generate(obj.par() ) )
                
        return points
    

    def prod(self, rhs):
        
        res = None
            
        lbgens = self.__lbgens.prod(rhs.__lbgens)
    
        ptgens = self.__ptgens.prod(rhs.__ptgens)
        
        res = Generator(lbgens, ptgens)
        
        return res
    
    
    def replace(self, rhs, dst_dims, src_offs):        
        
        self.__lbgens.replace(rhs.__lbgens, dst_dims, src_offs)
        
        self.__ptgens.replace(rhs.__ptgens, dst_dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
            
        lbgens = self.__lbgens.restrict(dims)
        
        ptgens = self.__ptgens.restrict(dims)
            
        res = Generator(lbgens, ptgens)
            
        return res
    
    
    def __loc(self, dims, did):
        
        res = did - dims[0]
        
        return res
    
    
    def check(self, incr, bounds, lod, did):
        #TODO
        return 'OK'
 
    
    def pjbnd(self, incr, sel_bnds, lod, points, did):

        pj_bnds = None
    
        dims = self.__lbgens.get_dims(did)
        
        obj = self.__lbgens.get_obj(dims)
        
        labels = points.get_labels()
        
        ldid = self.__loc(dims, did)
        
        rlabels = labels.restrict(dims)
        
        rlod = lod.restrict(dims)
        
        pj_bnds = obj.project_labels(obj.par(), incr, sel_bnds, rlod, rlabels, ldid) 
        
        return pj_bnds
    
    
    def label_trf(self, incr, bounds, lod, points, did):
        
        labels = None
        
        dims = self.__lbgens.get_dims(did)
        
        obj = self.__lbgens.get_obj(dims)
        
        labels = points.get_labels()
        
        ldid = self.__loc(dims, did)
        
        rlabels = labels.restrict(dims)
        
        rlod = lod.restrict(dims)
        
        rlabels = obj.label_trf(obj.par(), incr, bounds, rlod, rlabels, ldid) 
        
        labels.replace(rlabels, dims)
        
        return labels


    def point_trf(self, incr, bounds, lod, base_pt, did):
        
        points = None

        dims = self.__ptgens.get_dims(did)
        
        obj = self.__ptgens.get_obj(dims)
        
        points = base_pt.get_points()
        
        ldid = self.__loc(dims, did)
        
        rpoints = base_pt.get_points().restrict(dims)
        
        rlod = lod.restrict(dims)
        
        rpoints = obj.point_trf(obj.par(), incr, bounds, rlod, rpoints, ldid)
        
        points.replace(rpoints, dims)
        
        return points
    
