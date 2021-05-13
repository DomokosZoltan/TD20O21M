# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:44:45 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep

import mstruct as mst
import points as pt


class InvalidRequest(Exception):
    pass


#TODO: proper description for Grid
class Grid:
    
    def __init__(self, gen):
        #generators       
        self.__gen = gen
        
        #points
        labels = self.__gen.get_lod_zero_lb()
        points = self.__gen.get_lod_zero_pt()
        
        self.__base_points = pt.PointSet(labels, points)
        self.__curr_points = pt.PointSet(labels, points)
         
        #dims will be used for LODs and imits
        dimspan = range(0, self.__base_points.dims() )
        all_dims = [ (did, did) for did in dimspan ]
        
        #LODs
        lod_objs = [0 for did in all_dims]
        
        self.__lod = mst.MdimStruct(lod_objs, all_dims)
        
        #Limits - start
        it_base = self.__base_points.get_labels().__iter__('ND')
        
        start_objs = deep( it_base.nd_index() )
        
        self.__curr_abs_start = mst.MdimStruct(start_objs, all_dims)
            
        #Limits - end
        it_base.set_last()
            
        end_objs = deep( it_base.nd_index() )
        
        self.__curr_abs_end = mst.MdimStruct(end_objs, all_dims)
                
        return
    
    
    def __iter__(self, order='F', which='bs'):
        
        res_it = None
        
        if which == 'bs':
            
            res_it = self.__base_points.get_labels().__iter__(order)
        
        elif which == 'curr':
        
            res_it = self.__curr_points.get_labels().__iter__(order)
            
        else:
            #TODO: raise            
            pass
        
        return res_it
    
    
    def prod(self, rhs_grid, cplims = True):
        
        gen = self.__gen.prod(rhs_grid.__gen)
        
        res_grid = Grid(gen)        
       
        start = self.__curr_abs_start.prod(rhs_grid.__curr_abs_start)
        
        end = self.__curr_abs_end.prod(rhs_grid.__curr_abs_end)
        
        if cplims:
        
            dimspan = range(0, res_grid.__base_points.dims() )
            for did in dimspan:            
                
                dims = (did, did)
                
                select_bounds = [start.get_obj(dims), end.get_obj(dims) ]            
                
                res_grid.set_curr(0, select_bounds, did)
        
        return res_grid
    
    
    def replace(self, new_grid, dst_dims, src_offs = 0):
        
        #replace mstructs
        self.__gen.replace(new_grid.__gen, dst_dims, src_offs)
        
        #replace points
        self.__base_points.replace(new_grid.__base_points, dst_dims, src_offs)
        
        self.__curr_points.replace(new_grid.__curr_points, dst_dims, src_offs)
        
        #Set lods
        self.__lod.replace(new_grid.__lod, dst_dims, src_offs)
        
        #Finally the iterators
        self.__curr_abs_start.replace(new_grid.__curr_abs_start, dst_dims, src_offs)
        
        self.__curr_abs_end.replace(new_grid.__curr_abs_end, dst_dims, src_offs)
        
        return
    

    def restrict(self, dims):
        
        gen = self.__gen.restrict(dims)
        
        res = Grid(gen)
        
        #TODO: fix this workaround:  
        res.__curr_points = self.__curr_points.restrict(dims)
        
        res.__lod = self.__lod.restrict(dims)
        
        res.__curr_abs_start = self.__curr_abs_start.restrict(dims)
        
        res.__curr_abs_end = self.__curr_abs_end.restrict(dims)
        
        return res
    
    
    def __get_bounds(self, did):
        
        start = self.__curr_abs_start.get_obj(did)
        
        end = self.__curr_abs_end.get_obj(did)
        
        return [start, end]
    
    
    def set_curr(self, incr, sel_bnd, did):
        
        dst_dims = (did, did) 
        
        base = deep(self.__base_points)
        
        lod = deep(self.__lod)
        
        if sel_bnd is None:
        
                sel_bnd = self.__get_bounds(did)
                
        if incr != 0:
            
                try:
                    
                    self.__gen.check(incr, sel_bnd, lod, did)
                
                except InvalidRequest:
                    
                    pass
                
                pj_bnd = self.__gen.pjbnd(incr, sel_bnd, lod, base, did)
                
                labels = self.__gen.label_trf(incr, None, lod, base, did)
                
                points = self.__gen.point_trf(incr, None, lod, base, did)
                
                self.__base_points = pt.PointSet(labels, points)
                
                base = pt.PointSet(labels, points)
                
                new_bounds = pj_bnd
                
                new_lod = self.__lod.get_obj( dst_dims ) + incr
                
                self.__lod.replace(mst.MdimStruct( [new_lod], [ (0,0) ] ), dst_dims, 0)
                
        else:
            
            new_bounds = sel_bnd
        
        try:
            
            self.__gen.check(0, new_bounds, lod, did)
            
        except InvalidRequest:                        
            
            pass
        
        labels = self.__gen.label_trf(None, new_bounds, lod, base, did)
        
        points = self.__gen.point_trf(None, new_bounds, lod, base, did)
        
        #TODO: clean up workaround from here ...
        curr_labels = self.__curr_points.get_labels()
        
        labels = labels.restrict(dst_dims)
        
        curr_labels.replace(labels, dst_dims)
        
        self.__curr_points = pt.PointSet(curr_labels, points)
        #... to here.
        
        #TODO:
        #Limits - start
        #dims will be used for LODs and imits
        dimspan = range(0, self.__base_points.dims() )
        all_dims = [ (did, did) for did in dimspan ]
        
        it_curr = self.__curr_points.get_labels().__iter__('ND')
        
        start_objs = deep( it_curr.nd_index() )
        
        self.__curr_abs_start = mst.MdimStruct(start_objs, all_dims)
            
        it_curr.set_last()
            
        end_objs = deep( it_curr.nd_index() )
        
        self.__curr_abs_end = mst.MdimStruct(end_objs, all_dims)
        
        return
    
    
    def point(self, inds):
        #TODO: check 'inds' OK?
        return self.__curr_points.point(inds)
    
    
    def get_restr_points(self, dims):
        
        res = None
        
        res = self.__curr_points.restrict(dims)
        
        return res
    
    
    def get_lims(self):
         
        res = None
         
        dimspan = range(0, self.__base_points.dims() )
        all_dims = [ (did, did) for did in dimspan ]
        
        lims = []
        for dims in all_dims:
            
            rgrid = self.restrict(dims)
        
            curr_end = rgrid.get_end()
        
            next_end = curr_end.get_obj( (0,0) ) + 1
        
            rgrid.set_curr(0, next_end, 0)
        
            lims.append( [0.0, rgrid.point( [next_end] ) ] )
         
        for lim in lims:
        
            limstr = mst.MdimStruct( [lim], [ (0,0) ] )
            
            if res is None:
                
                res = limstr
                
            else:
                
                res = res.prod(limstr)
        
        return res

    
    def last_did(self):
        
        res = self.__base_points.last_did()
        
        return res
    
    
    def get_start(self):
        
        return self.__curr_abs_start
    
    
    def get_end(self):
        
        return self.__curr_abs_end


    def get_curr_cards(self):
        
        return  self.__curr_points.get_labels().get_cards()
    
    
    def get_card(self, did):
        
        return self.__base_points.get_labels().part_card(did)
    
    
    def get_curr_card(self, did):
        
        return self.__curr_points.get_labels().part_card(did)