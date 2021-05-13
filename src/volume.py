# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:54:22 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep

import mstruct as mst 


def cube(part_grid, inds):

    points = []
        
    for did in range(0, len(inds) ):
        #TODO: it could be made more efficient, if the first vertex point 
        #in the bounding volume was accessed only once
        start = part_grid.point(inds)            
    
        tmp_inds = deep(inds)
        
        start = part_grid.point(tmp_inds)
        
        tmp_inds[did] += 1
        
        end = part_grid.point(tmp_inds) 
        
        points.append( tuple( (start, end) ) )
        
    return points    


def sector(part_grid, inds):
    
    points = []
        
    tmp_inds = deep(inds)
    
    for did in range(0, len(inds) ):
        #TODO: it could be made more efficient, if the first vertex point 
        #in the bounding volume was accessed only once
        start = part_grid.point(tmp_inds)
        
        tmp_inds[did] += 1
        
        end = part_grid.point(tmp_inds)
        
        points.append( tuple( (start, end) ) )
        
    return points


class Volume():
    
    
    def __init__(self, volmstr = None): 
    
        if volmstr is None:
            
            volfcs = [cube]
        
            dims = [ (0,0) ]
            
        else:
            
            volfcs = deep(volmstr._MdimStruct__objs)
        
            dims = deep(volmstr._MdimStruct__dims)
            
        self.__volstr = mst.MdimStruct(volfcs, dims)
        
        return
    
    
    def volume_at(self, bounds, inds):
        
        res = []
        
        for obj in self.__volstr:
            
            dims = obj[0]
            
            fc = obj[1]
            
            rgrid = bounds.restrict(dims)
    
            sl = mst.dims_to_sl(dims)
            
            res += fc(rgrid, inds[sl] )
            
        return res


    def prod(self, rhs):
        
        res = None
        
        res = Volume(self.__volstr.prod(rhs.__volstr) )
        
        return res
            
    
    def replace(self, rhs, dims, src_offs = 0):
        
        self.__volstr.replace(rhs.__volstr, dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
        
        res = Volume(self.__volstr.restrict(dims) )
        
        return res
    
