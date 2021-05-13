# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:15:58 2021

@author: Domokos Zoltán
"""

from copy import deepcopy

import mstruct as mst
import labels as lb


def default_points():

    res = None
    
    def pointfc(inds):
        res=[]
        
        for x in inds:
            res.append(1.0 * x)
        return res

    res = mst.MdimStruct( [pointfc], [(0, 0)] )

    return res 

class PointSet():
    
    
    def __init__(self, labels = None, points = None):
        
        if labels is None:
            
            self.__labels = lb.LabelSet()
            
        else:
        
            self.__labels = labels
        
        if points is None:
            
            self.__points = default_points()
        
        else:
            
            self.__points = points

        return


    def last_did(self):
        
        res = None
        
        res = self.__points.last_did()
        
        return res

    
    def prod(self, rhs):
                
        res = None
        
        labels = self.__labels.prod(rhs.__labels)
        
        points = self.__points.prod(rhs.__points)
        
        res = PointSet(labels, points)
        
        return res
    
    
    def replace(self, rhs, dims, src_offs = 0):
        
        self.__points.replace(rhs.__points, dims, src_offs)
    
        self.__labels.replace(rhs.__labels, dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
        
        points = self.__points.restrict(dims)
        
        labels = self.__labels.restrict(dims)
        
        res = PointSet(labels, points)
        
        return res
    
    
    def point(self, inds, labels = None):
        
        res = None
        
        if labels is None:
        
            labels = self.__labels
            
        mapped_inds = labels.label(inds)
        
        if not mapped_inds is None:
        
            res = []
            
            for pt in self.__points:
                    
                dims = pt[0]
                fc = pt[1]
                
                sl = mst.dims_to_sl(dims)
                
                res += fc( mapped_inds[sl] )
        
        return res
    
    
    def part_point(self, dims, inds, restr_labels = None):
        
        res = None
    
        sl = mst.dims_to_sl(dims)
    
        res = self.restrict(dims).point( inds[sl], restr_labels )
    
        return res


    def get_labels(self):
        
        res = deepcopy(self.__labels)
        
        return res
    
    
    def get_points(self):
        
        res = deepcopy(self.__points)
        
        return res

    #TODO: remove these legacy getters as soon as they will not be needed
    def dims(self):
        return self.__points.last_did() + 1

        
    def full_card(self):
        return self.__labels.full_card()


    def part_card(self, dim):
        res_card = None
        try:
            res_card = self.__labels.part_card(dim)
        except IndexError:
            #TODO: handle it
            pass
        return res_card

    
    def select_cards(self, dims):
        res_card = None  
        try:
            res_card = self.__labels.select_cards(dims)        
        except IndexError:
            #TODO: handle
            pass
        return res_card
