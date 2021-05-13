# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 14:41:51 2021

@author: Domokos Zoltán
"""

from copy import deepcopy as deep


def dims_to_sl(dims):
    
    sl = None
    
    if dims[0] != dims[1]:
        
        sl = slice(dims[0], dims[1] + 1)
        
    else:
        
        sl = slice(dims[0], dims[0] + 1)
        
    return sl
            

class MdimStruct:
    
    
    def __init__(self, objs, dims):
        
        self.__objs = deep(objs)

        self.__dims = deep(dims) 
        
        return
     
        
    def __iter__(self, direction = 'fwd'):  
    
         return MdimStructIterator(self, direction)   
    
    
    def last_did(self):
        
        res = None
        
        try:
        
            res = self.__dims[-1][1]
        
        except IndexError:
            #TODO: handle
            pass
        
        return res
    
    
    def prod(self, rhs):
        
        res = None
        
        lhs_dims = self.last_did() + 1
        
        rhs_dims = map(lambda elm : tuple ( (elm[0] + lhs_dims,\
                                             elm[1] + lhs_dims) ),\
                          rhs.__dims )
        
        dims = self.__dims + list(rhs_dims)
        
        objs = self.__objs + rhs.__objs
        
        res = MdimStruct(objs, dims)
        
        return res
    
    
    def replace(self, rhs, dst_dims, src_offs = 0):
        
        new_objs = []
        new_dims = []
        
        #Before selected region
        lim = dst_dims[0]
            
        dst_did = 0
        while (dst_did < lim):
                
            curr_dst_dims = self.get_dims(dst_did)
            dst_fcid = self.get_objid(curr_dst_dims)
                
            new_dims.append(curr_dst_dims)
            new_objs.append( deep(self.__objs[dst_fcid] ) )
                
            dst_did += ( (curr_dst_dims[1] + 1) - curr_dst_dims[0] )
              
        lim = (dst_dims[1] + 1) - dst_dims[0] + src_offs
              
        #Selected region
        src_did = src_offs
        while (src_did < lim):
                
            curr_src_dims = rhs.get_dims(src_did)
            src_fcid = rhs.get_objid(curr_src_dims)
                
            curr_new_dims = (curr_src_dims[0] - src_offs + dst_dims[0],
                             curr_src_dims[1] - src_offs + dst_dims[0] )
                
            new_dims.append(curr_new_dims)
            new_objs.append( deep(rhs.__objs[src_fcid] ) )
                
            src_did += ( (curr_src_dims[1] + 1) - curr_src_dims[0] )
              
        lim = self.last_did() + 1
             
        #After selected region 
        dst_did = dst_dims[1] + 1
        while (dst_did < lim):
                
            curr_dst_dims = self.get_dims(dst_did)
            dst_fcid = self.get_objid(curr_dst_dims)
                
            new_dims.append(curr_dst_dims)
            new_objs.append( deep(self.__objs[dst_fcid] ) )
                
            dst_did += ( (curr_dst_dims[1] + 1) - curr_dst_dims[0] )
            
        #set new lists
        self.__dims = deep(new_dims)
        self.__objs = deep(new_objs)
        
        return
    
    
    def restrict(self, new_dims):
        
        tmp_objs = []        
        tmp_dims = []
        
        did = new_dims[0]
        while (did < (new_dims[1] + 1) ):
        
            dims_old = self.get_dims(did)
            
            ind = self.get_objid(dims_old)
            
            tmp_objs.append( deep(self.__objs[ind] ) )
                   
            tmp_dims.append( (dims_old[0] - new_dims[0], 
                              dims_old[1] - new_dims[0] ) )
    
            did += ( (dims_old[1] + 1) - dims_old[0] )
        
        res = MdimStruct(tmp_objs, tmp_dims)
        
        return res
    
    
    def get_obj(self, dims):
        objid = self.get_objid(dims)
        res = deep(self.__objs[objid] )
        return res


    def get_objid(self, dims):
        res = None
        for elm in self.__dims:
            if elm == dims:
                res = deep(self.__dims.index(elm) )
        return res
    

    def get_dims(self, did):
        res = None
        for elm in self.__dims:
            
            if ( (did == elm[0] ) and (did == elm[1] ) ) or\
               ( (elm[0]<=did ) and (did <= elm[1] ) ):
        
                res = deep(elm)
                break
        
        return res    


class MdimStructIterator:
    
    def __init__(self, mdimstruct, direction):
        
        self.__direction = direction
        
        self.__mdimstruct = mdimstruct
        
        self.__index = 0
        
        return
    
    def __next__(self):
        
        res = None
        
        if 'fwd' == self.__direction: 
            
            if self.__index < len(self.__mdimstruct._MdimStruct__dims):
                
                dims = self.__mdimstruct._MdimStruct__dims[self.__index]
                
                res = deep( (dims, self.__mdimstruct.get_obj(dims) ) )
        
                self.__index += 1
                
            else:
                
                raise StopIteration
        
        elif 'bwd':
            #TODO
            pass
        
        return res
    
    
    def at(self):        
        return self.__index