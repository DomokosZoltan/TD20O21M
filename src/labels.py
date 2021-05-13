# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:34:08 2021

@author: Domokos Zoltán
"""

import numpy as np
import copy as cp

import mstruct as mst


def default_labels(cards):
    
    res = None
    
    if not cards is None:
    
        labelfcs = [lambda x: [ x[0] ] if x[0] < card[1]
                    else None for card in cards]
        
        fcsdims = [card[0] for card in cards]
    
        res = mst.MdimStruct(labelfcs, fcsdims)

    else:
        #TODO: raise exception
        pass
      
    return res

#TODO: rename to LabelMap
class LabelSet():
    
    
    def __init__(self, cards = None, labels = None):
        
        if cards is None:
            
            self.__cards = mst.MdimStruct( [4], [ (0,0) ])
        
        else:
        
            self.__cards = cards
        
        if labels is None:
            
            self.__labels = default_labels(self.__cards)
            
            
        else:
            
            self.__labels = labels
        
        return


    def __iter__(self, order = 'F', direction = 'fwd', did = None):
          
        res = None
        
        if order != 'ND':
        
            res = LabelSetIterator(self, order, direction)
            
        else:
        
            if did is None:
                
                did = 0
            
            res = LabelSetNDIterator(self, did, direction)
    
        return res
    
    
    def prod(self, rhs):
                
        res = None
        
        labels = self.__labels.prod(rhs.__labels)
        
        cards = self.__cards.prod(rhs.__cards)
        
        res = LabelSet(cards, labels)
        
        return res
    
    
    def replace(self, rhs, dst_dims, src_offs = 0):
        
        self.__cards.replace(rhs.__cards, dst_dims, src_offs)
        
        self.__labels.replace(rhs.__labels, dst_dims, src_offs)
        
        return
    
    
    def restrict(self, dims):
        
        res = None
        
        labels = self.__labels.restrict(dims)
        
        cards = self.__cards.restrict(dims)
        
        res = LabelSet(cards, labels)
        
        return res 

    
    def dims(self):
        
        res = None
        
        res = self.__cards.last_did() + 1
        
        return res
    
   
    def full_card(self):
        
        res = None
    
        for card in self.__cards:
            
            if res is None:
    
                res = card[1]
                
            else:
                
                res *= card[1]
        
        return res
        
    
    #For convenience.
    def get_cards(self):
        
        res = None
        
        res = cp.deepcopy(self.__cards)
        
        return res
    
    def part_card(self, did):
            
        res = None
        
        res = self.select_cards( (did, did) )[0]
        
        return res
        
       
    def select_cards(self, dims):
        
        res = None
        
        restr_cards = self.__cards.restrict(dims)
        
        res = []
        for card in restr_cards:
            
            res += [ card[1] ]
        
        return res
        
        
    #TODO: generalize
    def is_in(self, inds):
        
        res = True
        
        did = 0
        for card in self.__cards:
            
            if card is None:
                
                res = False
                
            elif inds[did] >= card[1]:
                
                res = False
            
            if res == False:
                
                break
            
            did += 1
        
        return res


    def label(self, inds):
        
        res = None
        
        if self.is_in(inds):

            res = []
            
            for lb in self.__labels:
                
                sl = mst.dims_to_sl( lb[0] )
                    
                fc = lb[1]
                    
                res += fc( inds[sl] )
                
        return res


#TODO: merge with ND iterator than remove this
class LabelSetIterator:
    
    def __init__(self, labels, order, direction):
        
        self.__labels = labels
        
        self.__order = order
        
        self.__direction = direction
        
        self.__flind = 0

        return
    
    
    def __next__(self):
        
        res = None
        
        if self.__labels.full_card() > self.__flind:
            
            res = self.nd_index()
            
            self.__flind += 1
        
        else:
        
            raise(StopIteration)
        
        return res


    def set_direction(self, direction):
        self.__direction = direction
        return


    #TODO: can be removed later
    def __inds(self, order, k, num):
        dims = self.__labels.dims()
        if 'C' == order:
            m = dims-1-k
        else:
            m = k
        card = self.__labels.part_card(m)
        if k < (dims-1):
            #recursive cases
            succnum = np.floor_divide(num, card)
            succls = self.__inds(order, (k+1), succnum)
            if 'C' == order:
                resls = succls + [np.mod(num, card)]
            else:
                resls = [np.mod(num, card)] + succls
            #return resls
        else:
            #base case
            resls = [num]
        
        return resls
    
    
    def reset(self):
        self.__flind = 0
        return
    
    
    def set_to(self, ind):
        
        res = None
        
        self.reset()

        next_ind = 0
        while next_ind < ind:
        
            try:
                
                self.__next__()                
            
                next_ind = self.flat_index()
                
            except StopIteration:
               
                next_ind = None
                
                break
                
        res = next_ind
            
        return res
    
    
    #TODO: make this more efficient. Good for now.
    def set_last(self):
        
        res = None
        
        self.reset()
        
        last = self.__labels.full_card() - 1
        
        while last  > self.__flind:
        
            res = self.__next__()
        
        return res
    
    
    def flat_index(self):
        
        return cp.deepcopy(self.__flind)
    
    
    def nd_index(self):
        
        res = None
        
        inds = cp.deepcopy(self.__inds(self.__order, 0, self.__flind) )
        
        res = self.__labels.label(inds)
        
        return res


class LabelSetNDIterator:
    
    
    def __init__(self, labels, did, direction):
        
        self.__did = did
        
        self.__direction = 'fwd'
        
        self.__labels = labels
        
        self.__inds = [0 for did in range(0, self.__labels.dims() ) ]
        
        #TODO: remove this workaround
        self.__stopped = False
        
        return
    
    
    def __next__(self):
        
        res = None
        
        if 'fwd' == self.__direction: 
            
            lim = self.__labels.part_card(self.__did)
            
            if (not self.__stopped) and (self.__inds[self.__did] < lim):
                
                mapped_inds = self.__labels.label(self.__inds)
                
                self.__inds[self.__did] += 1
            
            else:
            
                self.__stopped = True
                
                raise(StopIteration)
        
        elif 'bwd':
            
            if (not self.__stopped) and (self.__inds[self.__did] >= 0):
            
                mapped_inds = self.__labels.label(self.__inds)
                
                self.__inds[self.__did] -= 1
            
            else:
            
                self.__stopped = True
                
                raise(StopIteration)           
        
        res = cp.deepcopy(mapped_inds)
        
        return res


    def set_did(self, did):
        self.__did = did
        return
    
    
    def set_direction(self, direction):
        self.__direction = direction
        return
      
        
    def reset(self, did = None):
        
        if did is None:
            
            self.__inds = [0 for did in range(0, self.__labels.dims() ) ]
        
        else:
        
            self.__inds[did] = 0
        
        return
    
    
    def set_to(self, inds):
        
        res_inds = None
        
        self.__inds = cp.deepcopy(inds)
        
        return res_inds
    
    
    #TODO: make this more efficient. Good for now.
    def set_last(self):
        
        self.reset()
        
        for did in range(0, len(self.__inds) ):
        
            self.__inds[did] = self.__labels.part_card(did) - 1 
        
        return
    
    
    def flat_index(self, order = 'F'):
        #TODO
        return
    
    
    def nd_index(self):
        
        return cp.deepcopy(self.__inds)
    
    