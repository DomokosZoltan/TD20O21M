# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:34:24 2021

@author: Domokos Zoltán
"""

import numpy as np

import mstruct as mst
import labels as lb

import unittest


class TestLabelSet(unittest.TestCase):
    
    def setUp(self):
        
        cardsA = mst.MdimStruct( [2], [ (0, 0) ] )        
        cardsB = mst.MdimStruct( [3], [ (0, 0) ] )
        
        self.lsA = lb.LabelSet(cardsA, None)
        self.lsB = lb.LabelSet(cardsB, None)
        
        return

    def testLabelSetDims(self):
        
        self.assertAlmostEqual(self.lsA.dims(), 1)
        self.assertAlmostEqual(self.lsB.dims(), 1)
        
        lsC = self.lsA.prod(self.lsB)
        
        self.assertAlmostEqual(lsC.dims(), 2)
        
        return
    
    def testLabelSetCards(self):
        
        self.assertAlmostEqual(self.lsA.full_card(), 2)
        self.assertAlmostEqual(self.lsA.part_card(0), 2)
        
        self.assertAlmostEqual(self.lsB.full_card(), 3)
        self.assertAlmostEqual(self.lsB.part_card(0), 3)
        
        lsC = self.lsA.prod(self.lsB)
        
        self.assertAlmostEqual(lsC.full_card(), 6)
        self.assertAlmostEqual(lsC.part_card(0), 2)
        self.assertAlmostEqual(lsC.part_card(1), 3)
        
        return


    def testLabelSetFlatIterator(self):
        #set up product label sets
        lsC = self.lsA.prod(self.lsB) 
        lsD = lsC.prod(self.lsB.prod(self.lsA) )
        #get an iterator
        itls = lsD.__iter__('F')
    
        #set up a reference numpy array
        refarr = np.arange(36).reshape(2, 3, 3, 2)
        #get an nditerator
        ndit = np.nditer(refarr, flags=['multi_index'], order='F')
    
        for x in ndit:
            
            with self.subTest():
                
                refls = ndit.multi_index
                currls = itls.__next__()
                
                self.assertAlmostEqual( refls[0], currls[0] )
                self.assertAlmostEqual( refls[1], currls[1] )
                self.assertAlmostEqual( refls[2], currls[2] )
                self.assertAlmostEqual( refls[3], currls[3] )
                
        itls = lsD.__iter__('C')
        
        ndit = np.nditer(refarr, flags=['multi_index'], order='C')
    
        for x in ndit:
            
            with self.subTest():
                
                refls = ndit.multi_index
                currls = itls.__next__()
                
                self.assertAlmostEqual( refls[0], currls[0] )
                self.assertAlmostEqual( refls[1], currls[1] )
                self.assertAlmostEqual( refls[2], currls[2] )
                self.assertAlmostEqual( refls[3], currls[3] )
                
        return


    def testLabelSetNDIterator(self):
        #set up product label sets
        lsC = self.lsA.prod(self.lsB) 
        #get an iterator
        itls = lsC.__iter__('ND')
    
        #set up a reference numpy array
        refarr = np.arange(6).reshape(2, 3)
        #get an nditerator
        ndit = np.nditer(refarr, flags=['multi_index'], order='F')
    
        ind = 0
        for x in ndit:
            
            with self.subTest():
                
                refls = ndit.multi_index                
                
                currls = itls.__next__()
                
                self.assertAlmostEqual( refls[0], currls[0] )
                self.assertAlmostEqual( refls[1], currls[1] )
                
                ind += 1
                
                if 0 == ind % 2:
                    itls.reset(0)
                    itls.set_did(1)
                    itls.__next__()
                    itls.set_did(0)
                    
        return


if __name__ == '__main__':
    unittest.main()
