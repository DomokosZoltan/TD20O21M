# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:16:20 2021

@author: Domokos Zoltán
"""

import mstruct as mst

import labels as lb
import points as pt

import unittest

class TestPointSet(unittest.TestCase):

    def setUp(self):
        
        cardsA = mst.MdimStruct( [2], [ (0, 0) ] )        
        cardsB = mst.MdimStruct( [3], [ (0, 0) ] )
        
        self.lsA = lb.LabelSet(cardsA, None)
        self.lsB = lb.LabelSet(cardsB, None)
        
        return
    
    
    def testPointSet(self):
        
        def pointfcA(inds):
            res=[]
            for x in inds:
                res.append(10 * x)
            return res
        
        fcs = [pointfcA]
        dims = [(0,0)]
        
        ptA = pt.PointSet(self.lsA, mst.MdimStruct(fcs, dims) )
               
        self.assertAlmostEqual( 0.0, ptA.point( [0] )[0] )
        self.assertAlmostEqual( 10.0, ptA.point( [1] )[0] )
        
        def pointfcB(inds):
            res=[]
            for x in inds:
                res.append(0.1*x)
            return res
        
        fcs = [pointfcB]
        dims = [(0,0)]
        
        ptB = pt.PointSet(self.lsB, mst.MdimStruct(fcs, dims) )
        
        self.assertAlmostEqual( ptB.point( [0] )[0], 0.0 )
        self.assertAlmostEqual( ptB.point( [1] )[0], 0.1 )
        self.assertAlmostEqual( ptB.point( [2] )[0], 0.2 )
        
        ptC = ptA.prod(ptB)
        
        self.assertAlmostEqual( ptC.point( [1,1] )[0], 10.0 )
        self.assertAlmostEqual( ptC.point( [1,1] )[1], 0.1 )
        
        #simple replace
        ptC.replace(ptB, (0,0), 0)
        
        self.assertAlmostEqual( ptC.point( [1,1] )[0], 0.1 )
        self.assertAlmostEqual( ptC.point( [1,1] )[1], 0.1 )
        
        #simple restrict
        ptC = ptC.prod(ptA)
        
        ptD = ptC.restrict( (1, 2) )
        
        self.assertAlmostEqual( ptD.point( [1,1] )[0], 0.1 )
        self.assertAlmostEqual( ptD.point( [1,1] )[1], 10.0 )
        
        return


if __name__ == '__main__':
    unittest.main()