# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 19:58:58 2021

@author: domok
"""

import domain as dm

import mstruct as mst
import generator as gn
import grid as gd
import cycdiv as cd

import unittest

class TestDomain(unittest.TestCase):

    def setUp(self):
        
        #linear SAMPLE objects        
        ct = 4

        lim_t = 31
        
        lbpar_t = cd.LbFixParams( [ct], [lim_t] )

        shf_t = 5.0
        
        band_t = 40.0
        
        dt = band_t / ct
        
        diff = [dt]
        shf = [shf_t]
        band = [band_t]
        
        ptpar = cd.PtFixParams(diff, shf, band)
        
        lb_t = cd.LbGen(lbpar_t)
        
        pt_t = cd.PtGen(ptpar)
        
        dim_lb = (0, 0)
        
        obj_lb = [lb_t]
        dims_lb = [dim_lb]
        
        mst_lb = mst.MdimStruct(obj_lb, dims_lb) 
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        genS = gn.Generator(mst_lb, mst_pt)
        
        self.gdS = gd.Grid(genS)
        
        #BOUNDS
        diff = [dt]
        shf = [shf_t - 5.0]
        band = [band_t]
        
        ptpar = cd.PtFixParams(diff, shf, band)
        
        pt_t = cd.PtGen(ptpar)
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        genB = gn.Generator(mst_lb, mst_pt) 
        
        self.gdB = gd.Grid(genB)
        
        return
    
    def testSamples_single_eqdst(self):
        
        bandstr = mst.MdimStruct( [ [0.0, 40.0] ], [ (0,0) ] )
        
        dmA = dm.Domain(self.gdS, self.gdB, None, None, True, bandstr)
        
        self.assertAlmostEqual(dmA.sample_at( [1] )[0], 15.0)
        
        self.assertAlmostEqual(dmA.measure_at( [0] ), 10.0)
        
        self.assertAlmostEqual(dmA.measure_at( [3] ), 10.0)
        
        return


    def testSamples_prod_eqdst(self):
        
        bandstr = mst.MdimStruct( [ [0.0, 40.0] ], [ (0,0) ] )
        
        dmA = dm.Domain(self.gdS, self.gdB, None, None, True, bandstr)
        
        dmP = dmA.prod(dmA)
        
        self.assertAlmostEqual(dmP.sample_at( [1, 1] )[0], 15.0)
        self.assertAlmostEqual(dmP.sample_at( [1, 1] )[1], 15.0)
        
        self.assertAlmostEqual(dmP.measure_at( [0, 1] ), 100.0)
        
        self.assertAlmostEqual(dmP.measure_at( [3, 3] ), 100.0)


if __name__ == '__main__':
    unittest.main()