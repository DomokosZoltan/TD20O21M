# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:45:08 2021

@author: Domokos Zoltán
"""

import numpy as np
import copy as cp

import cycdiv as cd
import polardiv as pl

import grid as gd
import mstruct as mst
import generator as gn

import unittest


class TestGrid(unittest.TestCase):

    def setUp(self):
        
        cr = 8
        ca = 8
        
        lim_r = 31
        lim_a = 7
        
        lbpar_r = pl.LbFixParams( [cr], [lim_r] )
        lbpar_a = pl.LbFixParams( [ca], [lim_a] )
        
        shf_r = 0.0
        shf_a = 0.0
        
        band_r = 8.0
        band_a = 2.0 * np.pi
        
        dr = band_r / cr
        da = band_a / ca
        
        diff = [dr, da]
        shf = [shf_r, shf_a]
        band = [band_r, band_a]
        
        ptpar = pl.PtFixParams(diff, shf, band)
        
        pol_lb_r = pl.PolarLbGen(lbpar_r)
        pol_lb_a = pl.PolarLbGen(lbpar_a)
        
        pol_pt = pl.PolarPtGen(ptpar)
        
        dim_lb_r = (0, 0)
        dim_lb_a = (1, 1)
        
        obj_lb = [pol_lb_r, pol_lb_a]
        dims_lb = [dim_lb_r, dim_lb_a]
        
        mst_lb = mst.MdimStruct(obj_lb, dims_lb) 
        
        dims_pt = [ (0, 1) ]
        obj_pt = [pol_pt]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        self.gdPol = gd.Grid(gen)

        #linear objects        
        ct = 4

        lim_t = 31
        
        lbpar_t = cd.LbFixParams( [ct], [lim_t] )

        shf_t = 0.0
        
        band_t = 4.0
        
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
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        self.gdA = gd.Grid(gen)
        
        self.gdB = gd.Grid(gen)
        
        return
    
    def test_polar_grid_bounds(self): 
        
        gdPol = cp.deepcopy(self.gdPol)
        
        #VALUES BEFORE BOUND UPDATE
        self.assertAlmostEqual(gdPol.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdPol.point( [0, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 0] )[0], 1.0)
        self.assertAlmostEqual(gdPol.point( [1, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 1] )[0], 1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 1] )[1], 1.0/np.sqrt(2.0) )
        
        self.assertAlmostEqual(gdPol.point( [1, 2] )[0], 0.0)
        self.assertAlmostEqual(gdPol.point( [1, 2] )[1], 1.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 3] )[0], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 3] )[1], 1.0/np.sqrt(2.0) )
        
        self.assertAlmostEqual(gdPol.point( [1, 7] )[0], 1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 7] )[1], -1.0/np.sqrt(2.0) )
        
        gdPol.set_curr(0, [10, 2], 1)
        
        #VALUES AFTER BOUND UPDATE
        self.assertAlmostEqual(gdPol.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdPol.point( [0, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 0] )[0], 0.0)
        self.assertAlmostEqual(gdPol.point( [1, 0] )[1], 1.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 1] )[0], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 1] )[1], 1.0/np.sqrt(2.0) )
        
        self.assertAlmostEqual(gdPol.point( [1, 2] )[0], -1.0)
        self.assertAlmostEqual(gdPol.point( [1, 2] )[1],  0.0)
        
        self.assertAlmostEqual(gdPol.point( [1, 3] )[0], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 3] )[1], -1.0/np.sqrt(2.0) )
        
        self.assertAlmostEqual(gdPol.point( [1, 7] )[0], 1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdPol.point( [1, 7] )[1], 1.0/np.sqrt(2.0) )
        
        self.assertAlmostEqual(gdPol.point( [1, 8] )[0], 0.0)
        self.assertAlmostEqual(gdPol.point( [1, 8] )[1], 1.0)
    
        return
    
    
#    def test_polar_grid_depth(self): 
#        
#        #TODO
#        
#        #    gdC.set_curr(0, [0,3], 1)
#        #    gdC.set_curr(1, [0,3], 1)
#        #    
#        #    print('DEPTH UPDATE, VALUES AFTER')
#        #    print(gdC.point([0,0]) )
#        #    print(gdC.point([1,0]) )
#        #    print(gdC.point([0,2]) )
#        #    print(gdC.point([1,1]) )
#        #    print(gdC.point([2,1]) )
#        #    
#        #    gdC.set_curr(-1, [0,4], 1)
#        #    
#        #    print('DEPTH UPDATE, REVERT IT')
#        #    print(gdC.point([0,0]) )
#        #    print(gdC.point([1,0]) )
#        #    print(gdC.point([0,2]) )
#        #    print(gdC.point([1,1]) )
#        #    print(gdC.point([2,1]) )
#        #    
#        #    gdC.set_curr(0, [0,3], 1)
#        #    
#        #    print('DEPTH UPDATE, REVERT IT')
#        #    print(gdC.point([0,0]) )
#        #    print(gdC.point([1,0]) )
#        #    print(gdC.point([0,3]) )
#        #    print(gdC.point([1,1]) )
#        #    print(gdC.point([2,1]) )
#    
#        return

    def test_cyc_grid_prod_and_bounds(self):
        
        gdA = cp.deepcopy(self.gdA)
        gdB = cp.deepcopy(self.gdB)
            
        self.assertAlmostEqual(gdA.point( [0] )[0], 0.0)
        self.assertAlmostEqual(gdB.point( [0] )[0], 0.0)
        
        self.assertAlmostEqual(gdA.point( [1] )[0], 1.0)
        self.assertAlmostEqual(gdB.point( [1] )[0], 1.0)
        
        self.assertAlmostEqual(gdA.point( [2] )[0], 2.0)
        self.assertAlmostEqual(gdB.point( [2] )[0], 2.0)
        
        self.assertAlmostEqual(gdA.point( [3] )[0], 3.0)
        self.assertAlmostEqual(gdB.point( [3] )[0], 3.0)
        
        gdC = gdA.prod(gdB)
        
        #VALUES BEFORE BOUND UPDATE
        self.assertAlmostEqual(gdC.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [1, 0] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [0, 1] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 1] )[1], 1.0)
        
        self.assertAlmostEqual(gdC.point( [1, 1] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 1] )[1], 1.0)

        self.assertAlmostEqual(gdC.point( [2, 1] )[0], 2.0)
        self.assertAlmostEqual(gdC.point( [2, 1] )[1], 1.0)
        
        #the cardinalities
        self.assertAlmostEqual(gdC.get_card(0), 4)
        self.assertAlmostEqual(gdC.get_card(1), 4)
            
            
        gdC.set_curr(0, [10,2], 1)
        
        #VALUES AFTER BOUND UPDATE
        self.assertAlmostEqual(gdC.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 0] )[1], 10.0)
        
        self.assertAlmostEqual(gdC.point( [1, 0] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 0] )[1], 10.0)
        
        self.assertAlmostEqual(gdC.point( [0, 6] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 6] )[1], -16.0)
        
        self.assertAlmostEqual(gdC.point( [1, 1] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 1] )[1], 11.0)

        self.assertAlmostEqual(gdC.point( [2, 1] )[0], 2.0)
        self.assertAlmostEqual(gdC.point( [2, 1] )[1], 11.0)
        
        gdC.set_curr(0, [0,3], 1)
        gdC.set_curr(1, [0,3], 1)
        
        #VALUES AFTER DEPTH UPDATE
        self.assertAlmostEqual(gdC.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [1, 0] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [0, 4] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 4] )[1], 2.0)
        
        self.assertAlmostEqual(gdC.point( [1, 1] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 1] )[1], 0.5)

        self.assertAlmostEqual(gdC.point( [2, 5] )[0], 2.0)
        self.assertAlmostEqual(gdC.point( [2, 5] )[1], 2.5)
        
        gdC.set_curr(-1, [0,15], 1)
        
        #DEPTH UPDATE, REVERT IT
        self.assertAlmostEqual(gdC.point( [0, 0] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [1, 0] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 0] )[1], 0.0)
        
        self.assertAlmostEqual(gdC.point( [0, 7] )[0], 0.0)
        self.assertAlmostEqual(gdC.point( [0, 7] )[1], 7.0)
        
        self.assertAlmostEqual(gdC.point( [1, 1] )[0], 1.0)
        self.assertAlmostEqual(gdC.point( [1, 1] )[1], 1.0)
        
        return


    def test_replace_dims(self):
        
        gdPol = cp.deepcopy(self.gdPol)
        gdB = cp.deepcopy(self.gdB)
        
        #assemble the destination grid
        gdDst = cp.deepcopy(self.gdA)
        
        gdDst = gdDst.prod(gdPol)
        
        gdDst = gdDst.prod(gdB)
        
        gdDst = gdDst.prod(gdPol)
        
        #check cards after product
        self.assertAlmostEqual(gdDst.get_card(0), 4)
        self.assertAlmostEqual(gdDst.get_card(1), 8)
        self.assertAlmostEqual(gdDst.get_card(2), 8)
        self.assertAlmostEqual(gdDst.get_card(3), 4)
        self.assertAlmostEqual(gdDst.get_card(4), 8)
        self.assertAlmostEqual(gdDst.get_card(5), 8)
        
        #get a point
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[0], 2.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[1], 1.0 )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[2], 0.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[3], 1.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[4], 1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[5], 1.0/np.sqrt(2.0) )
        
        #now assemble a source grid that has a dimansionality that matches
        #the dst dims when [ (1,2), (3,3) ] is replaced
        gdSrc = cp.deepcopy(self.gdPol)
        
        gdSrc = gdSrc.prod(gdB)
        
        #update source grid to see some difference
        gdSrc.set_curr(0, [2, 10], 1)
        gdSrc.set_curr(1, [0, 3], 2)
        
        #check a point
        self.assertAlmostEqual(gdSrc.point( [1, 0, 1] )[0], 0.0)
        
        self.assertAlmostEqual(gdSrc.point( [1, 0, 1] )[1], 1.0)
        self.assertAlmostEqual(gdSrc.point( [1, 0, 1] )[2], 0.5)
        
        #replace one 2d and a 1d subgrid
        gdDst.replace(gdSrc, (1, 3), 0)
        
        #get a point from the modified destination grid
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[0], 2.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[1], 0.0 )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[2], 1.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[3], 0.5 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[4], 1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[5], 1.0/np.sqrt(2.0) )
        
        #replace just a 1d subgrid with 1d
        gdSrc.set_curr(-1, [0, 7], 2)
        
        gdDst.replace(gdSrc, (3, 3), 2)
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[3], 1.0 )
        
        #replace just a 2d subgrid with a 2d
        
        gdDst.replace(gdSrc, (4, 5), 0)
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[0], 2.0 )
    
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[1], 0.0 )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[2], 1.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[3], 1.0 )
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[4], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdDst.point( [2, 1, 0, 1, 1, 1] )[5], 1.0/np.sqrt(2.0) )
        
        #replace 2d subgrid with two 1d subgrids
        gdSrc = gdSrc.prod(gdB)
        
        gdSrc.set_curr(1, [0, 3], 3)
        
        gdDst.replace(gdSrc, (1, 2), 2)
        
        #Note that we expect reduced dimensionality!
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[0], 2.0 )        
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[1], 2.0 )        
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[2], 0.5 )        
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[3], 1.0 )        
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[4], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdDst.point( [2, 2, 1, 1, 1, 1] )[5], 1.0/np.sqrt(2.0) )
        
        #replace two 1d subgrids with one 2d subgrid
        gdSrc.set_curr(0, [2, 10], 1)
        
        gdDst.replace(gdSrc, (1, 2), 0)
        
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[0], 2.0 )        
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[1], -1.0/np.sqrt(2.0) )        
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[2], 1.0/np.sqrt(2.0) )        
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[3], 1.0 )        
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[4], -1.0/np.sqrt(2.0) )
        self.assertAlmostEqual(gdDst.point( [2, 1, 1, 1, 1, 1] )[5], 1.0/np.sqrt(2.0) )
        
        return


if __name__ == '__main__':
    unittest.main()