# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 19:53:15 2021

@author: Domokos Zoltán
"""

import unittest

from copy import deepcopy as deep

import mstruct as mst
import generator as gn
import grid as gd
import cycdiv as cd
import domain as dm
import codomain as codm
import sglspace as sp
import sglutils as sgl

import view as vw

class TestView(unittest.TestCase):

    def setUp(self):

        samp_shf = 0.00
        
        #linear SAMPLE objects        
        ct = 4

        lim_t = 31
        
        lbpar_t = cd.LbFixParams( [ct], [lim_t] )

        shf_t = samp_shf
        
        band_t = 1.0
        
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
        shf = [shf_t - samp_shf]
        band = [band_t]
        
        ptpar = cd.PtFixParams(diff, shf, band)
        
        pt_t = cd.PtGen(ptpar)
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        genB = gn.Generator(mst_lb, mst_pt) 
        
        self.gdB = gd.Grid(genB)
        
        bandstr = mst.MdimStruct([ [0.0, 1.0] ], [(0,0)])
        
        dmA = dm.Domain(self.gdS, self.gdB, None, None, True, bandstr)
        
        self.dmP2d = dmA.prod(dmA)
        
        self.dmP3d = self.dmP2d.prod(dmA)       
        
        #set up codomains
        T = dmA.bandlims.get_obj((0,0))[1]
        n = dmA.ns.get_obj((0,0)) 
        dt = T/n
        shf = 2
        phi = 1.0 + 0.0j
        ampl = 1.0
        symb = 't'
    
        fcstr = sgl.cyc_stepfc(dt, T, shf, phi, ampl, symb)
        dims = (0,0)
    
        codmstr = mst.MdimStruct( [fcstr], [dims] )
        
        codmA = codm.Codomain(codmstr)
        
        shf = 1
        
        fcstr = sgl.cyc_stepfc(dt, T, shf, phi, ampl, symb)
        dims = (0,0)
        
        codmstr = mst.MdimStruct( [fcstr], [dims] )
        
        codmB = codm.Codomain(codmstr)
        
        self.codmP2d = codmA.prod(codmB)
        
        return


    def testHeatVec(self):
        
        codmP2d = deep(self.codmP2d)
        domain = deep(self.dmP2d)
        
        sigvec = sp.cdmvec(codmP2d.fcval_at, domain)[1]
        
        vw.plot_heatvec(domain, sigvec, md = 'cartes')
        
        return
    

if __name__ == '__main__':
    unittest.main()