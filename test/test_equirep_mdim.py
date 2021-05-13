# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 13:20:31 2021

@author: Domokos Zoltán
"""

import numpy as np
from copy  import deepcopy as deep

import mstruct as mst
import domain as dm
import sglspace as sp

import sglgen as sgg
import equirep as eq #TODO: Temporarily base transforms will be tested here. Review this later.
import view as vw

import operators as op

import unittest


class TestEquirepMdim(unittest.TestCase):

    def setUp(self):
        
        st = [0.1, 0.39, 0.47, 0.82]
        bt = [0.0, 0.27, 0.43, 0.7]
        
        self.nufdm_t = dm.nufactory( [ [0.0, 1.0] ], [4], st, bt, [0.0], 'cartes')
        
        self.dm_t = dm.factory( [ [0.0, 1.0] ], [4], [0.0], 'cartes')
        
        self.dm_f = dm.factory( [ [0.0, 4.0] ], [4], [0.0], 'cartes')
                    
        self.dm_ra = dm.factory([ [0, 1.0], [0, 2.0 * np.pi] ], [4, 8], [0.0, 0.0], 'polar')
        
        self.dm_ra_shf = dm.factory([ [0, 1.0], [0, 2.0 * np.pi] ], [4, 8], [0.001, 0.0], 'polar')
        
        self.dm_s = dm.factory( [ [0.0, 4.0] ], [4], [0.0], 'cartes')
        
        self.dm_p = dm.factory( [ [0.0, 4.0 / (2.0 * np.pi) ] ], [4], [0.0], 'cartes')
        
        self.symb_tt = mst.MdimStruct( ['t', 't'], [ (0,0), (1,1) ] )
        
        self.symb_ff = mst.MdimStruct( ['f', 'f'], [ (0,0), (1,1) ] )
        
        self.symb_ra = mst.MdimStruct( ['ra'], [ (0,1) ] )
        
        self.symb_sp = mst.MdimStruct( ['sp'], [ (0,1) ] )
        
        gendict = {'t' : {'t' : sgg.get_codm_gen('t', 't'),
                          'f' : sgg.get_codm_gen('f','t') },
    
                   'f' : {'t' : sgg.get_codm_gen('t', 'f'),
                          'f' : sgg.get_codm_gen('f', 'f') } }
        
        gens = mst.MdimStruct( [gendict], [ (0,0) ] )
        
        dmgendict = {'t' : {'t' : None,
                            'f' : eq.get_cartes_gen('t', 'f') },
    
                     'f' : {'t' : eq.get_cartes_gen('f', 't'),
                            'f' : None } }
    
        dmgens = mst.MdimStruct( [dmgendict], [ (0,0) ] )
        
        rlsdict = {'t' : {'t' : None,
                          'f' : eq.get_XX_rule('t', 'f'),
                          'n': eq.get_Xn_rule('t') },
                                            
                   'f' : {'t' : eq.get_XX_rule('f', 't'),
                          'f' : None, 
                          'n': eq.get_Xn_rule('f') },
                                            
                    'n' : {'t' : eq.get_nX_rule('t'),
                           'f' : eq.get_nX_rule('f'), 
                           'n': None } }
                   
        crrls = mst.MdimStruct( [rlsdict], [ (0,0) ] )
        
        symb_t = mst.MdimStruct( ['t'], [ (0,0) ] )
        
        self.crtab = eq.Table(self.dm_t, symb_t, symb_t, crrls, gens, dmgens)
        
        codmgens = mst.MdimStruct( [ {'ra' : {'ra' : sgg.get_sphero_gen('ra', 'ra'), 
                                              'sp' : sgg.get_sphero_gen('sp','ra') },
                    
                                    'sp' : {'ra' : sgg.get_sphero_gen('ra', 'sp'),
                                            'sp' : sgg.get_sphero_gen('sp', 'sp') } } ],
                                    [ (0,1) ] )
        
        dmgens = mst.MdimStruct( [ {'ra' : {'ra' : None,
                                            'sp' : eq.get_polar_gen('ra', 'sp') },
    
                                  'sp' : {'ra' : eq.get_polar_gen('sp', 'ra'),
                                          'sp' : None } } ],
                                    [ (0,1) ] )
    
        rules0 = {'r' : {'s' : eq.get_XX_rule('s', 'r', 0, ['s','p'] ), 
                         'n': eq.get_Xn_rule('r', 0, ['s','p'] ) },
               
                  's' : {'r' : eq.get_XX_rule('r', 's', 0, ['s','p'] ),
                         'n': eq.get_Xn_rule('s', 0, ['s','p'] ) },
    
                  'n' : {'r' : eq.get_nX_rule('r', 0, ['s','p'] ),
                         's' : eq.get_nX_rule('s', 0, ['s','p'] ) } }

        rules1 = {'p' : {'a' : eq.get_XX_rule('a', 'p', 0, ['s','p'] ) },
    
                  'n' : {'a' : eq.get_nX_rule('a', 1, ['s','p'] ), #angle fix -> OK
                         's' : eq.get_nX_rule('s', 0, ['s','p'] ) } }
    
        polrls = mst.MdimStruct( [rules0, rules1], [ (0,0), (1,1)] )
    
        dm_sp = self.dm_s.prod(self.dm_p)
    
        symbs = mst.MdimStruct( ['s', 'p'], [ (0, 0), (1, 1) ] )
        words = mst.MdimStruct( ['sp'], [ (0, 1) ] )
    
        self.poltab = eq.Table(dm_sp, symbs, words, polrls, codmgens, dmgens)
        
        return


    def test_bsgen(self):
        
        dm_t = deep(self.dm_t)
        dm_f = deep(self.dm_f)
        
        dm_tt = dm_t.prod(dm_t)
        dm_ff = dm_f.prod(dm_f)
        
        symb_tt = deep(self.symb_tt)
        symb_ff = deep(self.symb_ff)
        
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        inds = [1,2]
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_tt)
        bsfc = bsgen.gen(dm_tt, inds)
        
        vres = sp.cdmvec(bsfc.fcval_at, dm_tt)[1]
        vw.plot_heatvec(dm_tt, vres)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_ff)
        bsfc = bsgen.gen(dm_ff, inds)
        
        vres = sp.cdmvec(bsfc.fcval_at, dm_tt)[1]
        vw.plot_heatvec(dm_tt, vres)
        
        bsgen = pdtab.get_bsgen(symb_ff, symb_tt)
        bsfc = bsgen.gen(dm_tt, inds)
        
        vres = sp.cdmvec(bsfc.fcval_at, dm_ff)[1]
        vw.plot_heatvec(dm_ff, vres)
        
        bsgen = pdtab.get_bsgen(symb_ff, symb_ff)
        bsfc = bsgen.gen(dm_ff, inds)
        
        vres = sp.cdmvec(bsfc.fcval_at, dm_ff)[1]
        vw.plot_heatvec(dm_ff, vres)
        
        return
        

    def test_base_chg_unif(self):
        
        dm_t = deep(self.dm_t)
        dm_f = deep(self.dm_f)
        
        dm_tt = dm_t.prod(dm_t)
        dm_ff = dm_f.prod(dm_f)
        
        symb_tt = deep(self.symb_tt)
        symb_ff = deep(self.symb_ff)
        
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        labels = deep(pdtab.bsubsp_labels)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_ff)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[2,2] = 0.7 + 1.2j
        
        sig_tt_ff = sp.Signal(labels, coeffs, bsgen)
        
        sig_ff_ff = pdtab.base_change(symb_tt, symb_ff, dm_ff, dm_tt, sig_tt_ff)
        
        vres = sp.sigvec(sig_ff_ff.fcval_at, dm_ff, dm_ff)[1]
        vw.plot_heatvec(dm_ff, vres)
        
        sig_tt_tt = pdtab.base_change(symb_ff, symb_tt, dm_ff, dm_ff, sig_ff_ff)
        
        vres = sp.sigvec(sig_tt_tt.fcval_at, dm_tt, dm_tt)[1]
        vw.plot_heatvec(dm_tt, vres)
        
        norm = op.vec_hsp_sc_prod(vres, vres, dm_tt)
        
        self.assertAlmostEqual(norm, (0.7+1.2j) * (0.7-1.2j), places = 5)
        
        return


    def test_base_chg_no_unif(self):
        
        nudm_t = deep(self.dm_t)
        dm_f = deep(self.dm_f)
        
        nudm_tt = nudm_t.prod(nudm_t)
        dm_ff = dm_f.prod(dm_f)
        
        symb_tt = deep(self.symb_tt)
        symb_ff = deep(self.symb_ff)
        
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        labels = deep(pdtab.bsubsp_labels)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_ff)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[0,1] = 0.7 + 1.2j
        
        coeffs[1,2] = 1.5 + 0.1j 
        
        sig_tt_ff = sp.Signal(labels, coeffs, bsgen)
        
        sig_ff_ff = pdtab.base_change(symb_tt, symb_ff, dm_ff, nudm_tt, sig_tt_ff)
        
        vres = sp.sigvec(sig_ff_ff.fcval_at, dm_ff, dm_ff)[1]
        
        vw.plot_heatvec(dm_ff, vres)
        
        norm = op.vec_hsp_sc_prod(vres, vres, dm_ff)
    
        expect = (0.7+1.2j)*(0.7-1.2j)+(1.5+0.1j)*(1.5-0.1j)
    
        self.assertAlmostEqual(norm, expect, places = 5)
    
        return
    

#    def test_product_crcr(self):
#        
#        #get prelim.
#        symb_tt = deep(self.symb_tt)
#        symb_ff = deep(self.symb_ff)
#        
#        #get test signal
#        sigC = deep(self.sigcr)
#        
#        #get table w. 1d cartes. grid
#        crtab = deep(self.crtab)
#        
#        #do product
#        pdtab = crtab.product(crtab)
#            
#        #check some base change functionality and show result plots
#        sigRes = pdtab.base_change(symb_tt, symb_ff, sigC)
#            
#        vres = sigRes.sigvec()[1]
#        norm = op.vec_hsp_sc_prod(vres, vres, sigRes.domain)
#            
#        self.assertAlmostEqual(norm, 1.0 + 0.0j)
#            
#        vw.plot_heatvec(sigRes.domain, vres)
#            
#        #BACKWARD
#        sigRecon = pdtab.base_change(symb_ff, symb_tt, sigRes)
#            
#        vres = sigRecon.sigvec()[1]
#        norm = op.vec_hsp_sc_prod(vres, vres, sigRecon.domain)
#            
#        self.assertAlmostEqual(norm, 1.0 + 0.0j)
#            
#        vw.plot_heatvec(sigRecon.domain, vres)
#        
#        return
    
    
#    def test_replace_restrict(self):
#        
#        #get prelim.
#        symb_tt = deep(self.symb_tt)
#        symb_ff = deep(self.symb_ff)
#        
#        symb_ra = deep(self.symb_ra)
#        symb_sp = deep(self.symb_sp)
#        
#        #get table w. 1d cartes. grid
#        crtab = deep(self.crtab)
#        
#        #get table w 2d polar grid
#        poltab = deep(self.poltab)
#        
#        #do product
#        pdtab = crtab.product(crtab)
#        pdtab = pdtab.product(pdtab)
#        
#        #replace second cr-cr pair with polar
#        pdtab.replace(poltab, (2,3), 0)
#            
#        #restrict product to left and right substructures
#        lhtab = pdtab.restrict( (0,1) )
#        rhtab = pdtab.restrict( (2,3) )
#        
#        #get test signals
#        lhsig = lhtab.base_signal(symb_tt, symb_tt, [2,2] )
#        rhsig = rhtab.base_signal(symb_ra, symb_sp, [1,1] )
#        
#        vw.plot_heatmap(lhsig, 'cartes')
#        vw.plot_heatmap(rhsig, 'polar')
#        
#        #check some base change functionality and show result plots in cr-cr slice
#        symb_ttra = symb_tt.prod(symb_ra) 
#        symb_ttsp = symb_tt.prod(symb_sp) 
#        symb_ffsp = symb_ff.prod(symb_sp)
#        
#        sig = pdtab.base_signal(symb_ttra, symb_ttsp, [2,2,1,1])
#        sigRes = pdtab.base_change(symb_ttra, symb_ffsp, sig)
#            
#        vres = sigRes.sigvec()[1]
#        norm = op.vec_hsp_sc_prod(vres, vres, sigRes.domain)
#            
#        self.assertAlmostEqual(norm, 1.0 + 0.0j, places = 5)
#        
#        #BACKWARD
#        sigRecon = pdtab.base_change(symb_ffsp, symb_ttra, sigRes)
#            
#        vres = sigRecon.sigvec()[1]
#        norm = op.vec_hsp_sc_prod(vres, vres, sigRecon.domain)
#            
#        self.assertAlmostEqual(norm, 1.0 + 0.0j, places = 5)
#        
#        return


#TODO: design cases where rules are used with prod./restr./repl.
      
if __name__ == '__main__':
    unittest.main()
    
