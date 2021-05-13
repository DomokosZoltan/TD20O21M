# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 14:48:20 2021

@author: Domokos Zoltán
"""

import numpy as np
from copy  import deepcopy as deep

import mstruct as mst
import domain as dm
import sglspace as sg
import sglgen as sgg
import equirep as eq
import view as vw

import operators as op

import unittest


class TestOperators(unittest.TestCase):

    def setUp(self):
        
        st = [0.1, 0.39, 0.47, 0.82]
        bt = [0.0, 0.27, 0.43, 0.7]
        
        self.nufdm_t = dm.nufactory( [ [0.0, 1.0] ], [4], st, bt, [0.0], 'cartes')
        
        self.dm_t = dm.factory( [ [0.0, 1.0] ], [4], [0.0], 'cartes')
        
        self.dm_f = dm.factory( [ [0.0, 4.0] ], [4], [0.0], 'cartes')
                    
        self.dm_ra = dm.factory([ [0, 1.0], [0, 2.0 * np.pi] ], [4, 4], [0.0, 0.0], 'polar')
        
        self.dm_ra_shf = dm.factory([ [0, 1.0], [0, 2.0 * np.pi] ], [4, 4], [0.00001, 0.00001], 'polar')
        
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
        
        self.symb_t = mst.MdimStruct( ['t'], [ (0,0) ] )
        
        self.crtab = eq.Table(self.dm_t, self.symb_t, self.symb_t, crrls, gens, dmgens)
        
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
        
    
    def test_FFTN(self):
       
        symb_tt = deep(self.symb_tt)
        
        dm_t = deep(self.dm_t)
        dm_f = deep(self.dm_f)
        
        dm_tt = deep(dm_t)        
        dm_tt = dm_tt.prod(dm_t)
        
        dm_ff = deep(dm_f)        
        dm_ff = dm_ff.prod(dm_f)
        
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        labels = deep(pdtab.bsubsp_labels)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_tt)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[1,1] = 1.0 + 0.0j
        
        sig_tt_tt = sg.Signal(labels, coeffs, bsgen)
        
        vsign = sg.sigvec(sig_tt_tt.fcval_at, dm_tt, dm_tt)[1]
        
        vw.plot_heatvec(dm_tt, vsign)
        
        norm_tt = op.vec_hsp_sc_prod(vsign, vsign, dm_tt)
        
        sel_axs = [0,1]   
        vdfted = op.forward_fftn(vsign, dm_tt, dm_ff, sel_axs)
        
        vw.plot_heatvec(dm_ff, vdfted)
        
        #print measure
        norm_ff = op.vec_hsp_sc_prod(vdfted, vdfted, dm_ff)
        
        self.assertAlmostEqual(norm_tt, 1.0 + 0.0j)
        
        self.assertAlmostEqual(norm_ff, 1.0 + 0.0j)
        
        return
    
    
    def test_NUFFTN(self):
        
        symb_tt = deep(self.symb_tt)
        symb_ff = deep(self.symb_ff)
        
        nufdm_tt = deep(self.nufdm_t) 
        nufdm_tt = nufdm_tt.prod(self.nufdm_t)
        
        dm_t = deep(self.dm_t)
        dm_f = deep(self.dm_f)
   
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        labels = deep(pdtab.bsubsp_labels)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_ff)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[1,2] = 1.0 + 0.0j
        
        sig_tt_ff = sg.Signal(labels, coeffs, bsgen)
        
        dm_tt = dm_t.prod(dm_t)
        dm_ff = dm_f.prod(dm_f)
        
        vsign = sg.sigvec(sig_tt_ff.fcval_at, dm_ff, nufdm_tt)[1]
        
        vw.plot_heatvec(nufdm_tt, vsign)
        
        norm = op.vec_hsp_sc_prod(vsign, vsign, nufdm_tt)
        
        gendict = {'t' : {'t' : sgg.get_codm_gen('t', 't'),
                          'f' : sgg.get_codm_gen('f','t') },
    
                   'f' : {'t' : sgg.get_codm_gen('t', 'f'),
                          'f' : sgg.get_codm_gen('f', 'f') } }
        
        gens = mst.MdimStruct( [gendict, gendict], [ (0,0), (1,1) ] )
        
        dmgendict = {'t' : {'t' : None,
                            'f' : eq.get_cartes_gen('t', 'f') },
    
                     'f' : {'t' : eq.get_cartes_gen('f', 't'),
                            'f' : None } }
    
        dmgens = mst.MdimStruct( [dmgendict, dmgendict], [ (0,0), (1,1) ] )
        
        table = eq.Table(dm_tt, symb_tt, None, None, gens, dmgens)
        
        sigRes = table.base_change(symb_tt, symb_ff, dm_ff, nufdm_tt, sig_tt_ff)
        
        vres = sg.sigvec(sigRes.fcval_at, dm_ff, dm_ff)[1]
        norm = op.vec_hsp_sc_prod(vres, vres, dm_ff)
        
        self.assertAlmostEqual(norm, 1.0 + 0.0j, 5)
        
        vw.plot_heatvec(dm_ff, vres)
    
        sigRecon = table.base_change(symb_ff, symb_tt, dm_ff, dm_ff, sigRes)
    
        vres = sg.sigvec(sigRecon.fcval_at, dm_tt, dm_tt)[1]
        
        norm = op.vec_hsp_sc_prod(vres, vres, dm_tt)
        
        self.assertAlmostEqual(norm, 1.0 + 0.0j, 5)
        
        vw.plot_heatvec(dm_tt, vres)
    
        return
    
    
    def test_cartes_basefcs(self):
     
        symb_tt = deep(self.symb_tt)
        symb_ff = deep(self.symb_ff)
        
        dm_t = deep(self.dm_t)
        
        dm_tt = dm_t.prod(dm_t)
        
        dm_f = deep(self.dm_f)
        
        dm_ff = dm_f.prod(dm_f)
        
        crtab = deep(self.crtab)
        
        pdtab = crtab.product(crtab)
        
        labels = deep(pdtab.bsubsp_labels)
        
        bsgen = pdtab.get_bsgen(symb_tt, symb_ff)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[1,0] = 1.0 + 0.0j
        
        sig_tt_ff = sg.Signal(labels, coeffs, bsgen)
        
        vres = sg.sigvec(sig_tt_ff.fcval_at, dm_ff, dm_tt)[1]
        
        norm = op.vec_hsp_sc_prod(vres, vres, dm_tt)
        
        self.assertAlmostEqual(norm, 1.0 + 0.0j)
        
        vw.plot_heatvec(dm_tt, vres)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[2,2] = 1.0 + 0.0j
        
        sig_tt_ff_2 = sg.Signal(labels, coeffs, bsgen)
        
        vres2 = sg.sigvec(sig_tt_ff_2.fcval_at, dm_ff, dm_tt)[1]
        
        norm = op.vec_hsp_sc_prod(vres2, vres2, dm_tt)
        
        self.assertAlmostEqual(norm, 1.0 + 0.0j)
        
        vw.plot_heatvec(dm_tt, vres2)
        
        norm = op.vec_hsp_sc_prod(vres, vres2, dm_tt)
        
        self.assertAlmostEqual(norm, 0.0 + 0.0j)
        
        return
    
    
    def test_polar_basefcs(self):
        
        symb_ra = deep(self.symb_ra)
        symb_sp = deep(self.symb_sp)
        
        dm_s = deep(self.dm_s)     
        dm_p = deep(self.dm_p)
        
        dm_sp = dm_s.prod(dm_p)       

        dm_ra = deep(self.dm_ra)

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
    
        
        table = eq.Table(dm_sp, symb_sp, None, None, codmgens, dmgens)
        
        labels = deep(table.bsubsp_labels)
        
        bsgen = table.get_bsgen(symb_ra, symb_ra)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[3,3] = 1.0 + 0.0j
        
        sig_ra_ra = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_ra = sg.sigvec(sig_ra_ra.fcval_at, dm_ra, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, svec_ra_ra, 'polar')
        
        bsgen = table.get_bsgen(symb_sp, symb_sp)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[3,3] = 1.0 + 0.0j
        
        sig_sp_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_sp_sp = sg.sigvec(sig_sp_sp.fcval_at, dm_sp, dm_sp)[1]
        
        vw.plot_heatvec(dm_sp, svec_sp_sp, 'cartes')
        
        bsgen = table.get_bsgen(symb_sp, symb_ra)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[1,1] = 1.0 + 0.0j
        
        sig_ra_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_sp = sg.sigvec(sig_ra_sp.fcval_at, dm_ra, dm_sp)[1]
        
        vw.plot_heatvec(dm_sp, svec_ra_sp, 'cartes')
        
        bsgen = table.get_bsgen(symb_ra, symb_sp)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[1,1] = 1.0 + 0.0j
        
        sig_sp_ra = sg.Signal(labels, coeffs, bsgen)
        
        svec_sp_ra = sg.sigvec(sig_sp_ra.fcval_at, dm_sp, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, svec_sp_ra, 'polar')
  
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[3,1] = 1.0 + 0.0j
        
        return
    
    
    def test_polar_basechng(self):
        
        symb_ra = deep(self.symb_ra)
        symb_sp = deep(self.symb_sp)
        
        dm_s = deep(self.dm_s)     
        dm_p = deep(self.dm_p)       
 
        dm_sp = dm_s.prod(dm_p)
        
        dm_ra = self.dm_ra_shf
        
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
    
        symbs = mst.MdimStruct( ['s', 'p'], [ (0, 0), (1, 1) ] )
        words = mst.MdimStruct( ['sp'], [ (0, 1) ] )   
        
        table = eq.Table(dm_sp, symbs, words, None, codmgens, dmgens)
        
        labels = table.bsubsp_labels
        
        bsgen = table.get_bsgen(symb_ra, symb_sp)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[2,1] = 1.0 + 0.0j
        
        sig_ra_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_sp = sg.sigvec(sig_ra_sp.fcval_at, dm_sp, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, svec_ra_sp, 'polar')
        
        sigRes = table.base_change(symb_ra, symb_sp, dm_sp, dm_ra,  sig_ra_sp)
        
        vres = sg.sigvec(sigRes.fcval_at, dm_sp, dm_sp)[1]
        
        vw.plot_heatvec(dm_sp, vres, 'cartes')
    
        sigRecon = table.base_change(symb_sp, symb_ra, dm_sp, dm_sp, sigRes)
    
        vres = sg.sigvec(sigRecon.fcval_at, dm_ra, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, vres, 'polar')
        
        return
    
    
    def test_polar_rls_update(self):
        
        dm_ra_upd = dm.factory([ [0, 1.0], [0, 2.0 * np.pi] ], [4, 8], [0.0, 0.0], 'polar')
        
        symb_ra = deep(self.symb_ra)
        symb_sp = deep(self.symb_sp)
        
        dm_s = deep(self.dm_s)        
        dm_p = deep(self.dm_p)
        
        dm_sp = dm_s.prod(dm_p)
        
        dm_ra = self.dm_ra_shf
        
        symbs = mst.MdimStruct( ['s', 'p'], [ (0, 0), (1, 1) ] )
        words = mst.MdimStruct( ['sp'], [ (0, 1) ] )

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
    
        rls = mst.MdimStruct( [rules0, rules1], [ (0,0), (1,1)] )
    
        #table = eq.Table(sigSP.domain, coord_symb, coord_words, rls, codmgens, dmgens)
        table = eq.Table(dm_sp, symbs, words, rls, codmgens, dmgens)
        
        print("BEFORE UPDATE 'ra' in 'ra' : " )
        labels = table.bsubsp_labels
        
        bsgen = table.get_bsgen(symb_ra, symb_ra)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[3,2] = 1.0 + 0.0j
        
        sig_ra_ra = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_ra = sg.sigvec(sig_ra_ra.fcval_at, dm_ra, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, svec_ra_ra, 'polar')

        
#        base_ra_ra = table.base_signal(symb_ra, symb_ra, [3,2] )
#        vw.plot_heatmap(base_ra_ra, 'polar')
        
        print("BEFORE UPDATE 'sp' in 'sp' : " )
        bsgen = table.get_bsgen(symb_sp, symb_sp)
        
        coeffs = np.zeros(shape = (4,4), dtype = np.complex64)
        
        coeffs[3,2] = 1.0 + 0.0j
        
        sig_sp_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_sp_sp = sg.sigvec(sig_sp_sp.fcval_at, dm_sp, dm_sp)[1]
        
        vw.plot_heatvec(dm_sp, svec_sp_sp, 'cartes')
        
        #Now update ns:
        cmd = {'src':'nn','fix':'ra','dst':'sp'}
        table.update(8, 1, cmd)
        
        print("AFTER UPDATE 'ra' in 'ra' : " )
        labels = table.bsubsp_labels
        
        bsgen = table.get_bsgen(symb_ra, symb_ra)
        
        coeffs = np.zeros(shape = (4,8), dtype = np.complex64)
        
        coeffs[3,4] = 1.0 + 0.0j
        
        sig_ra_ra = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_ra = sg.sigvec(sig_ra_ra.fcval_at, dm_ra_upd, dm_ra_upd)[1]
        
        vw.plot_heatvec(dm_ra_upd, svec_ra_ra, 'polar')
        
        print("AFTER UPDATE 'sp' in 'sp' : " )
        bsgen = table.get_bsgen(symb_sp, symb_sp)
        
        coeffs = np.zeros(shape = (4,8), dtype = np.complex64)
        
        coeffs[3,4] = 1.0 + 0.0j
        
        sig_sp_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_sp_sp = sg.sigvec(sig_sp_sp.fcval_at, dm_sp, dm_sp)[1]
        
        vw.plot_heatvec(dm_sp, svec_sp_sp, 'cartes')
        
        return


    def test_bess_basechng(self):
        
        symb_ra = deep(self.symb_ra)
        symb_sp = deep(self.symb_sp)
        
        dm_s = deep(self.dm_s)     
        dm_p = deep(self.dm_p)       
 
        dm_sp = dm_s.prod(dm_p)
        
        dm_ra = self.dm_ra_shf
        
        #Then a table
        codmgens = mst.MdimStruct( [ {'ra' : {'ra' : sgg.get_bess_gen('ra', 'ra'), 
                                              'sp' : sgg.get_bess_gen('sp','ra') },
                    
                                    'sp' : {'ra' : sgg.get_bess_gen('ra', 'sp'),
                                            'sp' : sgg.get_bess_gen('sp', 'sp') } } ],
                                    [ (0,1) ] )
        
        dmgens = mst.MdimStruct( [ {'ra' : {'ra' : None,
                                            'sp' : eq.get_polar_gen('ra', 'sp') },
    
                                  'sp' : {'ra' : eq.get_polar_gen('sp', 'ra'),
                                          'sp' : None } } ],
                                    [ (0,1) ] )
    
        symbs = mst.MdimStruct( ['s', 'p'], [ (0, 0), (1, 1) ] )
        words = mst.MdimStruct( ['sp'], [ (0, 1) ] )   
        
        table = eq.Table(dm_sp, symbs, words, None, codmgens, dmgens)
        
        #Finally request a base function...
        labels = table.bsubsp_labels
        
        bsgen = table.get_bsgen(symb_ra, symb_sp)
        
        coeffs = np.zeros(shape = (4, 4), dtype = np.complex64)
        
        coeffs[1,2] = 1.0 + 0.0j
        
        #coeffs[3,0] = 1.2 + 0.0j
        
        sig_ra_sp = sg.Signal(labels, coeffs, bsgen)
        
        svec_ra_sp = sg.sigvec(sig_ra_sp.fcval_at, dm_sp, dm_ra)[1]
        
        vw.plot_heatvec(dm_ra, svec_ra_sp, 'polar')
        
        #Change the base
        sigRes = table.base_change(symb_ra, symb_sp, dm_sp, dm_ra,  sig_ra_sp, md = 'polar')
        
        vres = sg.sigvec(sigRes.fcval_at, dm_sp, dm_sp)[1]
        
        print('sp in sp, after base change')
        vw.plot_heatvec(dm_sp, vres, 'cartes')
    
        #... and revert it.
        sigRecon = table.base_change(symb_sp, symb_ra, dm_sp, dm_sp, sigRes)
    
        vres = sg.sigvec(sigRecon.fcval_at, dm_ra, dm_ra)[1]
        
        print('sp in ra, after base change inversion')
        vw.plot_heatvec(dm_ra, vres, 'polar')
        
        return


if __name__ == '__main__':
    unittest.main()
    
    