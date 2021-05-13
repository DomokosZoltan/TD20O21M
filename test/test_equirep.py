# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:56:50 2021

@author: Domokos Zoltán
"""


import numpy as np

from hypothesis import given, settings #, example
from hypothesis.strategies import sampled_from
from hypothesis.strategies import integers
from hypothesis.strategies import floats
from hypothesis.strategies import just
from hypothesis.strategies import data

import mstruct as mst
import generator as gn
import cycdiv as cd
import equirep as tb
import domain as dm
import grid as gd

import cmdloop as cl

import unittest


class TestTable(unittest.TestCase):
    
    def setUp(self):
        
        #linear SAMPLE objects        
        ct = 4

        lim_t = 31
        
        lbpar_t = cd.LbFixParams( [ct], [lim_t] )

        shf_t = 0.0
        
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
        
        gdS = gd.Grid(genS)
        
        #BOUNDS
        diff = [dt]
        shf = [shf_t]
        band = [band_t]
        
        ptpar = cd.PtFixParams(diff, shf, band)
        
        pt_t = cd.PtGen(ptpar)
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        genB = gn.Generator(mst_lb, mst_pt) 
        
        gdB = gd.Grid(genB)
                
        bandstr = mst.MdimStruct( [ [0.0, 1.0] ], [ (0,0) ] )
        
        symb_t = mst.MdimStruct( ['t'], [ (0,0) ] )
        
        dm_t = dm.Domain(gdS, gdB, None, None, True, bandstr)
        
        self.table = tb.Table(dm_t, symb_t, None, None, None, None, True)
        
        return
    
    
    cmd_list_all=['t','f','n','abc', None]
    
    
    #first verify command processing
    @given(sampled_from(cmd_list_all), sampled_from(cmd_list_all) )
    
    def test_proc_cmd(self, src, fix): 
        
        dmdef = cl.DmDef()
        
        cmd_list_valid=['t','f','n']
        
        if (src is None) or (fix is None):
            
            with self.assertRaises(cl.InvalidCommand):
                
                dmdef.cmdproc_depr(src, fix)
        
        elif (src not in cmd_list_valid) or (fix not in cmd_list_valid):
        
            with self.assertRaises(cl.InvalidCommand):
            
                dmdef.cmdproc_depr(src, fix)
        
        elif (src == fix):
        
            with self.assertRaises(cl.InvalidCommand):
            
                dmdef.cmdproc_depr(src, fix)
        
        else:
            
            cmd = dmdef.cmdproc_depr(src, fix)
            #we have exactly those keys that we need
            self.assertEqual(cmd.keys(), {'src','dst','fix'} )
            #src and fix are in the dict, and with proper keys
            self.assertEqual(cmd['src'], src)
            self.assertEqual(cmd['fix'], fix)
            #command 'dst' has a valid id
            self.assertTrue(cmd['dst'] in cmd_list_valid)
            #no key pairs made from {src,dst,fix} have equal ids
            #as values 
            self.assertNotEqual(cmd['src'], cmd['dst'] )
            self.assertNotEqual(cmd['src'], cmd['fix'] )
            self.assertNotEqual(cmd['dst'], cmd['fix'] )
            
        return
    
    
    #list of all valid command combinations
    cmd_list_cmb=[ {'src':'t','fix':'f','dst':'n'},
                   {'src':'t','fix':'n','dst':'f'},
                   {'src':'f','fix':'t','dst':'n'},
                   {'src':'f','fix':'n','dst':'t'},
                   {'src':'n','fix':'t','dst':'f'},
                   {'src':'n','fix':'f','dst':'t'} ]
    
    
    #get_table, one command step (from init)
    @given(just(cl.DmDef() ), floats(), sampled_from(cmd_list_cmb) )
    
    def test_update_table(self, dmdef, val, cmd):
        
        table = self.table
        
        thrown = False
        
        try:
        
            if cmd['dst'] == 'n':
            
                tb.icheck_X_Y(table.domains[cmd['fix'] ].X(0), val)
                nrule = tb.get_XX_rule(cmd['src'], cmd['fix'] )
                tb.ocheck_n(nrule(val, table.domains[cmd['fix'] ].X(0) )[0] )
            
            elif cmd['src'] == 'n':

                tb.icheck_X_n(table.domains[cmd['fix'] ].X(0), val)
                nXrule = tb.get_nX_rule(cmd['fix'] )
                tb.ocheck_Y(nXrule(val, table.domains[cmd['fix'] ].X(0) )[1][1] )
            
            else:
            
                tb.icheck_X_n(val, table.domains['t'].n(0) )
                Xnrule = tb.get_Xn_rule(cmd['src'] )
                tb.ocheck_Y(Xnrule(val, table.domains['t'].n(0) )[1][1] )
        
        except(ZeroDivisionError, ValueError):
        
            thrown=True
        
        if thrown:
            
            pass
            #with self.assertRaises(dm2.InvalidTableArgs):
            #    dmdef.set_cmd(table, val, cmd['src'], cmd['fix'])
        else:
        
            eps=0.0000001
            
            table = dmdef.set_cmd(table, val, cmd['src'], cmd['fix'])
            
            self.assertTrue(np.abs(1.0 - table.domains['t'].measure_at([0] )\
                / (table.domains['t'].X(0) / table.domains['t'].n(0) ) ) < eps)
            
            self.assertTrue(np.abs(1.0 - table.domains['t'].measure_at([0] )\
                                   / (1.0 / table.domains['f'].X(0) ) ) < eps)
            
            self.assertTrue(np.abs(1.0 - table.domains['f'].measure_at([0] )\
                / (table.domains['f'].X(0) / table.domains['t'].n(0) ) ) < eps)
            
            self.assertTrue(np.abs(1.0 - table.domains['f'].measure_at([0] )\
                                   / (1.0 / table.domains['t'].X(0) ) ) < eps)
            #TODO: check if initial table is unchanged
   
        return
    
    
    #get_table, sequence of steps (from init)
    @given(data() )
    
    def test_seq_get_table(self, data):
        
        cmd_list_cmb = [ {'src':'t','fix':'f','dst':'n'},
                         {'src':'t','fix':'n','dst':'f'},
                         {'src':'f','fix':'t','dst':'n'},
                         {'src':'f','fix':'n','dst':'t'},
                         {'src':'n','fix':'t','dst':'f'},
                         {'src':'n','fix':'f','dst':'t'} ]
        
        dmdef = cl.DmDef()
        table = self.table
        
        for k in range(0,10):
            
            cmd = data.draw(sampled_from(cmd_list_cmb) )
            #print(cmd)
            if 'n' == cmd['src']:
            
                val = data.draw(integers(max_value=1000) )
            
            else:
            
                val = data.draw(floats(max_value=100000.0) )
            
            thrown = False
            
            try:
                
                if cmd['dst'] == 'n':
                
                    tb.icheck_X_Y(table.domains[cmd['fix'] ].X(0), val)
                    nrule = tb.get_XX_rule(cmd['src'], cmd['fix'] )
                    tb.ocheck_n(nrule(val, table.domains[cmd['fix'] ].X(0) )[0] )
                
                elif cmd['src'] == 'n':
    
                    tb.icheck_X_n(table.domains[cmd['fix'] ].X(0), val)
                    nXrule = tb.get_nX_rule(cmd['fix'] )
                    tb.ocheck_Y(nXrule(val, table.domains[cmd['fix'] ].X(0) )[1][1] )
                
                else:
                
                    tb.icheck_X_n(val, table.domains['t'].n(0) )
                    Xnrule = tb.get_Xn_rule(cmd['src'] )
                    tb.ocheck_Y(Xnrule(val, table.domains['t'].n(0) )[1][1] )
            
            except(ZeroDivisionError, ValueError):
            
                thrown = True
            
            if thrown:
                
                pass
                #with self.assertRaises(dm2.InvalidTableArgs):
                #   dmdef.set_cmd(table, val, cmd['src'], cmd['fix'])
            else:
                
                eps = 0.0000001
                
                table = dmdef.set_cmd(table, val, cmd['src'], cmd['fix'])
                
                self.assertTrue(np.abs(1.0-table.domains['t'].measure_at([0] )\
                    / (table.domains['t'].X(0)/table.domains['t'].n(0)))<eps)
                
                self.assertTrue(np.abs(1.0-table.domains['t'].measure_at([0] )\
                                       / (1.0/table.domains['f'].X(0)))<eps)
                
                self.assertTrue(np.abs(1.0-table.domains['f'].measure_at([0] )\
                    / (table.domains['f'].X(0)/table.domains['t'].n(0)))<eps)
                
                self.assertTrue(np.abs(1.0-table.domains['f'].measure_at([0] )\
                                       / (1.0/table.domains['t'].X(0)))<eps)
                
        return


if __name__ == '__main__':
    unittest.main()