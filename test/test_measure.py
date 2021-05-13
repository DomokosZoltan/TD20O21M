# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:55:18 2021

@author: domok
"""

import unittest

import numpy as np

import domain as dm


class TestMeasure(unittest.TestCase):

    def setUp(self):
        pass


    def test_polar(self):
        
        nr = 16
        
        na = 16
        
        R = 2.0
        
        dm_ra = dm.factory([ [0, R], [0, 2.0 * np.pi] ], [nr, na], [0.0, 0.0], 'polar')
        
        m = 0.0
        
        for r in range(0, nr):
            
            for a in range(0, na):
                
                m += dm_ra.measure_at( [r, a] )
                
        self.assertAlmostEqual(m, (R**2.0 * np.pi) )


if __name__ == '__main__':
    unittest.main()