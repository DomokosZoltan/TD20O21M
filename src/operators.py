# -*- coding: utf-8 -*-
"""
Created on Sun May  9 21:57:42 2021

@author: Domokos Zoltán
"""

import numpy as np

import sglspace as sg


def vec_hsp_sc_prod(sgvec_lhs, sgvec_rhs, dm):
    
    res = 0.0

    with np.nditer(sgvec_rhs, flags = ['multi_index'],\
                   op_flags=['readonly'], order='C') as it:
        
        for elm in it:
        
            np_inds = it.multi_index
            inds = list(np_inds)
    
            sc = dm.measure_at(inds)
    
            res += sc * sgvec_lhs[np_inds] * np.conjugate( sgvec_rhs[np_inds] )
    
    return res


#DFT/IDFT wrappers with unitary scaling.
def forward_fft(sgvec, dm_src, dm_dst):
    
    unit_sign = sg.scalev(sgvec, dm_src, 'sqrt(dx)')
    
    unit_dfted = np.fft.fft(unit_sign, norm = 'ortho')
    
    dfted = sg.scalev(unit_dfted, dm_dst, 'sqrt(1/dx)')
    
    return dfted


def inverse_fft(dfted, dm_src, dm_dst):
    
    unit_dfted = sg.scalev(dfted, dm_src, 'sqrt(dx)')
    
    unit_sign = np.fft.ifft(unit_dfted, norm = 'ortho')
    
    sign = sg.scalev(unit_sign, dm_dst, 'sqrt(1/dx)')
    
    return sign


#Multidimensional generalizations. Caller have to fill symbs_src and symbs_dst.
#If for an idx: idx not in prd_axes -> symbs_src[idx] == symbs_dst[idx] shall
#always hold. Note: there is no restriction for the choice regarding the dmsymb
#in such cases.
def forward_fftn(vsign, dm_src, dm_dst, prd_axes):
    
    unit_vsign = sg.scalev(vsign, dm_src, 'sqrt(dx)')
    
    unit_vdfted = np.fft.fftn(unit_vsign, axes = prd_axes, norm = 'ortho')
    
    vdfted = sg.scalev(unit_vdfted, dm_dst, 'sqrt(1/dx)')
    
    return vdfted


def inverse_fftn(vdfted, dm_src, dm_dst, prd_axes):
    
    unit_vdfted = sg.scalev(vdfted,dm_src, 'sqrt(dx)')
    
    unit_vsign = np.fft.ifftn(unit_vdfted, axes = prd_axes, norm = 'ortho')
    
    vsign = sg.scalev(unit_vsign, dm_dst, 'sqrt(1/dx)')
    
    return vsign


def shift(self, shf):
     #TODO
     pass