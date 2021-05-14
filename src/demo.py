# -*- coding: utf-8 -*-
"""
Created on Thu May 13 22:17:00 2021

@author: Domokos Zoltán
"""

import numpy as np

import mstruct as mst
import domain as dm
import equirep as eq
import sglspace as sp
import sglutils as su
import view as vw

#This file briefly demonstartes some features of the framework. The goal was 
#to build a tool that simplifies the following regularly recurring tasks in
#studying the mathematical background of general DSP problems:

#1. Defining a multidimensional grid that can be used for sampling or tiling 
#   a given domain. (E.g. cartesian and polar grids are implemented.)

#2. Defining a domain that models a vector space of arbitrary dimensionality,
#   with a tiling, a sampling pattern and a measure associated to it. 

#3. Defining mapping of samples of a domain to elements of an associated codomain.

#4. Representing sampled functions (signals) as a linear comabination of signals
#   that form an orthonormed base. (E.g. complex exponentials, Bessel functions,
#   etc.)

#5. Support base change between pairs of orthonormed base function sets.

#6. Define relations that automatically maintain sufficient sampling rate when
#   a given domain is rescaled or adapt scaling when the sampling rate changes. 

#7. Implement utilities to analyze signals (4). (E.g. approximation of an inner
#   product, shift relations, etc.) 

#All the objects implementing the features above have to be able to handle:

#1. Arbitrary dimensionality.

#2. Concatenation (or 'product') instances of the underlying mathematical
#   structure, if applicable.

#3. Replacement of instances of the underlying mathematical structure in any
#   dimensions when it's possible.

#4. Restriction to subdimensions, if possible. (E.g. pick a 1d subgrid from a 
#   grid that is a product of N 1d grids. Counter example: a polar grid, where
#   such reduction is not supported.)

#5. All the data content of the objects have to be generated just in time, follo-
#   wing a functional programming style, and numpy arrays are used only in case
#   it is strictly necessary. (E.g. plotting, storing irregular patterns of 
#   non-uniform sampling, using a sampled period of a signal to perform FFT, etc.)

#6. Easy integration of special base function sets. Therefore base change is 
#   generally done by solving a system of linear equations, and faster methods
#   are used only in case they are well-known. (E.g. FFT/IFFT. )

#7. Partial transforms in arbitrary dimensions.

#The following examples highlight those features that have been implemented
#so far. Execution time is relatively high, due to reasons explained above. (See 5.)

##EXAMPLE 1: domain definition, codomain definition, setting up a uniformly 
##sampled signal that is bandlimited and discrete in the frequency domain.
##Transformation to frequency space and back is done using the equation system
##approach. (There is a scipy FFTN wrapper in 'operators.py' that is much faster.)
n1 = 8
n2 = 4

T1 = 1.0
T2 = 1.0

F1 = n1
F2 = n2

band_offs = 0.0

band_t = [0.0, T1]

band_f1 = [0.0, F1]
band_f2 = [0.0, F2]

#Initially one dimensional domains.
dm_tt = dm.factory( [band_t], [n1], [band_offs], 'cartes')
dm_t2 = dm.factory( [band_t], [n2], [band_offs], 'cartes')

#Helper object that we will need later. 
pdtab = eq.Table(dm_tt)
pdtab = pdtab.product(eq.Table(dm_t2) )

#Apply product.
dm_tt = dm_tt.prod(dm_t2)

#Same for the frequency domain.       
dm_ff = dm.factory( [band_f1], [n1], [band_offs], 'cartes')

dm_ff = dm_ff.prod( dm.factory( [band_f2], [n2], [band_offs], 'cartes') )

#Aquire the object that will label the elements of the set of orthonormed 
#base functions.
labels = pdtab.bsubsp_labels
       
#Symbol structures
dim1 = (0,0)
dim2 = (1,1)

symb_tt = mst.MdimStruct( ['t', 't'], [ dim1, dim2 ] )
symb_ff = mst.MdimStruct( ['f', 'f'], [ dim1, dim2 ] )

#Aquire generator object that will be used to generate base function values 
#needed for the linear combination. Symbols are used to tell which generator
#is wanted. Important: tt - ff means that coordinatization is done in time
#space and we will query the elements of the frequency space basis, thus the 
#signal will be defined as a linear combination of plane waves in time.   
bsgen = pdtab.get_bsgen(symb_tt, symb_tt)
  
#Coefficients for the linear combination.    
coeffs = np.zeros(shape = (n1, n2), dtype = np.complex64)
    
#Here tt - tt means that the signal is a linear combination of time domain base 
#functions and coordinatization is also done in this space. Base functions are app-
#roximated currently by step functions. 
bsgen = pdtab.get_bsgen(symb_tt, symb_tt)

#Coefficients are set in a way that the time domain signal form will be a planar
#wave with frequencies (2,1).
phase = 1.0 + 0.0j
ampl = 1.0

f1 = 2
f2 = 1

exp1 = su.expfc(T1, f1, phase, ampl, 't')
exp2 = su.expfc(T2, f2, phase, ampl, 't')

coeffs = np.zeros(shape = (n1, n2), dtype = np.complex64)

for i in range(0,n1):
    for k in range(0,n2):
        coeffs[i,k] = exp1( [i*T1/n1] ) * exp2( [k*T2/n2] )
        
sgl_tt_tt = sp.Signal(labels, coeffs, bsgen)

#Vectorization is needed for the plot. Plotting is done by a simple Matplotlib
#wrapper.
vsign = sp.sigvec(sgl_tt_tt.fcval_at, dm_tt, dm_tt)[1]  
vw.plot_heatvec(dm_tt, vsign)

#Get representation in frequency space. In this example DFT could be done, but
#the more generic (and slower) method of solving a system of linear equations 
#is used.
sgl_res = pdtab.base_change(symb_tt, symb_ff, dm_tt, dm_tt, sgl_tt_tt)

vres = sp.sigvec(sgl_res.fcval_at, dm_ff, dm_ff)[1]
vw.plot_heatvec(dm_ff, vres)

#To get a planar wave decomposition in time domain, now it's sufficient to 
#take the coefficients of the transformed signal. ff - tt below means that 
#cordinatization is done in time domain, but the signal is constructed as a 
#linear combination  of members of the other base function set.
bsgen = pdtab.get_bsgen(symb_tt, symb_ff)

#Create and plot the signal.
sgl_tt_ff = sp.Signal(labels, sgl_res.coeffs, bsgen)

vres = sp.sigvec(sgl_tt_ff.fcval_at, dm_ff, dm_tt)[1]
vw.plot_heatvec(dm_tt, vres)

#Now interpolation can be done, by a simple update of the sampling rate in time
#domain.
n1u = 32
n2u = 32

dm_tt.update(n1u, band_t, 0)
dm_tt.update(n2u, band_t, 1)

vres = sp.sigvec(sgl_tt_ff.fcval_at, dm_ff, dm_tt)[1]
vw.plot_heatvec(dm_tt, vres)

#EXAMPLE 2: as a more interesting example, zero order Bessel functions of the 
#first kind multiplied by angular waves can provide an orthonormed base function
#set when a polar grid defines the set of samples over the domain. 
nr = 4
na = 8

eps = 0.00001

R = 1.0
A = 2.0 * np.pi

band_rad = [0, R]
band_ang = [0, A]

S = nr
P = na / (2.0 * np.pi)

band_s = [0, S]
band_p = [0, P]

band_offs_sp = 0.0

dm_ra = dm.factory( [band_rad, band_ang], [nr, na], [eps, eps], 'polar')
        
dm_s = dm.factory( [band_s], [nr], [band_offs_sp], 'cartes')      
dm_p = dm.factory( [band_p], [na], [band_offs_sp], 'cartes')

dm_sp = dm_s.prod(dm_p)

symb_ra = mst.MdimStruct( ['ra'], [ (0,1) ] )       
symb_sp = mst.MdimStruct( ['sp'], [ (0,1) ] )

tab = eq.get_bessel_tab(dm_s, dm_p)

phase = 1.0 + 0.0j
ampl = 1.0

s = 1
fa = 2 

bessfc = su.spher_bessfc(s, fa, R, A, phase, ampl, 'ra', alph = 0.0)

coeffs = np.zeros(shape = (nr, na), dtype = np.complex64)

for i in range(0, nr):
    for k in range(0, na):
        x = i*R/nr * np.cos(k*2.0*np.pi/na) 
        y = i*R/nr * np.sin(k*2.0*np.pi/na) 
        coeffs[i,k] = bessfc([x, y])
        
labels = tab.bsubsp_labels
        
bsgen = tab.get_bsgen(symb_ra, symb_ra)

sig_ra_ra = sp.Signal(labels, coeffs, bsgen)
        
svec_ra_ra = sp.sigvec(sig_ra_ra.fcval_at, dm_ra, dm_ra)[1]
vw.plot_heatvec(dm_ra, svec_ra_ra, 'polar')

sgl_res = tab.base_change(symb_ra, symb_sp, dm_ra, dm_ra,  sig_ra_ra, md = 'polar')
        
vres = sp.sigvec(sgl_res.fcval_at, dm_sp, dm_sp)[1]        
vw.plot_heatvec(dm_sp, vres, 'cartes')

bsgen = tab.get_bsgen(symb_ra, symb_sp)

#Create and plot the signal.
sgl_ra_sp = sp.Signal(labels, sgl_res.coeffs, bsgen)

vres = sp.sigvec(sgl_ra_sp.fcval_at, dm_sp, dm_ra)[1]
vw.plot_heatvec(dm_ra, vres, 'polar')

nru = 32
nau = 32

dm_ra = dm.factory( [band_rad, band_ang], [nru, nau], [eps, eps], 'polar')

vres = sp.sigvec(sgl_ra_sp.fcval_at, dm_sp, dm_ra)[1]
vw.plot_heatvec(dm_ra, vres, 'polar')

















