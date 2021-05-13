# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 19:58:42 2021

@author: Domokos Zoltán
"""


from copy import deepcopy as deep
import numpy as np

import mstruct as mst
import generator as gn
import grid as gd
import measure as ms
import volume as vol

import cycdiv as cd
import nufcycdiv as ncd
import polardiv as pl


def nufactory(bandlims, ns, sdt = None, bdt = None, shf = [0.0], gdtype = 'cartes'):

    res = None
    
    if 'cartes' == gdtype: 

        ct = ns[0]

        lim_t = ct * 4 - 1
        
        lbpar_t = ncd.LbFixParams( [ct], [lim_t] )

        shf_t = shf[0]
        
        T = bandlims[0][1]
        
        if sdt is None:
        
            diff = [T / ct for dt in range(0, ct) ]
            
        else:
            
            diff = deep(sdt)
        
        shf = [shf_t]
        band = [T]
        
        ptpar = ncd.PtFixParams(diff, shf, band)
        
        lb_t = ncd.LbGen(lbpar_t)
        
        pt_t = ncd.PtGen(ptpar)
        
        dim_lb = (0, 0)
        
        obj_lb = [lb_t]
        dims_lb = [dim_lb]
        
        mst_lb = mst.MdimStruct(obj_lb, dims_lb) 
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        samps = gd.Grid(gen)
        
        shf = [0.0]
        
        if bdt is None:
        
            diff = [T / ct for dt in range(0, ct) ]
            
        else:
            
            diff = deep(bdt)
        
        ptpar = ncd.PtFixParams(diff, shf, band)

        pt_t = ncd.PtGen(ptpar)
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        bounds = gd.Grid(gen)
                    
        bandstr = mst.MdimStruct( [ [0, T] ], [ (0,0) ] )
        
        res = Domain(samps, bounds, None, None, True, bandstr)
        
    return res
        
    
def factory(bandlims, ns, shf = [0.0], gdtype = 'cartes'):
    
    res = None
    
    if 'cartes' == gdtype: 

        ct = ns[0]

        lim_t = ct * 4 - 1
        
        lbpar_t = cd.LbFixParams( [ct], [lim_t] )

        shf_t = shf[0]
        
        T = bandlims[0][1]
        
        dt = T / ct
        
        diff = [dt]
        shf = [shf_t]
        band = [T]
        
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
        
        samps = gd.Grid(gen)
        
        shf = [0.0]
        
        ptpar = cd.PtFixParams(diff, shf, band)

        pt_t = cd.PtGen(ptpar)
        
        dims_pt = [ (0, 0) ]
        obj_pt = [pt_t]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        bounds = gd.Grid(gen)
                    
        bandstr = mst.MdimStruct( [ [0, T] ], [ (0,0) ] )
        
        res = Domain(samps, bounds, None, None, True, bandstr)
        
    elif 'polar' == gdtype:
        
        R = bandlims[0][1]
        A = bandlims[1][1]
        
        cr = int( ns[0] )
        ca = int( ns[1] )
        
        lim_r = cr * 4 - 1
        lim_a = ca * 2 - 1
        
        lbpar_r = pl.LbFixParams( [cr], [lim_r] )
        lbpar_a = pl.LbFixParams( [ca], [lim_a] )
        
        shf_r = shf[0]
        shf_a = shf[1]
        
        dr = R / cr
        da = A / ca
        
        diff = [dr, da]
        shf = [shf_r, shf_a]
        band = [R, A]
        
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
        
        samps = gd.Grid(gen)
        
        shf_r = 0.0
        shf_a = 0.0
        
        ptpar = pl.PtFixParams(diff, shf, band)
        
        pol_pt = pl.PolarPtGen(ptpar)
        
        dims_pt = [ (0, 1) ]
        obj_pt = [pol_pt]
        
        mst_pt = mst.MdimStruct(obj_pt, dims_pt)
        
        gen = gn.Generator(mst_lb, mst_pt)
        
        bounds = gd.Grid(gen)

        volume = vol.Volume( mst.MdimStruct( [vol.sector], [ (0,1) ] ) )
        
        meas = ms.Measure( mst.MdimStruct( [ms.sector], [ (0,1) ] ) )
        
        bandstr = mst.MdimStruct( [ [0, R], [0, A] ], [ (0,0), (1,1) ] )
        
        res = Domain(samps, bounds, volume, meas, True, bandstr)
        
    else:
       
        pass 
    
    return res


class Domain:

    def __init__(self, samps = None,
                 bounds = None,
                 volume = None,
                 measure = None,
                 is_eq = True,
                 bandlims = None):
        
        #TODO: remove this member
        self.is_eq = is_eq
        
        if samps is None:
        
            self.samps = gd.Grid()
            self.bounds = gd.Grid()
        
        else:
        
            self.samps = samps
            self.bounds = bounds
        
        self.ns = self.samps.get_curr_cards()
        
        if volume is None:
            
            self.volume = vol.Volume()
        
        else:
            
            self.volume = volume
        
        if measure is None:
            
            self.measure = ms.Measure()
        
        else:
        
            self.measure = measure
             
        
        if bandlims is None:
            
            bounds = deep(self.bounds)
            
            self.bandlims = deep(bounds.get_lims() )
            
        else:
            
            self.bandlims = deep(bandlims)
            
        return
    
    
    #TODO: remove these (X,n) deprecated getters
    def X(self, did):
        
        return self.bandlims.get_obj( (did, did) )[1]
    
    
    def n(self, did):
        
        return self.ns.get_obj( (did, did) )
    
    
    def __iter__(self, order='F', which='bs'):
        
        res_it = None
        
        res_it = self.samps.__iter__(order, which)
        
        return res_it


    def get_start(self):
        
        return self.bounds.get_start()
    
    
    def get_end(self):
        
        return self.bounds.get_end()
    
    
    def volume_at(self, inds):
        
        bounds = self.bounds
        
        res = deep( self.volume.at(bounds, inds) )
        
        return res
    
    
    def measure_at(self, inds):
        
        bounds = deep(self.bounds)
        
        for end in bounds.get_end():
        
            bounds.set_curr(0, [0, end[1] + 1], end[0][0] )
        
        vol = self.volume
        
        res = deep( self.measure.measure_at(vol, bounds, inds) )
        
        return res
    
    
    def sample_at(self, inds):
        
        res = deep( self.samps.point(inds) )
        
        return res
    
    
    def bounds_at(self, inds):
        
        res = deep( self.bounds.point(inds) )
        
        return res
    
    
    def prod(self, rhs):    
        
        #1. prod. Grids
        res_samps = self.samps.prod(rhs.samps)
        res_bounds = self.bounds.prod(rhs.bounds)
        
        #2. prod. funcs
        res_volumes = self.volume.prod(rhs.volume)
        res_measure = self.measure.prod(rhs.measure)
        
        bandlims = self.bandlims.prod(rhs.bandlims)
        
        res = Domain(res_samps,
                     res_bounds,
                     res_volumes,
                     res_measure,
                     True,
                     bandlims)
        
        return res
   
    
    def replace(self, rhs, dims, src_offs):

        self.ns.replace(rhs.ns, dims, src_offs)
        self.bandlims.replace(rhs.bandlims, dims, src_offs)
        
        #1.Use Grid.replace
        self.samps.replace(rhs.samps, dims, src_offs)
        self.bounds.replace(rhs.bounds, dims, src_offs)
        
        #2.Use Volume.replace and Measure.replace
        self.volume.replace(rhs.volume, dims, src_offs)
        self.measure.replace(rhs.measure, dims, src_offs)
        
        return
    

    def restrict(self, dims):

        bandlims = self.bandlims.restrict(dims)
        
        samps = self.samps.restrict(dims)
        bounds = self.bounds.restrict(dims)
        
        volume = self.volume.restrict(dims)
        measure = self.measure.restrict(dims)
        
        res = Domain(samps, bounds, volume, measure, True, bandlims)
        
        return res
        

    def adapt_scale(self, band_old, n_old, band_new, n_new, did):
        
        dims = (did, did)
        
        X = band_new[1] - band_new[0]
        dx = X / n_new
        shf = band_new[0]
        
        par_lb = cd.LbFixParams( [n_new], [n_new * 4 - 1] )
        par_pt = cd.PtFixParams( [dx], [shf], [X] )
        
        lb_obj = cd.LbGen(par_lb)        
        pt_obj = cd.PtGen(par_pt)
        
        mst_lb = mst.MdimStruct( [lb_obj], [(0,0)] ) 
        mst_pt = mst.MdimStruct( [pt_obj], [(0,0)] )
        
        gen = gn.Generator(mst_lb, mst_pt)      
        rgrid = gd.Grid(gen)
        
        self.samps.replace(rgrid, dims)
        self.bounds.replace(rgrid, dims)
        
        return
        

    def update(self, n_new, band_new, did):
  
        dims = (did, did)
        
        band_old = self.bandlims.get_obj(dims)
        n_old = self.ns.get_obj(dims)
        
        self.ns.replace(mst.MdimStruct( [n_new], [ (0,0) ] ), dims, 0)
        self.bandlims.replace(mst.MdimStruct( [band_new], [ (0,0) ] ), dims, 0)
        
        self.adapt_scale(band_old, n_old, band_new, n_new, did)
        
        return


    def get_mesh(self, seldims = None, md = None):
         
        #1. Extract selected 1d subdomains from domain.
        subdomains = [ self.restrict( seldims[0] ), self.restrict( seldims[1] ) ] 
        
        #2. Convert subdomains to numpy arrays.
        dmsp_x = int(subdomains[0].ns.get_obj( (0,0) ) )
        dmsp_y = int(subdomains[1].ns.get_obj( (0,0) ) )
        
        if 'color' == md:
        
            dmsp_x += 1
            dmsp_y += 1
        
        #Create the numpy array that fits dims in limits        
        subdm_x = np.zeros(shape = dmsp_x, dtype = np.float32)
        subdm_y = np.zeros(shape = dmsp_y, dtype = np.float32)

        dmls_x = []
        dmls_y = []
       
        if 'color' == md:
            
            for end in subdomains[0].bounds.get_end():
        
                subdomains[0].bounds.set_curr(0, [0, end[1] + 1], end[0][0] )
                
            for end in subdomains[1].bounds.get_end():
        
                subdomains[1].bounds.set_curr(0, [0, end[1] + 1], end[0][0] )
        
        #Fill the array with domain point coordinates
        with np.nditer(subdm_x, flags = ['multi_index'], op_flags=['writeonly'], order='C') as it:
            for elm in it:
                inds = list(it.multi_index)
                if 'color' == md:
                    dmls_x.append(subdomains[0].bounds_at(inds) )
                else:
                    dmls_x.append(subdomains[0].sample_at(inds) )
                
        with np.nditer(subdm_y, flags = ['multi_index'], op_flags=['writeonly'], order='C') as it:
            for elm in it:
                inds = list(it.multi_index)
                if 'color' == md:
                    dmls_y.append(subdomains[1].bounds_at(inds) )
                else:
                    dmls_y.append(subdomains[1].sample_at(inds) )
        
        subdm_x = np.asarray(dmls_x)
        subdm_y = np.asarray(dmls_y)
        
        #3. Create meshgrid from arrays
        res = np.meshgrid(subdm_x, subdm_y, indexing='ij')
        
        #4. Return meshgrids
        return res
    
    
    def get_polar_mesh(self, dims = None, md = None):
        
        res = None
        
        subdomain = self.restrict(dims)
        
        dmsp_r = subdomain.ns.get_obj( (0,0) )
        dmsp_a = subdomain.ns.get_obj( (1,1) )
        
        if 'color' == md:
            
            dmsp_r = int(dmsp_r + 1)                
            dmsp_a = int(dmsp_a + 1)
        
        else:
        
            dmsp_r = int(dmsp_r)                
            dmsp_a = int(dmsp_a)
        
        #Create the numpy array that fits dims in limits
        subdm_r = np.zeros(shape = (dmsp_r), dtype = np.float32)
        subdm_a = np.zeros(shape = (dmsp_a), dtype = np.float32)
        
        #go through indices and fill arrays
        for ind in range(0, dmsp_r):
            
            dr = subdomain.bandlims.get_obj( (0,0) )[1] / subdomain.ns.get_obj( (0,0) )
            
            subdm_r[ind] = ind * dr
            
        for ind in range(0, dmsp_a):
            
            da = subdomain.bandlims.get_obj( (1,1) )[1] / subdomain.ns.get_obj( (1,1) )
            
            subdm_a[ind] = ind * da
        
        res = np.meshgrid(subdm_a, subdm_r)
        
        return res