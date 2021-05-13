# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:14:21 2021

@author: Domokos Zoltán
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_partindsfc(fix_inds):
    
    #TODO: adapt to lod and lims updates
    def fc(lod, lims, sel_inds):
        return fix_inds
    
    return fc


def plot_heatmap(signal, md = 'cartes'):
    
    if 'polar' == md:
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar') )
        
    else:
        
        fig, ax = plt.subplots()
    
    vs = signal.sigvec() 
    
    if 'polar' == md:
        
        vx, vy = signal.get_polar_mesh( (0, 1), md = 'color' )
    
    else:
        
        vx, vy = signal.get_mesh( [ (0, 0), (1, 1) ], md = 'color' )
    
    ax.pcolormesh(vx, vy, vs[1].real)
    
    plt.show()
    
    return


def plot_heatvec(domain, sigvec, md = 'cartes'):
    
    if 'polar' == md:
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar') )
        
    else:
        
        fig, ax = plt.subplots()
    
    vs = sigvec
    
    if 'polar' == md:
        
         vx, vy = domain.get_polar_mesh( (0, 1), md = 'color' )
        
    else:
        
        vx, vy = domain.get_mesh( [ (0, 0), (1, 1) ], md = 'color' )
    
    ax.pcolormesh(vx, vy, vs.real)
    
    plt.show()
    
    return


def plot_sigvec(domain, sigvec):
    
    fig = plt.figure()
    
    ax = Axes3D(fig)
    
    vs = sigvec
    
    vx, vy = domain.get_mesh( [ (0, 0), (1, 1) ] )
    
    ax.plot_surface(vx, vy, vs.real, rstride = 1, cstride = 1)
    
    plt.show()
    
    return  


def plot_sigfc(signal):
    
    fig = plt.figure()
    
    ax = Axes3D(fig)
    
    vs = signal.sigvec() 
    
    vx, vy = signal.get_mesh( [ (0, 0), (1, 1) ] )
    
    ax.plot_surface(vx, vy, vs[1].real, rstride = 1, cstride = 1)
    
    plt.show()
    
    return  


def plot_sigfc2(signal, seldims, fix_inds):
    
    fig = plt.figure()
    
    ax = Axes3D(fig)
    
    new_partindsfc = get_partindsfc(fix_inds)
    
    vs = signal.sigvec2(None, None, None, seldims, new_partindsfc)
    
    vx, vy = signal.get_mesh(seldims)
    
    ax.plot_surface(vx, vy, vs[1].real, rstride = 1, cstride = 1)
    
    plt.show()
    
    return
    

def plot_heatmap2(signal, seldims, fix_inds):
    
    fig, ax = plt.subplots()
    
    new_partindsfc = get_partindsfc(fix_inds)
    
    vs = signal.sigvec2(None, None, None, seldims, new_partindsfc)
    
    vx, vy = signal.get_mesh( seldims, md = 'color' )
    
    ax.pcolormesh(vx, vy, vs[1].real)
    
    plt.show()
    
    return

