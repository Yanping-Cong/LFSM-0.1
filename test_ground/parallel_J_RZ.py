#!/usr/bin/env python
# coding: utf-8

import scipy
import h5py
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
import healpy as hp

import numpy as np
from numpy import sin,cos,pi
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
import scipy.constants as C
import healpy as hp
import h5py
import scipy.optimize as optimize
from scipy.integrate import quad

#from matplotlib import cm
from pylab import cm
import time

#python wrapping fortran code about ne2001 model 
import pyne2001

#here produce the hangqizhi diffuse sky map kelvin value after smooth
# import diffuse map from diffuse.hdf5 produced by index_ssm.py by huangqz

#read catalog
from caput import mpiutil

import sys
sys.path.append("../..")


#from Smooth.save_fit_params import free_free





import ctypes as ct
import numpy as np

# import the dll
libNE2001 = ct.CDLL('/public/home/wufq/congyanping/Software/NE2001_4python/src.NE2001/libNE2001.so')
# max integrated distance
dist = 50.


class absorption_JRZ(object):
    
    def __init__(self, v, nside, clumping_factor):
        self.v = v
        self.nside = nside
        self.clumping_factor = clumping_factor

    def pyne2001_optical_deepth(self, l, b, Te = 7000, v = 1e6):
        
        rad=57.2957795
        #radian per degree
        #distance equals 50kpc
        dist=50.0
        step = 0.01
        N =np.int(dist/step)

        print 'N',N

        nd = ct.pointer( ct.c_int(N) )          # setup the pointer

        em1D = np.arange(0, N, dtype=np.float32)  # setup the N-long

        l = l / rad #now its radian unit
        b = b / rad
        _ = libNE2001.dmdsm1_(nd, ct.pointer( ct.c_float(l) ), ct.pointer( ct.c_float(b) ), ct.pointer( ct.c_float(dist) ), np.ctypeslib.as_ctypes(em1D))
        #EM = pyne2001.get_dm_full(l, b, r)['EM']	
        Tao_mw = 3.28*1e-7 * (Te/1e4)**-1.35 * (v * 1e-9)**-2.1 * em1D
        return Tao_mw 

    def integrate_by_hand(self, f, a, b, args = [], dx=0.01):
        tao = self.clumping_factor * self.pyne2001_optical_deepth(args[0], args[1])
        
        step = 0.01

        i = a 
        s = 0
        while i <= b:
            index_ = np.int(i / step - 1)
            s += f(i,args[0],args[1],args[2],args[3]) * np.exp(-tao[index_]) * dx
            i += dx
        return s

    def split_array(self, container, count):
        #return [container[_i::count] for _i in range(count)]
        return np.split(container, count)

    def gaussian(self, x, mu = 8.5, sigma = 1.33333):
        f =  1./np.sqrt(2*np.pi*sigma**2)* np.exp(-(x-mu)**2 / (2*sigma**2))
        return f


    def _new(self, r, l, b, delt_m, params):
        param = params
        A_v = param[0]
        R_0 = param[1]
        alpha = param[2]
        R_1 = param[3]
        beta = param[4]
        Z_0 = param[5]
        gamma = param[6]
        r0 = 8.5 

        l_rad = l * np.pi/180.
        b_rad = b * np.pi/180.

        x = r * np.sin(np.pi/2. - b_rad) * np.cos(l_rad)
        y = r * np.sin(np.pi/2. - b_rad) * np.sin(l_rad)
        z = r * np.cos(np.pi/2. - b_rad)

        x_1 = x - 8.5
        y_1 = y
        z_1 = z

        r_1 = np.sqrt(np.square(x_1) + np.square(y_1) + np.square(z_1))
        b_1 = np.pi/2.0 - np.arccos(z_1/r_1)
        l_1 = np.arctan(y_1/x_1)

        R = r_1
        Z = r_1 * np.sin(b_1)
        
        #ne = (R/(R_0+0.1))**alpha * a * np.exp(-np.abs(Z) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
        emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma) 
        #tao = self.clumping_factor * self.pyne2001_optical_deepth(r, l, b)
        #step = 0.01
        #index_ = np.int(r / step)
        
        j_RZ = (emissivity + delt_m/dist) #* np.exp(-tao[index])
        
        #if (b >-4 and b<4):
        #    if (l>0 and l<45.1) or (l>330. and l<360):
        #        j_RZ = (g * ne + delt_m *gaussian(r)) * np.exp(-tao)
         
        return j_RZ


    def mpi(self):
        rank = mpiutil.rank
        size = mpiutil.size

        if rank == 0:
            
            #g = free_free(v = self.v, nside = self.nside)
            #delt_m, params = g.delta_m()
            delt_m = np.ones(np.int(12*self.nside**2))
            params = [1,2,3,4,5,6,7]
            
        else:
            delt_m = None
        #local_delt_m = mpiutil.mpilist(delt_m, method = 'con',comm = MPI.COMM_WORLD)
        local_range = mpiutil.mpirange(0,hp.nside2npix(self.nside))

        delt_m = mpiutil.bcast(delt_m, root = 0)
        params = mpiutil.bcast(params, root = 0)
        result = []
        for pix_number in local_range:
            a = time.time()
            l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
            #pix_value = delt_m[pix_number] + 100.
            pix_value =self.integrate_by_hand(self._new, 0.01, dist, args=(l, b, delt_m[pix_number], params)) 
            b = time.time()
            
            print 'pix_number', pix_number, 'delta_time', b-a
            result.append([pix_number, pix_value])
        result = mpiutil.gather_list(result, root = None)
        if rank == 0:
            with h5py.File('text_0to1.hdf5','w') as f:
                f.create_dataset('result', data = result)
                print 'end'

if __name__ == '__main__':
    cla = absorption_JRZ(v = 1., nside = 2**5, clumping_factor = 6.74)
    cla.mpi()

