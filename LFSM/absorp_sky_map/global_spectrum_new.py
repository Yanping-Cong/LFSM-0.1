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
import matplotlib.pyplot as plt
import scipy.constants as C
import healpy as hp
import h5py
import scipy.optimize as optimize
from scipy.integrate import quad

#from matplotlib import cm
from pylab import cm
import time

#python wrapping fortran code about ne2001 model 
#import pyne2001

#here produce the hangqizhi diffuse sky map kelvin value after smooth
# import diffuse map from diffuse.hdf5 produced by index_ssm.py by huangqz

#read catalog
from caput import mpiutil



from LFSM.Smooth.save_fit_params import free_free
#from Smooth.least_sq_fit_params import free_free
#import "./F2py_file" 
from LFSM.I_E_term.I_E_equation import I_E
from LFSM.Index.produce_index import produce_index


import ctypes as ct
import numpy as np
import os
_path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append(_path)
_path = os.path.split(_path)[0]
file_dir = _path
print ('_path',file_dir + '/NE2001_4python/src_NE2001/libNE2001.so')
# import the dll
libNE2001 = ct.CDLL('/public/home/wufq/congyanping/Software/NE2001_4python/src.NE2001/libNE2001.so')
# max integrated distance (kpc)
dist = 50.


class absorption_JRZ(object):
    
    def __init__(self, v, nside, clumping_factor, index_type, distance,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,test, only_fit_Anu):
        self.v = v
        self.nside = nside
        self.clumping_factor = clumping_factor
        self.test = test
        self.index_type = index_type
        self.dist = distance
        self.emi_form = emi_form
        self.I_E_form = I_E_form
        self.R0_R1_equal = R0_R1_equal
        self.using_raw_diffuse = using_raw_diffuse
        self.only_fit_Anu = only_fit_Anu
        Te = 8000 
 
        f = produce_index(Nside = self.nside, freq = self.v, index_type = self.index_type, I_E_form = self.I_E_form)
        self.Beta_G = f.pixel_dependence_index_minus_I_E()
        f = produce_index(Nside = self.nside, freq = self.v, index_type = self.index_type, I_E_form = self.I_E_form)
        self.Beta_G_constant = f.constant_index_minus_I_E()

    def Fortran2Py_optical_deepth(self, l, b, Te = 8000):
        v = self.v * 1e6 #v in MHz
        rad=57.2957795
        #radian per degree
        #distance equals 50kpc
        #dist=50.0
        
        if self.test == True:
            
            step = 0.1
        else:
            step = 0.01
        N =np.int(dist/step)

        #print 'N',N

        nd = ct.pointer( ct.c_int(N) )          # setup the pointer

        em1D = np.arange(0, N, dtype=np.float32)  # setup the N-long

        l_rad = l / rad #now its radian unit
        b_rad = b / rad
        _ = libNE2001.dmdsm1_(nd, ct.pointer( ct.c_float(l_rad) ), ct.pointer( ct.c_float(b_rad) ), ct.pointer( ct.c_float(dist) ), np.ctypeslib.as_ctypes(em1D))
        #EM = pyne2001.get_dm_full(l, b, r)['EM']	
        Tao_mw = 3.28*1e-7 * (Te/1e4)**-1.35 * (v * 1e-9)**-2.1 * em1D
        #print 'Tao_mw',Tao_mw 
        return Tao_mw
 
    def raw_pyne2001_optical_deepth(self, r, l, b, Te = 8000):
        v = self.v * 1e6
        EM = pyne2001.get_dm_full(l, b, r)['EM']	
        Tao_mw = 3.28*1e-7 * (Te/1e4)**-1.35 * (v * 1e-9)**-2.1 * EM
        return Tao_mw 
         
    def integrate_by_hand(self, f, a, b, args = [], dx=0.01):
        if self.test == True:
            dx = 0.1
            step = dx
        else:
            dx = 0.01
            step = dx

        tao = self.Fortran2Py_optical_deepth(args[0], args[1])
        

        i = a 
        s = 0

        ##I_E = args[3][-1]
        I_E = self.I_E(self.v)

        while i <= b:
            index_ = np.int(i / step - 1)
            s += (f(i,args[0],args[1],args[2],args[3]) * np.exp(-tao[index_])) * dx
            i += dx
        #here find the bug
        s = s + I_E*np.exp(-tao[-1])
        return s

    def Quad(self, f, a, b, args = [], dx=0.01):
        #the different to integrate_by_hand is not including I_E
        if self.test == True:
            dx = 0.1
            step = dx
        else:
            dx = 0.01
            step = dx

        tao = self.Fortran2Py_optical_deepth(args[0], args[1])
        

        i = a 
        s = 0


        while i <= b:
            index_ = np.int(i / step - 1)
            s += (f(i,args[0],args[1],args[2],args[3]) * np.exp(-tao[index_])) * dx
            i += dx
        #here find the bug
        s = s 
        return s


    def integrate_by_hand_unabsorb(self, f, a, b, args = [], dx=0.01):

        i = a 
        s = 0
        while i <= b:
            s += f(i,args[0],args[1],args[2],args[3]) * dx
            i += dx
        return s

    def integrate_by_hand_low_resolution(self, f, a, b, args = [], dx=0.1):

        i = a 
        s = 0
        while i <= b:
            s += f(i,args[0],args[1],args[2],args[3]) * dx
            i += dx
        return s

    def split_array(self, container, count):
        #return [container[_i::count] for _i in range(count)]
        return np.split(container, count)

    def gaussian(self, x, mu = 8.5, sigma = 1.33333):
        f =  1./np.sqrt(2*np.pi*sigma**2)* np.exp(-(x-mu)**2 / (2*sigma**2))
        return f

    def sech2(self,x):
        return np.square(2/(np.exp(x) + np.exp(-x)))

    def I_E(self, v):
        f = I_E(v,self.I_E_form)
        result = f.I_E()
        return result

    def _new(self, r, l, b, delt_m, params):
        if self.R0_R1_equal == True:

            param = params
            A_v = param[0]
            R_0 = param[1]
            R_2 = 0.1
            alpha = param[2]
            R_1 = param[1]
            #beta = param[3]
            beta = 1
            Z_0 = param[3]
            gamma = param[4]

        if self.R0_R1_equal == False:
            param = params
            A_v = param[0]
            R_0 = param[1]
            alpha = param[2]
            R_1 = param[3]
            beta = param[4]
            Z_0 = param[5]
            gamma = param[6]

        if self.only_fit_Anu == True:
            param = params
            A_v = param[0]
            R_0 = param[1]
            R_2 = 0.1
            alpha = param[2]
            R_1 = param[1]
            #beta = param[3]
            beta = 1
            Z_0 = param[3]
            gamma = param[4]

        #I_E = param[7]

        r0 = 8.5 

        l_rad = l * np.pi/180.
        b_rad = b * np.pi/180.
        
        """
        x = r * np.sin(np.pi/2. - b_rad) * np.cos(l_rad)
        y = r * np.sin(np.pi/2. - b_rad) * np.sin(l_rad)
        z = r * np.cos(np.pi/2. - b_rad)

        x_1 = x - 8.5
        y_1 = y
        z_1 = z

        r_1 = np.sqrt(np.square(x_1) + np.square(y_1) + np.square(z_1))
        b_1 = np.pi/2.0 - np.arccos(z_1/r_1)
        l_1 = np.arctan(y_1/x_1)

        #R = r_1
        R = np.sqrt(r_1**2 - z**2) 
        Z = r_1 * np.sin(b_1)
        """
        R = np.sqrt(8.5**2 + (r*np.cos(b_rad))**2 -2*8.5*(r*np.cos(b_rad))*np.cos(l_rad))
        Z = r * np.sin(b_rad)

        ########ne = (R/(R_0+0.1))**alpha * a * np.exp(-np.abs(Z) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
        #emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
        if self.index_type == 'pixel_dependence_index_minus_I_E':
            pix_number = hp.ang2pix(self.nside, l, b, lonlat = True)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G[pix_number]
        elif self.index_type == 'constant_index_minus_I_E':
            if int(self.v) == int(408):
                if self.emi_form == 'exp':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                if self.emi_form == 'sech2':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)
            else:
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G_constant[0]
             
        else:
            if self.emi_form == 'exp':

                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
            if self.emi_form == 'sech2':
                    
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)
        j_RZ = emissivity #+ delt_m/dist) #* np.exp(-tao[index])
        return j_RZ

    def critical_distance(self,l,b,delt_m,params):
        import scipy.optimize as so
        #import scipy.integrate as integrate
        #bug report : the lower limit is from 0.01 not 0
        value = 0.5 * self.Quad(self._new, 0.01, 50,args=(l,b,delt_m,params)) 
        def func(x,l,b,delt_m,params):
            return self.Quad(self._new, 0.01, x,args=(l,b,delt_m,params)) - value
        #sol = so.fsolve(func,np.array([1]),args=(l,b,delt_m,params),xtol=1,maxfev=1000)
        sol = 0
        Y = []
        for i in np.arange(0.01,50,0.01):
            result = self.Quad(self._new,0.01,i,args=(l,b,delt_m,params)) - value
            Y.append(result)
        Y = list(np.abs(Y))
        container = np.arange(0.01,50,0.01)
        index = Y.index(min(Y))
        sol = container[index]  
            #if np.abs(result) < 100:
            #    sol = i
            #    break
        #print 'begin_crital', func(sol,l,b,delt_m,params),'end_critical','sol',sol,'min',min(np.abs(Y)),'index',index
        return sol

    def _new_unabsorb(self, r, l, b, delt_m, params):
        param = params
        A_v = param[0]
        R_0 = param[1]
        alpha = param[2]
        R_1 = param[3]
        beta = param[4]
        Z_0 = param[5]
        gamma = param[6]
        I_E = param[7]
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
        
        ########ne = (R/(R_0+0.1))**alpha * a * np.exp(-np.abs(Z) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
        #emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
        if self.emi_form == 'exp':

            emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma) + I_E
        if self.emi_form == 'sech2':
                    
            emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma) + I_E
        j_RZ = emissivity #+ delt_m/dist) #* np.exp(-tao[index])
        return j_RZ

    def raw_new_absorb(self, r, l, b, delt_m, params):
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
        
        if self.emi_form == 'exp':

            emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
        if self.emi_form == 'sech2':
                    
            emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)
        tao = self.clumping_factor * self.raw_pyne2001_optical_deepth(r, l, b)
        j_RZ = emissivity * np.exp(-tao) 
         
        return j_RZ

    def mpi(self):
        rank = mpiutil.rank
        size = mpiutil.size

        if rank == 0:
             
            g = free_free(v = self.v, nside = self.nside,index_type = self.index_type,dist = self.dist,emi_form = self.emi_form,I_E_form = self.I_E_form,R0_R1_equal = self.R0_R1_equal,using_raw_diffuse = self.using_raw_diffuse,only_fit_Anu = self.only_fit_Anu)
            delt_m, params = g.delta_m()
            
        else:
            delt_m = None
            params = None
        #local_delt_m = mpiutil.mpilist(delt_m, method = 'con',comm = MPI.COMM_WORLD)
        local_range = mpiutil.mpirange(0,hp.nside2npix(self.nside))

        delt_m = mpiutil.bcast(delt_m, root = 0)
        params = mpiutil.bcast(params, root = 0)
        result_absorb = []
        for pix_number in local_range:
            a = time.time()
            l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
            if self.test == True:
                pix_value =self.integrate_by_hand(self._new, 0.1, dist, args=(l, b, delt_m[pix_number], params)) 
            else:

                pix_value =self.integrate_by_hand(self._new, 0.01, dist, args=(l, b, delt_m[pix_number], params)) 
                distance = self.critical_distance(l,b,delt_m[pix_number],params)
                l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
            b = time.time()
            
            if self.test == True:

                result_absorb.append([pix_number, pix_value])
            else:
                result_absorb.append([pix_number, pix_value, distance])
        if self.test == True:

            result_absorb = mpiutil.gather_list(result_absorb, root = None)
        else:
            result_absorb = mpiutil.gather_list(result_absorb, root = None)
        if rank == 0:
            if self.test == True:
                with h5py.File('./' + str(self.emi_form)+str(self.v) + 'F2py_absorb.hdf5', 'w') as f:
                    f.create_dataset('F2py_absorb', data = result_absorb)
            else:
                with h5py.File('./' + 'exp'+str(self.v)+'Mhz_delt_m_and_unabsorb_and_delt_m_percentage.hdf5','r') as f:
                    #print f.keys()
                    unabsorb = f['integrated_temperature_total_m'][:]
                    diffuse_raw = f['diffuse_raw'][:]
                result_absorb = np.array(result_absorb)
                absorb = result_absorb[:,1]
                I_E = self.I_E(self.v)
                result = []
                print ('in the beginning')
                for pix_number in range(unabsorb.size):
                    print ('left number of pixel',unabsorb.size - pix_number) 
                    X = unabsorb[pix_number] - self.I_E(self.v)
                    l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
                    tao = self.Fortran2Py_optical_deepth(l, b)
                    Y = absorb[pix_number] - I_E * np.exp(-tao[-1])
                    mean_exptao = Y / X
                    pix_value = diffuse_raw[pix_number] * mean_exptao + I_E*np.exp(-tao[-1])
                    result.append([pix_number,pix_value])
                    #print 'pixel_number',pix_number
                with h5py.File('./' + str(self.emi_form)+str(self.v)+'MHz_global_spectrum_with_perterbation.hdf5','w') as h:
                    h.create_dataset('result',data = result)
                    #h.create_dataset('smooth_result',data = result_absorb)
                    print ('end, good job!, you are the best')

        return result

#if __name__ == '__main__':
#    for v in np.arange(0.1,1.0,0.1):
#        v = round(v,1)
#        nside = 2**6
#        # step integrate = 0.1, only calculate F2py absorb result for test = True
#        cla = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False,only_fit_Anu = False)
#        cla.mpi()
#        #cla.Fortran2Py_optical_deepth(l=60,b=10.)
#index_type = 'constant_index_minus_I_E'            I_E_form:"no_seiffert","seiffert","Dowell"
#             'freq_dependence_index_minus_I_E'     I_E_form:"extra_freq_dependence","seiffert_freq_depend"
#             'pixel_dependence_index_minus_I_E'
