#!/usr/bin/env python
# coding: utf-8
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
from numpy import sin,cos,pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
import scipy.constants as C
import healpy as hp
import h5py
import scipy.optimize as optimize
from scipy.integrate import quad

from pylab import cm
import time

from LFSM.fitting_params.how_to_smooth_SSM import smooth
from LFSM.I_E_term.I_E_equation import I_E
from LFSM.interpolate_sky.interpolate_sky_map import produce_index
class free_free(object):

    def __init__(self, v, nside, index_type, dist,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,only_fit_Anu):
        self.v = v#20Mhz frequency in hz
        self.nside = nside       
        self.index_type = index_type
        self.dist = dist
        self.I_E_form = I_E_form
        self.emi_form = emi_form
        self.R0_R1_equal = R0_R1_equal
        self.using_raw_diffuse = using_raw_diffuse
        self.only_fit_Anu = only_fit_Anu
        self.params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77])

    def plot_mollview(self,total,filename = '',log=True):
        cool_cmap = cm.jet
        cool_cmap.set_under("w") # sets background to white
        plt.figure(1)
        if log==True:
            m = np.log10(total)
            Min = None
            Max = None
            hp.mollview(m,title="The frequency in"+' '+str(self.v)+'MHz', min = Min, max = Max,cmap = cool_cmap)
        if log==False:
            m = np.log10(total)
            Min = None
            Max = None
            hp.mollview(total,title="The frequency in"+' '+str(self.v)+'MHz', min = Min, max = Max,cmap = cool_cmap)
        plt.savefig(str(self.v)+'MHz'+ filename +'.eps',format = 'eps')
        return 
    
    
    def produce_xyz(self):
        #v in Mhz
        result = smooth(self.nside,self.v, self.index_type,self.I_E_form,self.using_raw_diffuse).add_5()
        return result

    def diffuse_raw(self):
        diffuse_x = smooth(self.nside, self.v, self.index_type,self.I_E_form,self.using_raw_diffuse).produce_data()
        return diffuse_x
 
    def fun_for_curvefit(self, xyz, A_v, R_0, alpha, R_1, beta, Z_0, gamma, form = 'n_HII'):
        #produce for input for curve_fit
        result = []
        for l,b in xyz:
            self.l = l * np.pi / 180.0
            self.b = b * np.pi / 180.0

            def fun_inte(r):
                R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
                Z = r * np.sin(self.b)
                #integrate along the sight direction
                #return h * np.square(a * np.exp(-np.abs(Z)/b - np.square(R/c)) + d * np.exp(-np.abs(Z)/e - np.square(R/f - g)))
                #return a * (R/b)**c * np.exp(-d*(R-e/e) - np.abs(Z)/f)
                #integrate along sight direction
                #return e * np.square(a * np.exp(-R/b) * np.square(2/(np.exp((Z-d)/c)+ np.exp(-(Z-d)/c))))
                #ne = f * np.exp(-np.abs(d)/1e3 - (r_1/(2e4*A +1))**2) + g * np.exp(-np.abs(d)/(0.15*1e3) - (r_1/(2e3*B+1) - 2)**2)
                #ne = f * np.exp(-np.abs(d)/(1e3*B+1) - (r_1/(2e4*A +1))**2) 
                #ne = f * np.exp(-np.abs(d)/(1e3*B+1) - (r_1/(2e4*A +1))**4) 
                #ne = np.square(a * np.exp(-R/b) * np.square(2.0/(np.exp(Z/(1e3*c+1))+ np.exp(-Z/(1e3*c+1)))))
            	#ne = a * np.exp(-np.abs(d) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D * np.exp(-np.abs(d)*2/(e * 0.15 + 0.01) - 2*(r_1/(2*f+0.1))**2)
            	#ne = a * np.exp(-np.abs(d) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
            	#ne = (R/(R_0+0.1))**alpha * a * np.exp(-np.abs(Z) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
                emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                #get rid of square 
                return emissivity

            result.append(quad(fun_inte, 0, self.dist)[0])
        return np.array(result)

    def fun_for_curvefit_only_R0(self, xyz, A_v, R_0, alpha, Z_0, gamma, form = 'n_HII'):
        #produce for input for curve_fit
        beta = 1
        result = []
        for l,b in xyz:
            self.l = l * np.pi / 180.0
            self.b = b * np.pi / 180.0

            def fun_inte(r):
                #x = r * np.sin(np.pi/2. - self.b) * np.cos(self.l)
                #y = r * np.sin(np.pi/2. - self.b) * np.sin(self.l)
                #z = r * np.cos(np.pi/2. - self.b)

                #x_1 = x - 8.5
                #y_1 = y
                #z_1 = z

                #r_1 = np.sqrt(np.square(x_1) + np.square(y_1) + np.square(z_1))
                #b_1 = np.pi/2.0 - np.arccos(z_1/r_1)
                #l_1 = np.arctan(y_1/x_1)

                ##R = r_1
                #R = np.sqrt(r_1**2 - z**2)
                #Z = r_1 * np.sin(b_1)
                R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
                Z = r * np.sin(self.b)
                emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_0)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                #get rid of square 
                return emissivity

            result.append(quad(fun_inte, 0, self.dist)[0])
        return np.array(result)

    def fun_for_curvefit_R0R2(self, xyz, A_v, R_0, alpha, Z_0, gamma, form = 'n_HII',R_2 = 0.1):
        #produce for input for curve_fit
        beta = 1
        R_2 = R_2
        result = []
        for l,b in xyz:
            self.l = l * np.pi / 180.0
            self.b = b * np.pi / 180.0

            def fun_inte(r):
                R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
                Z = r * np.sin(self.b)
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_0)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                #get rid of square 
                return emissivity

            result.append(quad(fun_inte, 0, self.dist)[0])
        return np.array(result)
    
    def fun_for_curvefit_only_fit_Anu(self, xyz, A_v):
        #produce for input for curve_fit relatively error
        ##params = np.array([1.66983971e+06, 5.50771193e+00, -4.82064326e-02, 6.42470460e-01, 5.20197943e-01])
        #fitting params using absolute error
        ##params = np.array([1.70973109e+08, 3.13073265e+00, 6.75317432e-01, 8.80354845e-01, 1.10987383e+00])
        params = self.params_408
        R_0 = params[1]
        alpha = params[2]
        Z_0 = params[3]
        gamma = params[4]

        beta = 1
        R_2 = 0.1
        result = []
        for l,b in xyz:
            self.l = l * np.pi / 180.0
            self.b = b * np.pi / 180.0

            def fun_inte(r):
                R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
                Z = r * np.sin(self.b)
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_0)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                #get rid of square 
                return emissivity

            result.append(quad(fun_inte, 0, self.dist)[0])
        return np.array(result)
    def curve_fit_for_Anu(self):
        #A_v
        guess = [1]

        func = self.fun_for_curvefit_only_fit_Anu
        xyz = self.produce_xyz()
        #print 'xyz.shape',xyz.shape
        # bug report using absolute error
        result, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0]),np.array([1e10])), method='trf')

        ##params = np.array([1.66983971e+06, 5.50771193e+00, -4.82064326e-02, 6.42470460e-01, 5.20197943e-01])
        ##params = np.array([1.70973109e+08, 3.13073265e+00, 6.75317432e-01, 8.80354845e-01, 1.10987383e+00])
        params = self.params_408
        params[0] = result[0]

        with h5py.File(str(self.v)+'Mhz_fitted_param.hdf5','w') as f:
            f.create_dataset('params',data = params)
            f.create_dataset('v',data = self.v)
            f.create_dataset('pcov', data = pcov)
        #print 'frequency',self.v
        #print 'params',params
        #print 'pcov',pcov
        return params

    def model_m3(self,l,b,abcz0,form='two_componant',R_2=0.1):

        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = abcz0
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1

        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
            #get rid of square 
            return emissivity 

        return quad(fun_inte, 0, self.dist)[0]

    
    def curve_fit(self,form='CRs'):
        #if self.R0_R1_equal == False:
        if False:
            #R_0, alpha, a,B,c,D,g
            #A_v, R_0, alpha, R_1, beta, Z_0, gamma
            guess = [100,8.5,2,4,2,3,1]

            func = self.fun_for_curvefit
            xyz = self.produce_xyz()
            #print 'xyz.shape',xyz.shape
            params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess,  bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,100,3.1,20,3.1])), method='trf')


       # if self.R0_R1_equal == True:
       #     #A_v, R_0, alpha, Z_0, gamma
       #     guess = [1e7,5,2,1.6,1]

       #     func = self.fun_for_curvefit_only_R0
       #     xyz = self.produce_xyz()
       #     print 'xyz.shape',xyz.shape
       #     params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess,sigma=xyz[:,2],  bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,20,3.1])), method='trf')

        #if self.R0_R1_equal == True:
        if True:
            #A_v, R_0,R_2=0.1, alpha, Z_0, gamma
            guess = [80,5,2,1.6,1]

            func = self.fun_for_curvefit_R0R2
            xyz = self.produce_xyz()
            #print 'xyz.shape',xyz.shape
            #params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess,sigma = xyz[:,2], bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,20,3.1])), method='trf')
            #params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,20,3.1])), method='trf')
            params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0,1e-5,1e-5,1e-5,1e-5]),np.array([150,5,3.1,2,3.1])), method='trf')

        with h5py.File(str(self.v)+'Mhz_fitted_param.hdf5','w') as f:
            f.create_dataset('params',data = params)
            f.create_dataset('v',data = self.v)
            f.create_dataset('pcov', data = pcov)
        #print 'frequency',self.v
        #print 'params',params
        #print 'pcov',pcov
        return params

    def model_m2(self,l,b,abcz0,form='two_componant',R_2=0.1):

        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        if self.R0_R1_equal == False:

            A_v, R_0, alpha, R_1, beta, Z_0, gamma = abcz0
            R_2 = 0
        if self.R0_R1_equal == True:
            #A_v, R_0, alpha, Z_0, gamma = abcz0
            #R_1 = R_0.copy()
            #beta = 1
            A_v, R_0,alpha, Z_0, gamma = abcz0
            R_1 = R_0.copy()
            R_2 = R_2
            beta = 1
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
            #get rid of square 
            return emissivity 

        return quad(fun_inte, 0, self.dist)[0]

    def model_m4(self,l,b,params_408,pix_number,R_2 = 0.1):
        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = params_408
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1
        #pix_number = hp.ang2pix(self.nside, l, b, lonlat = True)
        f = produce_index(Nside = self.nside, freq = self.v, index_type = self.index_type,I_E_form = self.I_E_form)
        Beta_G = f.pixel_dependence_index_minus_I_E()
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**Beta_G[pix_number]
            #get rid of square 
            return emissivity 

        return quad(fun_inte, 0, self.dist)[0]
    def model_m5(self,l,b,params_408,R_2 = 0.1):
        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = params_408
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1
        #pix_number = hp.ang2pix(self.nside, l, b, lonlat = True)
        f = produce_index(Nside = self.nside, freq = self.v, index_type = self.index_type,I_E_form = self.I_E_form)
        Beta_G = f.constant_index_minus_I_E()
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**Beta_G[0]
            #get rid of square 
            return emissivity 
        return quad(fun_inte, 0, self.dist)[0]

    def I_E(self, v):
        f = I_E(v,self.I_E_form)
        result = f.I_E()
        #print 'I_E result',result
        return result
 
    def delta_m(self):
        try:
            with h5py.File('./'+str(self.v)+'Mhz_fitted_param.hdf5','r') as f:
                abcz0 = f['params'][:]
            print ('this is using fitted params')
            self.produce_xyz()
        except:
            print ('i am here doing fitting now')
            if self.index_type == 'pixel_dependence_index_minus_I_E':
               #params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77])
               abcz0 = self.params_408
            elif self.index_type == 'constant_index_minus_I_E':
                if int(self.v) == int(408):
                    try:
                        with h5py.File('408Mhz_fitted_param.hdf5','r') as f:
                            abcz0 = f['params'][:]
                    except:
                        abcz0 = self.curve_fit()
                    abcz0 = self.params_408
                else:
                    abcz0 = self.params_408
            else:
                if self.only_fit_Anu == True:
                    abcz0 = self.curve_fit_for_Anu()
                else:
                    abcz0 = self.curve_fit()
        #print 'params',abcz0 
        nside = self.nside
        m = np.zeros(hp.nside2npix(nside))

        I_E = self.I_E(self.v)
        for pix_number in range(0,hp.nside2npix(nside)):
            l,b = hp.pixelfunc.pix2ang(nside, pix_number, nest = False, lonlat = True)
            #the emissivity after absorption and plus the term of extragalaxy
            a = time.time()
            if self.index_type == 'pixel_dependence_index_minus_I_E':
                pix_value = self.model_m4(l,b,abcz0,pix_number) + I_E
            elif self.index_type == 'constant_index_minus_I_E':
                if int(self.v) == int(408):
                    pix_value = self.model_m2(l,b,abcz0) + I_E
                else:
                    pix_value = self.model_m5(l,b,abcz0) + I_E
                    
            else: 
                if self.only_fit_Anu == True:
            
                    pix_value = self.model_m3(l,b,abcz0) + I_E
                else:
                    pix_value = self.model_m2(l,b,abcz0) + I_E
            m[pix_number] = pix_value
            

            b = time.time()
            #print 'delt_m delay',b-a
        #print 'produce delt_m and params is over'
        #self.plot_mollview(m,filename= 'integrated_temperature')

        diffuse_raw = self.diffuse_raw() 
        delt_m = diffuse_raw + I_E - m

 
        #residual between fitted emissivity and raw input not the smoothed raw input 
        I_E_value = self.I_E(self.v)
        delt_m_percentage = delt_m / diffuse_raw * 100 
        #self.plot_mollview(delt_m,filename = 'delt_m',log=False)
        #self.plot_mollview(delt_m_percentage,filename = 'delt_m_percentage',log=False)
        
        with h5py.File('./'+str(self.emi_form)+str(self.v)+'Mhz' + '_delt_m_and_unabsorb_and_delt_m_percentage.hdf5','w') as f:
            f.create_dataset('delt_m',data = delt_m)
            f.create_dataset('delt_m_percentage',data = delt_m_percentage)
            f.create_dataset('integrated_temperature_total_m', data = m)
            f.create_dataset('diffuse_raw', data = diffuse_raw) 
            f.create_dataset('I_E_value', data = I_E_value) 
        #back to delt_m and params     
        return delt_m, abcz0 



