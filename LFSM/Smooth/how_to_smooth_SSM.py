import numpy as np
import healpy as hp
import h5py
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import sys
import os
from LFSM.Index.produce_index import produce_index
class smooth(object):
    def __init__(self, nside, v, index_type,I_E_form,using_raw_diffuse):
        self.nside = nside
        self.v = v
        self.index_type = index_type
        self.I_E_form = I_E_form
        self.using_raw_diffuse = using_raw_diffuse
 
    def produce_data(self):
        
        f = produce_index(Nside = self.nside, freq = self.v, index_type = self.index_type, I_E_form = self.I_E_form)
        diffuse_x = f.diffuse_x(freq = self.v)
        return diffuse_x

    def masked_smoothing(self, U, rad=10.0):     
        V=U.copy()
        V[U!=U]=0
        VV=hp.smoothing(V, fwhm=np.radians(rad))    
        W=0*U.copy()+1
        W[U!=U]=0
        WW=hp.smoothing(W, fwhm=np.radians(rad))    
        return VV/WW

    def mask_spur_loopI(self, arr):
        
        nside = self.nside
        l1,b1 = hp.pix2ang(nside,100,lonlat = True)
        l2,b2 = hp.pix2ang(nside,101, lonlat = True)
        #print (l1,l2,b1,b2)
        delt_l = np.abs(l2-l1)
        #print (delt_l)
        mask = np.zeros_like(arr,dtype = np.bool)
        #print (mask)
        for l in np.arange(-30,45,0.01):
            for b in np.arange(20,90,delt_l * 0.01):
                pix_number = hp.ang2pix(nside,l,b,lonlat = True)
                mask[pix_number] = True
        arr[mask] = np.nan
        return arr

    #def galprop_simulate_data(self):
    #    with h5py.File('/public/home/wufq/congyanping/Software/LFSAM/galprop_data/galprop_diffuse_x.hdf5','r') as f:
    #        diffuse_x = f['diffuse_x'][:]
    #    diffuse_x = hp.ud_grade(diffuse_x,nside_out = self.nside)
    #    result = []
    #    for index in range(diffuse_x.size):
    #        ll,bb = hp.pix2ang(self.nside, index, lonlat = True)
    #        result.append([ll,bb,diffuse_x[index]*1e22])
    #    result = np.array(result)
    #    return result,diffuse_x

    def add_5(self):
        diffuse = self.produce_data()

        
        self.mask_spur_loopI(diffuse) 

        if self.using_raw_diffuse == True:
            smooth_diffuse = diffuse
        
            self.mask_spur_loopI(smooth_diffuse)

        if self.using_raw_diffuse == False:
            if self.index_type == 'pixel_dependence_index_minus_I_E':
                smooth_diffuse = self.masked_smoothing(diffuse,rad=1.0)

            if self.index_type == 'constant_index_minus_I_E':
                smooth_diffuse = self.masked_smoothing(diffuse,rad=1.0)
            if self.index_type == 'freq_dependence_index_minus_I_E':
                smooth_diffuse = self.masked_smoothing(diffuse,rad=1.0)
                
        
            self.mask_spur_loopI(smooth_diffuse)

        result = []
        for index in range(smooth_diffuse.size):
        #for index in np.where(smooth_diffuse != np.nan)[0]:
            ll,bb = hp.pix2ang(self.nside, index, lonlat = True)
            #remove small scale terbulance
            #CAS-a 111.73,-2.13
            #Cygnus 71,4
            #Cen-A 309.51589573409,19.41727341133
            #SMC 302.8084 -44.3277
            #LMC 280.4652 -32.8884
            #Vela 263.9390 -03.3683
            nside = 2**4
            a = (hp.get_all_neighbours(nside,111.73,-2.13,lonlat = True))
            b = (hp.get_all_neighbours(nside,71,4,lonlat = True))
            c = (hp.get_all_neighbours(nside,309.51589573409,19.41727341133,lonlat = True))
            d = (hp.get_all_neighbours(nside,302.8084,-44.3277,lonlat = True))
            e = (hp.get_all_neighbours(nside,280.4652,-32.8884,lonlat = True))
            f = (hp.get_all_neighbours(nside,263.9390,-3.3683,lonlat = True))
            g = (hp.get_all_neighbours(nside,0,0,lonlat = True))
            total = list(a) + list(b) + list(c)+list(d)+list(e)+list(f)+list(g)

            a = (hp.ang2pix(nside,111.73,-2.13,lonlat = True))

            b = (hp.ang2pix(nside,71,4,lonlat = True))
            c = (hp.ang2pix(nside,309.51589573409,19.41727341133,lonlat = True))
            d = (hp.ang2pix(nside,302.8084,-44.3277,lonlat = True))
            e = (hp.ang2pix(nside,280.4652,-32.8884,lonlat = True))
            f = (hp.ang2pix(nside,263.9390,-3.3683,lonlat = True))
            g = (hp.ang2pix(nside,0,0,lonlat = True))
            total2 = [a,b,c,d,e,f,g]

            pass_pix_num = np.array(total+total2) 
            mask = np.ones_like(smooth_diffuse,dtype = np.bool)
            mask[pass_pix_num] = False
            if mask[index] == False:
                print ('pass_pix_number',index)
            if mask[index] == True:
                pix_value = smooth_diffuse[index]
                if ~np.isnan(pix_value):
                    result.append([ll,bb,pix_value])
        result = np.array(result)
        with h5py.File(str(self.v)+'MHz_Smooth_'+'after_mask_smooth_diffuse_x.hdf5','w') as f:
            f.create_dataset('smooth_diffuse',data = smooth_diffuse)
            f.create_dataset('masked_diffuse',data = diffuse)
            f.create_dataset('diffuse',data = diffuse)
        print (result[-1,:]) 
        return result

    def raw_and_smooth(self):
        diffuse = self.produce_data()
        self.mask_spur_loopI(diffuse) 
        smooth_diffuse = self.masked_smoothing(diffuse)
        self.mask_spur_loopI(smooth_diffuse)

        with h5py.File(str(self.v)+'Mhz_diffuse_x_and_smooth_diffuse_x.hdf5','w') as f:
            f.create_dataset('diffuse_x',data = diffuse)
            f.create_dataset('smooth_diffuse_x', data = smooth_diffuse)
        return diffuse, smooth_diffuse























 
