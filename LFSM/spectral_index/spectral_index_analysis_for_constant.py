# coding: utf-8
#!/usr/bin/env python

import matplotlib
#get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import numpy as np;import healpy as hp;import h5py;
import scipy.constants as Cons
K = Cons.k
C = Cons.c
from scipy.special import erf
import sys
#from LFSM.DIY_colorbar.diy_colorbar import diy_colorbar
import os

class reborn(object):
    def __init__(self):
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        self.file_dir = _path +'/input'
    def unit(self,v):
        return np.square(v*1e6) * 2*K/C**2
    
    def I_E(self,v):
        result = 24.4 *(v*1e-3/0.31)**-2.58
        #print (result)

        return result
    def change_coord(self, m, coord):
        """ Change coordinates of a HEALPIX map

        Parameters
        ----------
        m : map or array of maps
          map(s) to be rotated
        coord : sequence of two character
          First character is the coordinate system of m, second character
          is the coordinate system of the output map. As in HEALPIX, allowed
          coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

        Example
        -------
        The following rotate m from galactic to equatorial coordinates.
        Notice that m can contain both temperature and polarization.
        >>>> change_coord(m, ['G', 'C'])
        """
        # Basic HEALPix parameters
        npix = m.shape[-1]
        nside = hp.npix2nside(npix)
        ang = hp.pix2ang(nside, np.arange(npix))
    
        # Select the coordinate transformation
        rot = hp.Rotator(coord=reversed(coord))

        # Convert the coordinates
        new_ang = rot(*ang)
        new_pix = hp.ang2pix(nside, *new_ang)

        return m[..., new_pix]

    
    def IndexToDeclRa(self,index,downgrade_to):
        theta,phi=hp.pix2ang(downgrade_to,index)
        return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    
    def read_file(self,downgrade_to= 256, resolution = 5):

        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/drao_10MHz.hdf5','r') as f:
        #with h5py.File('drao_10MHz.hdf5','r') as f:
            #print (list(f.keys()))
            hpmap_10 = f['hpmap'][:]
        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/drao_22MHz.hdf5','r') as g:
            #print (list(g.keys()))
            hpmap_22 = g['hpmap'][:]
        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/wlb_45MHz.hdf5','r') as h:
            hpmap_45_old = h['hpmap'][:]
        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/parkes-aus_85MHz.hdf5','r') as i:
            hpmap_85 = i['hpmap'][:]
        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/parkes-aus_150MHz.hdf5','r') as j:
            hpmap_150 = j['hpmap'][:]
        hpmap_408 = hp.read_map(self.file_dir + '/sky_map_10-2300MHz/haslam408_dsds_Remazeilles2014.fits')
        #downgrade to 256 may be, and change the coordinate from galaxy to equatorial
        hpmap_408 = hp.ud_grade(hpmap_408,downgrade_to)
        hpmap_408 = self.change_coord(hpmap_408,['G', 'C'])
        
        with h5py.File(self.file_dir + '/sky_map_10-2300MHz/reich_1420MHz.hdf5','r') as y:
            hpmap_1420 = y['hpmap'][:]
        #####
        #the data from LWA    
        hpmap_35 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-35.fits')
        hpmap_38 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-38.fits')
        hpmap_40 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-40.fits')
        hpmap_45 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-45.fits')
        hpmap_50 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-50.fits')
        hpmap_60 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-60.fits')
        hpmap_70 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-70.fits')
        hpmap_74 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-74.fits')
        hpmap_80 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-80.fits')
       
        return hpmap_10,hpmap_408,hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80

    
    def masked_smoothing(self,U, rad=4.8):     
        V=U.copy()
        V[U!=U]=0
        VV=hp.smoothing(V, fwhm=np.radians(rad))    
        W=0*U.copy()+1
        W[U!=U]=0
        WW=hp.smoothing(W, fwhm=np.radians(rad))    
        return VV/WW
    

    def nan_helper(self,y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]
    
    def smooth(self,map_in, fwhm = 5.0):
        fwhm = np.radians(fwhm)
        return hp.smoothing(map_in,fwhm)
    
    def smoothing_data(self,downgrade_to):
        hpmap_10,hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        
        
        pix_number_45_old = []
        for dec in np.arange(70,90,0.1):
            for ra in np.arange(0,360,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45_old.append(pix_num)
        
        pix_number_45 = []
        for dec in np.arange(-90,-40,0.1):
            for ra in np.arange(0,360,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45.append(pix_num)

        Dict = {'hpmap_10':hpmap_10,'hpmap_408':hpmap_408,'hpmap_22':hpmap_22,'hpmap_45_old':hpmap_45_old,'hpmap_35':hpmap_35,'hpmap_38':hpmap_38,'hpmap_40':hpmap_40,'hpmap_45':hpmap_45,'hpmap_50':hpmap_50,'hpmap_60':hpmap_60,'hpmap_70':hpmap_70,'hpmap_74':hpmap_74,'hpmap_80':hpmap_80}
        for key,X in Dict.items():
            X = self.change_coord(X,['C','G'])
            nans, x= self.nan_helper(X)
            X[nans]= np.interp(x(nans), x(~nans), X[~nans])
            X = self.smooth(X)
            X = hp.ud_grade(X,downgrade_to)
            if key == 'hpmap_45_old':
                X = self.change_coord(X,['G','C'])
                X[pix_number_45_old] = 0.
                X = self.change_coord(X,['C','G'])

            if key == 'hpmap_35' or 'hpmap_38' or 'hpmap_40' or 'hpmap_45' or 'hpmap_50' or 'hpmap_60' or 'hpmap_70' or 'hpmap_74' or 'hpmap_80':
                X = self.change_coord(X,['G','C'])
                X[pix_number_45] = 0.
                X = self.change_coord(X,['C','G'])
            #print key,np.min(X),np.max(X) 
            if key == 'hpmap_35':
                Mask_missing_region_lwa = nans.copy()
                Mask_missing_region_lwa = hp.ud_grade(Mask_missing_region_lwa,downgrade_to)

            Dict[key] = X
            
        #the output map with coordinate of galaxy
        return Dict,Mask_missing_region_lwa
        
        
    def plot(self,downgrade_to):
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
       
        spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa = self.calculate_index(downgrade_to)
        
        hpmap_45_old = self.change_coord(index_45_old_and_408,['G','C']).copy()
        hpmap_45 = self.change_coord(spectral_index_lwa_and_408,['G','C']).copy()
        
        #print (np.where(Mask_missing_region_lwa)[0].shape,'plot:np.where(Mask_missing_region_lwa)[0].shape')
        spectral_index_lwa_and_408[Mask_missing_region_lwa] = hp.UNSEEN
        
        spectral_index_lwa_and_408 = self.change_coord(spectral_index_lwa_and_408,['G','C'])
        plt.figure(2)
        hp.mollview(spectral_index_lwa_and_408,cmap = plt.cm.jet)
        hp.graticule(10,coord='C')
        plt.title('spectral_index_lwa_and_408')
        plt.show()
        
        #hpmap_45_old = self.change_coord(hpmap_45_old,['C','G'])
        
        mask_total = list(np.where(mask)[0]) + list(Mask)
        #index_45_old_and_408[mask] = hp.UNSEEN
        #index_45_old_and_408[Mask] = hp.UNSEEN
        index_45_old_and_408[mask_total] = hp.UNSEEN
        #mask2 = list(np.isnan(index_45_old_and_408)) + list(np.where(index_45_old_and_408>-1.9)[0])
        #index_45_old_and_408[mask2] = np.nan
        #nans, x = self.nan_helper(index_45_old_and_408)
        #index_45_old_and_408[nans]=np.interp(x(nans), x(~nans), index_45_old_and_408[~nans])
        
        #index_45_old_and_408[Mask] = index_45_old_and_408[Mask-2000]
        #unregular_value_index = np.where(index_45_old_and_408> -1.)[0]
        #index_45_old_and_408[unregular_value_index] = hp.UNSEEN
        #print ('unregular_value_index',unregular_value_index)
        index_45_old_and_408 = self.change_coord(index_45_old_and_408,['G','C'])
        plt.figure(3)
        hp.mollview(index_45_old_and_408,cmap = plt.cm.jet)
        hp.graticule(10,coord='C')
        plt.title('spectral_index_45_old_and_408')
        plt.show()
       
       
        new_map = []
        #LWA = -40; Guzman = 67
        Dec_0 = -5;LWA_bottom_limit = -30
        A =2. / ((Dec_0 - LWA_bottom_limit)/Dec_0)
        #for pix_number in range(12*downgrade_to**2):
        for pix_number in range(12*256**2):
            Dec,Ra = self.IndexToDeclRa(pix_number,downgrade_to) 
            pix_value = 0.5*(1 + erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45[pix_number] + 0.5*(1 - erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45_old[pix_number]
            new_map.append(pix_value)
        new_map = np.array(new_map)
        #new_map = self.change_coord(new_map,['C','G'])
        plt.figure(678)
        hp.mollview(new_map,cmap = plt.cm.jet)
        hp.graticule(10,coord='C')
        plt.show()
        return
    def func1(self,beta, x1, x2):
        return (x2-self.I_E(408.)) * (x1/408.)**beta

    def error1(self,beta, x1, x2, y):
        #print 'error1','x1,x2,y.shape',x1.shape,x2.shape,y.shape
        return (self.func1(beta,x1,x2) - (y-self.I_E(x1)))/(y) 
        #return (func1(beta,x1,x2) - (y-I_E(x1)))
        
    def calculate_index(self,downgrade_to):
        Dict,Mask_missing_region_lwa = self.smoothing_data(downgrade_to)
        #print 'in calculate index'
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = Dict['hpmap_408'],Dict['hpmap_22'],Dict['hpmap_45_old'],Dict['hpmap_35'],Dict['hpmap_38'],Dict['hpmap_40'],Dict['hpmap_45'],Dict['hpmap_50'],Dict['hpmap_60'],Dict['hpmap_70'],Dict['hpmap_74'],Dict['hpmap_80']
        
        #hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        hpmap_40 = self.change_coord(hpmap_40,['G','C'])
        hpmap_45 = self.change_coord(hpmap_45,['G','C'])
        hpmap_45_old = self.change_coord(hpmap_45_old,['G','C'])
        hpmap_50 = self.change_coord(hpmap_50,['G','C'])
        hpmap_60 = self.change_coord(hpmap_60,['G','C'])
        hpmap_70 = self.change_coord(hpmap_70,['G','C'])
        pix_number_45 = []
        pix_number_45_old = []
        for dec in np.arange(-40,0,0.1):
            #for ra in np.arange(-100,-30,0.1):
            for ra in np.arange(30,100,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45.append(pix_num)

        for dec in np.arange(10,60,0.1):
            for ra in np.arange(130,170,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45_old.append(pix_num)
 
        hpmap_40[pix_number_45] = 0. 
        hpmap_45[pix_number_45] = 0. 
        hpmap_45_old[pix_number_45_old] = 0. 
        hpmap_50[pix_number_45] = 0. 
        hpmap_60[pix_number_45] = 0.
        hpmap_70[pix_number_45] = 0. 
        hpmap_40 = self.change_coord(hpmap_40,['C','G'])
        hpmap_45 = self.change_coord(hpmap_45,['C','G'])
        hpmap_45_old = self.change_coord(hpmap_45_old,['C','G'])
        hpmap_50 = self.change_coord(hpmap_50,['C','G'])
        hpmap_60 = self.change_coord(hpmap_60,['C','G'])
        hpmap_70 = self.change_coord(hpmap_70,['C','G'])
        hpmap_45_old = np.array(hpmap_45_old,dtype='float64') 
        #freq = np.array([22,35,38,40,45,50,60,70,74,80])
        #freq = np.array([45,35,38,40,45,50,60,70,74,80])
        freq = np.array([45,35,38,40,50,60,70,74,80])
        #spectral_index_lwa_and_408 = []
        X1 = [];X2=[];Y =[]
        for i in range(12*downgrade_to**2):
            mask_condition = np.array([hpmap_45_old[i]-2.725-self.I_E(45), hpmap_35[i]-2.725-self.I_E(35), hpmap_38[i]-2.725-self.I_E(38), hpmap_40[i]-2.725-self.I_E(40), hpmap_50[i]-2.725-self.I_E(50), hpmap_60[i]-2.725-self.I_E(60), hpmap_70[i]-2.725-self.I_E(70), hpmap_74[i]-2.725-self.I_E(74), hpmap_80[i]-2.725-self.I_E(80)])
            mask_ = np.where(mask_condition>0)[0]
            value_freq = np.array([hpmap_45_old[i] - 2.725, hpmap_35[i]-2.725, hpmap_38[i]-2.725, hpmap_40[i]-2.725, hpmap_50[i]-2.725, hpmap_60[i]-2.725, hpmap_70[i]-2.725, hpmap_74[i]-2.725, hpmap_80[i]-2.725])
            value_freq = value_freq[mask_]
            
            
            value_408 = np.ones_like(value_freq) * hpmap_408[i]
            x1 = freq[mask_].copy()
            x2 = value_408.copy()
            y = value_freq.copy()
            X1 = X1 + list(x1)
            X2 = X2 + list(x2)
            Y = Y + list(y)
            beta = [-2.6]
            #if int(value_freq.shape[0])==0:
            #    Para_constant = np.nan
            #    print ('i',i)
            #else:  
            #    Para_constant=leastsq(self.error1,beta,args=(x1,x2,y))[0][0]
            #spectral_index_lwa_and_408.append(Para_constant)
        
        X1 = np.array(X1)
        X2 = np.array(X2)
        Y = np.array(Y)
        Para_constant=leastsq(self.error1,beta,args=(X1,X2,Y))[0][0]
        print ('Para_constant',Para_constant)
        #spectral_index_lwa_and_408 = np.array(spectral_index_lwa_and_408).reshape(-1)
        #index_45_old_and_408,Mask = self.index_between_45_and_408(hpmap_45_old,hpmap_408)
        #return spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa
        return Para_constant

    def IndexToDeclRa(self,index,downgrade_to):
        theta,phi=hp.pix2ang(downgrade_to,index)
        return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    
    def DeclRaToIndex(self,decl,RA,downgrade_to):
        return hp.pixelfunc.ang2pix(downgrade_to,np.radians(-decl+90.),np.radians(360.-RA))
    
        
        




