
# coding: utf-8

# In[2]:

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

class reborn_direction(object):
    def __init__(self):
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        self.file_dir = _path +'/input'
    def unit(self,v):
        return np.square(v*1e6) * 2*K/C**2
    
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
        
        #hpmap_35 = self.change_coord(hpmap_35,['C','G'])
        #nans, x= self.nan_helper(hpmap_35)
        #hpmap_35[nans]= np.interp(x(nans), x(~nans), hpmap_35[~nans])
        #index = np.where(np.isnan(hpmap_35))[0]
        #hpmap_35[index] = hp.UNSEEN
        #hpmap_35 = hp.ma(hpmap_35)
        #ang = hp.pix2ang(256,index,lonlat = True)
        #print ('hpmap_35',hpmap_35)
        #val = hp.get_interp_val(hpmap_35,ang[0],ang[1],lonlat=True)
        #print ('val',val)     
        #plt.figure(2)
        #plt.plot(np.arange(val.shape[0]),val,'o')
        #plt.show()
        #hpmap_35[index] = val
        #hpmap_38 = self.change_coord(hpmap_38,['C','G'])
        #Mask_missing_region_lwa = np.isnan(hp.ud_grade(hpmap_38,downgrade_to))
        
        
        
        Dict = {'hpmap_10':hpmap_10,'hpmap_408':hpmap_408,'hpmap_22':hpmap_22,'hpmap_45_old':hpmap_45_old,'hpmap_35':hpmap_35,'hpmap_38':hpmap_38,'hpmap_40':hpmap_40,'hpmap_45':hpmap_45,'hpmap_50':hpmap_50,'hpmap_60':hpmap_60,'hpmap_70':hpmap_70,'hpmap_74':hpmap_74,'hpmap_80':hpmap_80}
        for key,X in Dict.items():
            X = self.change_coord(X,['C','G'])
            nans, x= self.nan_helper(X)
            X[nans]= np.interp(x(nans), x(~nans), X[~nans])
            X = self.smooth(X)
            if key == 'hpmap_10':
                pix_number = np.where(nans==True)[0]
                pix_vec = hp.pix2vec(256,pix_number)
                pix_number = hp.vec2pix(2**6,pix_vec[0],pix_vec[1],pix_vec[2])
                nans = np.zeros(12*64**2,dtype = np.bool)
                nans[pix_number] = True

            X = hp.ud_grade(X,downgrade_to)
            if key == 'hpmap_35':
                print (np.where(nans)[0].shape,'smoothing_data:np.where(nans)[0].shape')
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
        
        print (np.where(Mask_missing_region_lwa)[0].shape,'plot:np.where(Mask_missing_region_lwa)[0].shape')
        spectral_index_lwa_and_408[Mask_missing_region_lwa] = hp.UNSEEN
        
        spectral_index_lwa_and_408 = self.change_coord(spectral_index_lwa_and_408,['G','C'])
        _min = np.round(np.nanmin(spectral_index_lwa_and_408),1)
        _max = np.round(np.nanmax(spectral_index_lwa_and_408),1)
        print '_min',_min,'_max',_max
        plt.figure(1234)
        #hp.mollview(spectral_index_lwa_and_408,min = _min, max = _max, cbar = False,cmap = plt.cm.jet)
        #hp.mollview(spectral_index_lwa_and_408,min = -10, max = 10,cmap = plt.cm.jet)
        plt.plot(np.arange(spectral_index_lwa_and_408.size),spectral_index_lwa_and_408)
        #diy_colorbar(_min, _max, unit=' ')
        #hp.graticule(10,coord='C')
        #plt.title('spectral_index_lwa_and_408')
        plt.savefig('lwa_408.png',ppi = 600)
        #plt.show()
        raise 'ERROR' 
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
    def I_E(self,v):
        result = 24.4 *(v*1e-3/0.31)**-2.58
        #print (result)
        return result
    def func1(self,beta, x1, x2):
        return (x2-self.I_E(408.)) * (x1/408.)**beta

    def error1(self,beta, x1, x2, y):
        #print 'error1','x1,x2,y.shape',x1.shape,x2.shape,y.shape
        return (self.func1(beta,x1,x2) - (y-self.I_E(x1)))/(y) 
        #return (func1(beta,x1,x2) - (y-I_E(x1)))

    def plot_two_intermediate_index(self,downgrade_to): 
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        mask_45_old = np.isnan(hpmap_45_old)
        mask_45 = np.isnan(hpmap_45)
        hpmap_45_old[mask_45_old] = 0.
        mask_3 = list(np.where(hpmap_45_old - 2.725 - self.I_E(45) <350)[0])
        print 'mask_3.len',len(mask_3)
        Mask_3 = []
        for num in mask_3:
            dec,ra = self.IndexToDeclRa(num,downgrade_to)
            if (15<dec<60):# and 120<ra<180):
                Mask_3.append(num)
        print 'Mask_3.len',len(Mask_3)
            
        spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa = self.calculate_index(downgrade_to)
        cmap_reversed = matplotlib.cm.get_cmap('jet_r') 
        plt.figure(1324)
        spectral_index_lwa_and_408 = self.change_coord(spectral_index_lwa_and_408,['G','C'])
        spectral_index_lwa_and_408[mask_45] = np.nan
        #_min = round(np.nanmin(spectral_index_lwa_and_408),1)
        #_max = round(np.nanmax(spectral_index_lwa_and_408),1)
        _min = -3.1
        _max = -2.1
        hp.mollview(spectral_index_lwa_and_408,min = _min, max = _max, cbar=False,cmap=cmap_reversed)
        plt.title(' ')
        diy_colorbar(_min,_max,unit = ' ')
        plt.savefig('45_lwa_index.eps',format='eps')

        plt.figure(1327)
        index_45_old_and_408 = self.change_coord(index_45_old_and_408,['G','C'])
        index_45_old_and_408[Mask_3] = np.nan
        index_45_old_and_408[mask_45_old] = np.nan
        #_min = round(np.nanmin(index_45_old_and_408),1)
        #_max = round(np.nanmax(index_45_old_and_408),1)
        hp.mollview(index_45_old_and_408, min = _min, max = _max, cbar = False, cmap = cmap_reversed)
        #plt.plot(np.arange(index_45_old_and_408.size),index_45_old_and_408)
        plt.title(' ')
        diy_colorbar(_min,_max,unit = ' ')
        plt.savefig('45_guzman_index.eps',format='eps')
        return

    def calculate_index(self,downgrade_to):
        
        Dict,Mask_missing_region_lwa = self.smoothing_data(downgrade_to)
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = Dict['hpmap_408'],Dict['hpmap_22'],Dict['hpmap_45_old'],Dict['hpmap_35'],Dict['hpmap_38'],Dict['hpmap_40'],Dict['hpmap_45'],Dict['hpmap_50'],Dict['hpmap_60'],Dict['hpmap_70'],Dict['hpmap_74'],Dict['hpmap_80']
        
        #hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        hpmap_40 = self.change_coord(hpmap_40,['G','C'])
        hpmap_45 = self.change_coord(hpmap_45,['G','C'])
        hpmap_50 = self.change_coord(hpmap_50,['G','C'])
        hpmap_60 = self.change_coord(hpmap_60,['G','C'])
        hpmap_70 = self.change_coord(hpmap_70,['G','C'])
        pix_number_45 = []

        for dec in np.arange(-40,0,0.1):
            for ra in np.arange(-100,-30,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45.append(pix_num)
        
        #a = np.zeros_like(hpmap_70)
        #a[pix_number_45] = 1.
        #plt.figure(323)
        #hp.mollview(a)
        #plt.savefig('a.png',ppi=600)
        #raise 'ERROR'

        hpmap_40[pix_number_45] = 0. 
        hpmap_45[pix_number_45] = 0. 
        hpmap_50[pix_number_45] = 0. 
        hpmap_60[pix_number_45] = 0.
        hpmap_70[pix_number_45] = 0. 
        hpmap_40 = self.change_coord(hpmap_40,['C','G'])
        hpmap_45 = self.change_coord(hpmap_45,['C','G'])
        hpmap_50 = self.change_coord(hpmap_50,['C','G'])
        hpmap_60 = self.change_coord(hpmap_60,['C','G'])
        hpmap_70 = self.change_coord(hpmap_70,['C','G'])
        
        #freq = np.array([22,35,38,40,45,50,60,70,74,80])
        freq = np.array([35,38,40,45,50,60,70,74,80])
        spectral_index_lwa_and_408 = []
        for i in range(12*downgrade_to**2):
            #mask_condition = np.array([hpmap_22[i] - 2.725-self.I_E(22), hpmap_35[i]-2.725-self.I_E(35), hpmap_38[i]-2.725-self.I_E(38), hpmap_40[i]-2.725-self.I_E(40), hpmap_45[i]-2.725-self.I_E(45), hpmap_50[i]-2.725-self.I_E(50), hpmap_60[i]-2.725-self.I_E(60), hpmap_70[i]-2.725-self.I_E(70), hpmap_74[i]-2.725-self.I_E(74), hpmap_80[i]-2.725-self.I_E(80)])
            mask_condition = np.array([hpmap_35[i]-2.725-self.I_E(35), hpmap_38[i]-2.725-self.I_E(38), hpmap_40[i]-2.725-self.I_E(40), hpmap_45[i]-2.725-self.I_E(45), hpmap_50[i]-2.725-self.I_E(50), hpmap_60[i]-2.725-self.I_E(60), hpmap_70[i]-2.725-self.I_E(70), hpmap_74[i]-2.725-self.I_E(74), hpmap_80[i]-2.725-self.I_E(80)])
            mask_ = np.where(mask_condition>0)[0]
            #value_freq = np.array([hpmap_22[i] - 2.725, hpmap_35[i]-2.725, hpmap_38[i]-2.725, hpmap_40[i]-2.725, hpmap_45[i]-2.725, hpmap_50[i]-2.725, hpmap_60[i]-2.725, hpmap_70[i]-2.725, hpmap_74[i]-2.725, hpmap_80[i]-2.725])
            value_freq = np.array([hpmap_35[i]-2.725, hpmap_38[i]-2.725, hpmap_40[i]-2.725, hpmap_45[i]-2.725, hpmap_50[i]-2.725, hpmap_60[i]-2.725, hpmap_70[i]-2.725, hpmap_74[i]-2.725, hpmap_80[i]-2.725])
            value_freq = value_freq[mask_]
            
            
            value_408 = np.ones_like(value_freq) * hpmap_408[i]
            x1 = freq[mask_].copy()
            x2 = value_408.copy()
            y = value_freq.copy()
            
            beta = [-2.6]
            if int(value_freq.shape[0])==0:
                Para_constant = np.nan
                print ('i',i)
            else:  
                Para_constant=leastsq(self.error1,beta,args=(x1,x2,y))[0][0]
            spectral_index_lwa_and_408.append(Para_constant)
        spectral_index_lwa_and_408 = np.array(spectral_index_lwa_and_408).reshape(-1)
        index_45_old_and_408,Mask = self.index_between_45_and_408(hpmap_45_old,hpmap_408)
        return spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa
    
    def IndexToDeclRa(self,index,downgrade_to):
        theta,phi=hp.pix2ang(downgrade_to,index)
        return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    
    def DeclRaToIndex(self,decl,RA,downgrade_to):
        return hp.pixelfunc.ang2pix(downgrade_to,np.radians(-decl+90.),np.radians(360.-RA))
    
    def combined_index(self,downgrade_to):
        spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa = self.calculate_index(downgrade_to)
        
        
        hpmap_45_old = self.change_coord(index_45_old_and_408,['G','C']).copy()
        hpmap_45 = self.change_coord(spectral_index_lwa_and_408,['G','C']).copy()
        
        new_map = []
        #LWA = -40; Guzman = 67
        Dec_0 = -15;LWA_bottom_limit = -40
        A =2. / ((Dec_0 - LWA_bottom_limit)/Dec_0)
        #for pix_number in range(12*downgrade_to**2):
        for pix_number in range(12*256**2):
            Dec,Ra = self.IndexToDeclRa(pix_number,downgrade_to) 
            pix_value = 0.5*(1 + erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45[pix_number] + 0.5*(1 - erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45_old[pix_number]
            new_map.append(pix_value)
        new_map = np.array(new_map)
        new_map = self.change_coord(new_map,['C','G'])
        return new_map    
    def paper_plot(self):
        with h5py.File('./new_map.hdf5','r') as f:
            new_map = f['new_map'][:]
        plt.figure(789)
        hp.mollview(new_map,min = np.round(np.min(new_map),2), max = np.round(np.max(new_map),2), cbar = False, cmap = plt.cm.jet)
        plt.title(' ')
        diy_colorbar(unit=' ')
        plt.savefig('new_map.eps',format = 'eps')
        
    def index_between_45_and_408(self,map_45_old,map_408):
        #interpolate the pixel less than I_E(45)
        mask = np.where(map_45_old - 2.725-self.I_E(45)<=0)[0]
        map_45_old[mask] = np.nan
        nans, x = self.nan_helper(map_45_old)
        map_45_old[nans]=np.interp(x(nans), x(~nans), map_45_old[~nans])
        
        mask_45_less_I_E = mask.copy()
        #mask = np.where(map_408-self.I_E(408)<0.1)[0]
        #map_408[mask] = np.nan
        #nans, x = self.nan_helper(map_408)
        #map_408[nans]=np.interp(x(nans), x(~nans), map_408[~nans])
        #print ('mask_408_index',mask)
        
        index = np.log10((map_45_old-2.725-self.I_E(45))/(map_408-self.I_E(408)))/np.log10(45./408.)
        index = index.reshape(-1)
        return index,mask_45_less_I_E
        

    def test2_realtive_difference_between_two_map_of_45Mhz(self,downgrade_to,plot_=True):
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        M_lwa = np.isnan(hpmap_35)
        M_45_old = np.isnan(hpmap_45_old)
        #hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        Dict,Mask_missing_region_lwa = self.smoothing_data(downgrade_to)
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = Dict['hpmap_408'],Dict['hpmap_22'],Dict['hpmap_45_old'],Dict['hpmap_35'],Dict['hpmap_38'],Dict['hpmap_40'],Dict['hpmap_45'],Dict['hpmap_50'],Dict['hpmap_60'],Dict['hpmap_70'],Dict['hpmap_74'],Dict['hpmap_80']
        if plot_ == True:
            hpmap_45_old = self.change_coord(hpmap_45_old,['G','C'])
            hpmap_45_old[M_45_old] = hp.UNSEEN
            _min = round(np.nanmin(np.log10(hpmap_45_old)),1)
            _max = round(np.nanmax(np.log10(hpmap_45_old)),1) - 0.1#stay with LWA upper limit
            print '_min',_min,'_max',_max
            plt.figure(1)
            hp.mollview(np.log10(hpmap_45_old), min= _min, max=_max, cbar = False, cmap=plt.cm.jet)
            plt.title(' ')
            diy_colorbar_C(_min, _max, unit=' ')
            plt.savefig('hpmap_45Mhz_guzman.eps',format = 'eps')
            hpmap_45 = self.change_coord(hpmap_45,['G','C'])
            hpmap_45[M_lwa] = hp.UNSEEN
            _min = round(np.nanmin(np.log10(hpmap_45)),1)
            _max = round(np.nanmax(np.log10(hpmap_45)),1)
            print '_min',_min,'_max',_max
            plt.figure(1)
            hp.mollview(np.log10(hpmap_45), min = _min, max = _max, cbar = False, cmap=plt.cm.jet)
            plt.title(' ')
            diy_colorbar_C(_min, _max, unit=' ')
            plt.savefig('hpmap_45Mhz_lwa.eps',format = 'eps')

        #
        #plt.figure(1)
        #hp.mollview(np.log10(hpmap_408-self.I_E(408)),cmap=plt.cm.jet)
        #plt.title('408')
        #plt.show()
        #
        #index = np.log10((hpmap_80-2.725-self.I_E(80))/(hpmap_408-self.I_E(408)))/np.log10(80./408.)
        #index = self.change_coord(index,['G','C'])
        #index[M_lwa] = hp.UNSEEN
        #plt.figure(2)
        #hp.mollview(index,cmap=plt.cm.jet)
        #plt.title('index 80 and 408')
        #plt.show()
        #
        #index = np.log10((hpmap_70-2.725-self.I_E(70))/(hpmap_408-self.I_E(408)))/np.log10(70./408.)
        #index = self.change_coord(index,['G','C'])
        #index[M_lwa] = hp.UNSEEN
        #plt.figure(2)
        #hp.mollview(index,cmap=plt.cm.jet)
        #plt.title('index 70 and 408')
        #plt.show()
        #
        mask_1 = np.where(M_lwa==True)[0]
        print 'mask_1.shape',mask_1.shape
        mask_2 = np.where(M_45_old==True)[0]
        print 'mask_2.shape',mask_2.shape
        mask = np.array(list(set(list(mask_1) + list(mask_2))))
        print 'mask.shape',mask.shape,mask[0]
        MASK = np.zeros_like(M_lwa,dtype = np.bool)
        MASK[mask] = True
         
        T_input = hpmap_45_old[~MASK]
        T_fit = hpmap_45[~MASK]
        print 'min(T_input)',np.min(T_input) 
        print 'min(T_fit)',np.min(T_fit) 
        N_pix = hpmap_45.size - mask.size
        #T_input = hpmap_45_old
        #T_fit = hpmap_45
        total_relative_difference = np.sqrt(1./N_pix * np.sum(np.square( (T_fit - T_input)/T_input )))*100.
        print 'total_relative_difference_Xuelei',total_relative_difference
        total_relative_difference = 1./(1./N_pix*np.sum(T_input)) * np.sqrt(1./N_pix*np.sum(np.square(T_input-T_fit)))*100.
        print 'total_relative_difference_Binyue',total_relative_difference
        print 'np.mean(hpmap_45)',np.mean(T_fit)
        print 'np.mean(hpmap_45_old)',np.mean(T_input)
        return
    
    def test(self,downgrade_to):
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        hpmap_45_old = self.change_coord(hpmap_45_old,['C','G'])
        M = np.isnan(hpmap_45_old)
        M_lwa = np.isnan(hpmap_35)
        print ('M',np.where(M)[0])
        
        Dict,Mask_missing_region_lwa = self.smoothing_data(downgrade_to)
        hpmap_408, hpmap_22,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = Dict['hpmap_408'],Dict['hpmap_22'],Dict['hpmap_45_old'],Dict['hpmap_35'],Dict['hpmap_38'],Dict['hpmap_40'],Dict['hpmap_45'],Dict['hpmap_50'],Dict['hpmap_60'],Dict['hpmap_70'],Dict['hpmap_74'],Dict['hpmap_80']
        mask = np.where(hpmap_45_old - 2.725-self.I_E(45)<=0)[0]
        
         
        index = np.log10((hpmap_35-2.725-self.I_E(35))/(hpmap_408-self.I_E(408)))/np.log10(35./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('35Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_38-2.725-self.I_E(38))/(hpmap_408-self.I_E(408)))/np.log10(38./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('38Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        
        index = np.log10((hpmap_40-2.725-self.I_E(40))/(hpmap_408-self.I_E(408)))/np.log10(40./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('40Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_45-2.725-self.I_E(45))/(hpmap_408-self.I_E(408)))/np.log10(45./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('45Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_50-2.725-self.I_E(50))/(hpmap_408-self.I_E(408)))/np.log10(50./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('50Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_60-2.725-self.I_E(60))/(hpmap_408-self.I_E(408)))/np.log10(60./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('60Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_70-2.725-self.I_E(70))/(hpmap_408-self.I_E(408)))/np.log10(70./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('70Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_74-2.725-self.I_E(74))/(hpmap_408-self.I_E(408)))/np.log10(74./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('74Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        index = np.log10((hpmap_80-2.725-self.I_E(80))/(hpmap_408-self.I_E(408)))/np.log10(80./408.)     
        index = self.change_coord(index,['G','C'])
        index[M_lwa] = hp.UNSEEN
        plt.figure(1)
        hp.mollview(index,cmap = plt.cm.jet)
        plt.title('80Mhz and 408Mhz calculated spectral index')
        plt.show()
        
        bool_ = np.zeros_like(hpmap_45_old)
        bool_[mask]=1
        plt.figure(2)
        hp.mollview(bool_)
        plt.title('hpmap_45_old - 2.725-self.I_E(45)<=0')
        plt.show()
        
        plt.figure(3)
        hp.mollview(hpmap_45_old - 2.725-self.I_E(45),cmap = plt.cm.tab20c)
        plt.title('T_of_45MHz - CMB - I_E')
        plt.show()
        '''
        index = np.log10((hpmap_45_old-2.725-self.I_E(45))/(hpmap_408-self.I_E(408)))/np.log10(45./408.)
        index = index.reshape(-1)
        index = self.change_coord(index,['G','C'])
        plt.figure(4)
        hp.mollview(index,cmap =plt.cm.tab20b)
        plt.title('raw 45Mhz and 408Mhz to calculate index')
        plt.show()
        plt.figure(44)
        hp.mollview(index,max=-1,cmap =plt.cm.tab20b)
        plt.title('raw 45Mhz and 408Mhz to calculate index')
        plt.show()
        
        '''
        index_45_old_and_408,Mask = self.index_between_45_and_408(hpmap_45_old,hpmap_408)
        index_45_old_and_408[Mask]=hp.UNSEEN
        
        index_45_old_and_408[M] = hp.UNSEEN
        index_45_old_and_408 = self.change_coord(index_45_old_and_408,['G','C'])
        plt.figure(5)
        hp.mollview(index_45_old_and_408,cmap =plt.cm.jet)
        plt.title('Guzman 45Mhz and 408Mhz to calculate index')
        hp.graticule(10,coord='C')
        plt.show()
        
        plt.figure(6)
        hp.mollview(index_45_old_and_408,max = -2.04,cmap =plt.cm.jet)
        hp.graticule(10,coord='C')
        plt.title('Guzman 45Mhz and 408Mhz to calculate index')
        plt.show()
        
        
        
#if __name__ == '__main__':
#    f = reborn()
#    f.combined_index(256)
    



    #f.paper_plot(256)
    #f.test2(256)
    #f.test(256)
    #f.read_file()
    #f.plot(downgrade_to = 256)
    #f.smoothing_data(2**6)
