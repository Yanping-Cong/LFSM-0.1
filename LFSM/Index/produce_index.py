import h5py
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from pylab import cm
from scipy.optimize import curve_fit
import sys
#sys.path.append('../')
#from Constant_Index.improved_leastsq_method import value_of_index
import os
#sys.path.append("..")
#from Index.produce_index import produce_index
from LFSM.I_E_term.I_E_equation import I_E

class produce_index(object):
    
    def __init__(self, Nside,freq,index_type, I_E_form, beta_1 = 0.7, v_1 = 1.):
        self.nside = Nside
        self.freq = freq
        #self.catalog = catalog
        #self.smin = smin
        self.index_type = index_type
        #self.ssmpath = ssmpath
        #self.input_freq = input_freq
        self.I_E_form = I_E_form
        self.beta_1 = beta_1
        self.v_1 = v_1
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        self.file_dir = _path +'/input'
        #self.index_type = index_type
    #def SSM_(self, freq):
    #    diffuse, source, catalog = SSM(freq = freq, nside = self.nside, catalog = self.catalog, Smin = self.smin, name = "useless.hdf5",inSSMpath = self.ssmpath) 
    #    return diffuse, source, catalog
    #
    #def read_file(self):
    #    dict = {}
    #    for freq in self.input_freq:
    #        diffuse, source, catalog = self.SSM_(freq)
    #        dict.update({np.float(freq):diffuse})
    #    print 'dict',dict
    #    return dict


    #    #f = h5py.File('../init_input/10_22_45_85_150MHz.hdf5','r')
    #    #g = h5py.File('../init_input/diffuse_408MHz.hdf5','r')
    #    #diffuse_10 = f['diffuse_10'][:]
    #    #diffuse_22 = f['diffuse_22'][:]
    #    #diffuse_45 = f['diffuse_45'][:]
    #    #diffuse_85 = f['diffuse_85'][:]
    #    #diffuse_150 = f['diffuse_150'][:]
    #    #diffuse_408 = g['diffuse_408'][:]
    #    #
    #    ##print ('raw_difffuse',diffuse_408)
    #    ##plt.figure(222)
    #    ##hp.mollview(np.log10(diffuse_408),cmap = plt.cm.jet)
    #    ##plt.title('raw_diffuse408' + ' mean=' + str(np.mean(diffuse_408)))
    #    ##plt.show()
    #    #
    #    #f.close()
    #    #g.close()
    #    #return {10:diffuse_10,45:diffuse_45,85:diffuse_85,150:diffuse_150,408:diffuse_408}

    #def diff_index(self):
    #    #for item in read_file.item():
    #    #index : array like [array[].shape = 12*nside**2,...].shape = len(self.input_freq) -1, the 408MHz as the based frequency.
    #    index = []
    #    data = self.read_file()
    #    for (k1,v1) in  data.items():
    #        for (k_base,v_base) in data.items():
    #            if (k_base == float(408)) and (k1 != k_base):
    #                print ('k1,k_base',k1,k_base)
    #                print (np.where(v1==0)[0],np.where(v_base==0)[0],np.where(k1 == 0),np.where(k_base==0))
    #                alpha = np.log10(v1/v_base) / np.log10(k1/k_base)
    #                index.append(alpha)
    #    print 'index', index
    #    print 'len(index)', len(index)
    #    print 'index[0].shape', index[0].shape
    #    return index

    #def _45_408_induce_index(self):
    #    data = self.read_file()
    #    v1 = np.float(45)
    #    v_base = np.float(408)
    #    k1 = data[v1]
    #    k_base = data[v_base]
    #    alpha = np.log10(k1/k_base)/np.log10(v1/v_base) 
    #    return alpha

    #def plot_index(self):
    #    #index = self.mean_index()
    #    index = self._45_408_induce_index()
    #    plt.figure(2)
    #    hp.mollview(index, cmap = plt.cm.jet)
    #    #plt.savefig('fitted_index.eps', format = 'eps')
    #    plt.savefig('45_408_induce_index.eps', format = 'eps')
    #    return

    #def plot_fixed_pixel_index(self):
    #    index = self.diff_index()
    #    import random
    #    n = random.randint(0,index[0].shape[0]-1)
    #    Y = []
    #    for i in index:
    #        Y.append(i[n])
    #    plt.figure(3)
    #    plt.plot(np.arange(len(Y)),Y)
    #    plt.show()
    #    return

    #def mean_index(self):
    #    index = self.diff_index()
    #    index = np.array(index)
    #    fitted_index = []
    #    def func_(x, d):
    #        return d
    #    for i in range(index.shape[1]):
    #        Y = index[:,i]
    #        X = np.arange(Y.size)
    #        popt, pcov = curve_fit(func_, X, Y)
    #        fitted_index.append(popt[0])

    #    print '(X,Y)', (X,Y)
    #    fitted_index = np.array(fitted_index)
    #    print 'fitted_index and fitted_index.shape',fitted_index,fitted_index.shape
    #    #void = np.zeros_like(index[0])
    #    #for i in index:
    #    #    void += i
    #    #void = void / len(index)
    #    with h5py.File('Nside_'+str(self.nside)+'fitted_index.hdf5', 'w') as f:
    #        f.create_dataset('fitted_index',data = fitted_index)
    #    return fitted_index

    ##def constant_index(self):
    #    beta = -2.4697003
    #    index_ = np.array(12*self.nside**2 * [beta])
    #    return index_  
    def constant_index_minus_I_E(self):
        #beta = -2.47900839
        #f = value_of_index(I_E_form = self.I_E_form)
        #beta = f.Para_constant_Para_freq_dependence()[0]
        #beta = float(beta)
        beta = float(-2.4894)
        #print 'beta',beta
        #('Para_constant', array([2.47900839]), 'Para_0', array([2.73479401, 0.25710172])) 
        index_ = np.array(12*self.nside**2 * [beta])
        return index_

    def freq_dependence_index_minus_I_E(self,freq):
        #beta_0, beta_1 = np.array([2.73479401, 0.25710172])
        f = value_of_index(I_E_form = self.I_E_form)
        #beta0 = f.Para_constant_Para_freq_dependence()[0]

        #beta0 = -2.48044886
        #beta_1 = self.beta_1
        #v_1 = self.v_1
        #beta = beta0 + beta_1 * np.exp(-freq/v_1)
        beta_1 = self.beta_1;v_1 = self.v_1
        #beta0 = f.Para_constant_Para_freq_dependence()[0]
        beta0 = float(-2.4894)
        beta = beta0 + beta_1 * np.exp(-freq/v_1)
        #print 'beta0,beta_1,v_1,freq',beta0,beta_1,v_1,freq
        #print 'beta dependence beta:',beta
        #beta_0, beta_1 = f.Para_constant_Para_freq_dependence()[1]
        #beta_0 = float(beta_0)
        #beta_1 = float(beta_1)
        #print 'beta_0',beta_0
        #print 'beta_1',beta_1
        #beta = -(beta_0 + beta_1 * np.log10(self.freq/408.))
        index_ = beta
        return index_

    def pixel_dependence_index_minus_I_E(self):
        with h5py.File('/public/home/wufq/congyanping/Software/LFSAM_new/Direction_dependence_Index/new_map.hdf5','r') as f:
           index = f['new_map'][:]
        return hp.ud_grade(index,self.nside)
    
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
 
    def diffuse_x(self, freq):
        #if self.index_type == 'constant_index':
        #    index = self.constant_index()

        if self.index_type == 'constant_index_minus_I_E':
            index = self.constant_index_minus_I_E()
        if self.index_type == 'freq_dependence_index_minus_I_E':
            index = self.freq_dependence_index_minus_I_E(freq)

        if self.index_type == 'pixel_dependence_index_minus_I_E':
            index = self.pixel_dependence_index_minus_I_E()
            #plt.figure(1)
            #hp.mollview(index,cmap = plt.cm.tab20b)
            #plt.savefig(str(freq)+'MHz_'+'pixel_dependence_index.eps',format = 'eps')
            #freq is in Mhz
            #try:
            #    with h5py.File('Nside_'+str(freq)+'fitted_index.hdf5', 'r') as g:
            #        index = f['fitted_index'][:]
            #except:
            #    index = self.mean_index()

        
        data_freq = 408.
        #the based data of 408MHz coming from SSM model output
        #data_diffuse, source, catalog = self.SSM_(data_freq)

        ##the based data of 408MHz coming from HS14 data set
        data_diffuse = hp.read_map(self.file_dir + '/sky_map_10-2300MHz/haslam408_dsds_Remazeilles2014.fits')
        data_diffuse = hp.ud_grade(data_diffuse, self.nside) 
        #print 'min(data_diffuse_408',min(data_diffuse),'I_E(408)',I_E(data_freq, self.I_E_form).I_E()
        Mask = np.where(data_diffuse -I_E(data_freq, self.I_E_form).I_E() <0)[0]
        data_diffuse[Mask] = np.nan
        
        nans, x= self.nan_helper(data_diffuse)
        data_diffuse[nans]= np.interp(x(nans), x(~nans), data_diffuse[~nans])
        #plt.figure(1)
        #hp.mollview(np.log10(data_diffuse - I_E(data_freq, self.I_E_form).I_E()))
        #plt.savefig('linshi_408_G0rC.png',ppi=600)
        #data_diffuse[Mask] = data_diffuse[Mask-50]
        #print np.where(data_diffuse -I_E(data_freq, self.I_E_form).I_E() <0)[0],'np.where(data_diffuse -I_E(data_freq, self.I_E_form).I_E() <0)[0]'

        data_diffuse = data_diffuse - I_E(data_freq, self.I_E_form).I_E()
        #for freq  in [1,2,3,4,5,10,16,32,100,408]: 
        diffuse_x = np.multiply(data_diffuse, (freq/data_freq)**index)
        value = I_E(freq, self.I_E_form).I_E() 
        with h5py.File('Index_'+'diffuse_'+str(freq)+'MHz.hdf5', 'w') as f:
            f.create_dataset('freq', data = freq)
            f.create_dataset('nside',data = self.nside)
            f.create_dataset('diffuse_x', data = diffuse_x)
            f.create_dataset('diffuse_408',data = data_diffuse)
            f.create_dataset('index',data = index)
            f.create_dataset('I_E',data = np.array([value,np.nan]))
        print ('np.isnan(diffuse_x)',np.where(np.isnan(diffuse_x))) 
        return diffuse_x

