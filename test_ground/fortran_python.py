import ctypes as ct
import numpy as np

# import the dll
libNE2001 = ct.CDLL('../src.NE2001/libNE2001.so')
rad=57.2957795
#radian per degree

#distance equals 50kpc
dist=50.0

# setup the data
N =np.int(dist/0.01)

print N

nd = ct.pointer( ct.c_int(N) )          # setup the pointer

em1D = np.arange(0, N, dtype=np.float32)  # setup the N-long

l = 10.0
b = 0.0

import pyne2001
EM = pyne2001.get_dm_full(l, b, 50)['EM']
print 'EM',EM

l = l * np.pi/180. #now its radian unit
b = b * np.pi/180.


_ = libNE2001.dmdsm1_(nd, ct.pointer( ct.c_float(l) ), ct.pointer( ct.c_float(b) ), ct.pointer( ct.c_float(dist) ), np.ctypeslib.as_ctypes(em1D))



#i = 0
#while i < len(em1D):
#    print i,(i+1)*dist/N,em1D[i]
#    i = i + 1
print 'em1D',em1D[-1],em1D.shape
# call the function by passing the ctypes pointer using the numpy function:

#_ = fortlib.sqr_1d_arr_(nd, np.ctypeslib.as_ctypes(pyarr))

