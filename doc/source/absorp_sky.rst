LFSM.absorp_sky_map.produce_healpix_sky_map.absorption_JRZ
===========================================================
class::

   class LFSM.absorp_sky_map.produce_healpix_sky_map.absorption_JRZ(object):

method:

.. py:function:: __init__(self,v,nside,clumping_factor,index_type,distance,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,test,only_fit_Anu)
   
   initial parameter function

   :param v float: The frequency of output sky map
   :param int nside: The Nside value one choose in healpix mode, must be 2^N.
   :param int clumping_factor: the clumping factor influnce the absorption, the value set to one in our model. 
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param int distance: the maximux integrated distance of galaxy, normally setting to 50kpc.
   :param emi_form:  ['exp','sech'] the distribution form of emissivity, normally choosen 'exponantial'.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool R0_R1_equal: fixed True
   :param bool using_raw_diffuse: if True, using the raw input data without smoothing.
   :test bool test: if True, one can do some test with different parameter using in the code. normally fixed False.
   :param bool only_fit_Anu: if True, only fitting the first params of emissivity. normally setting to False.
   
.. py:function:: mpi()

   return [pixel_number, pixel_value] in healpix mode 

   :return: a numpy array shape as [12*Nside^2, 2] will be return.
   :rtype: np.array

   
   
   


   
   
