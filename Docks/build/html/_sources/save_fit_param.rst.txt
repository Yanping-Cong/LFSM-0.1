LFSM.fitting_params.save_fit_params.free_free
=======================================================
class::

   class LFSM.fitting_params.save_fit_params.free_free(object):

method:

.. py:function:: __init__(self,v,nside,index_type,dist,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,only_fit_Anu)
   
   initial parameter function

   :param v float: The frequency of output sky map
   :param int nside: The Nside value one choose in healpix mode, must be 2^N.
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param int dist: the maximux integrated distance of galaxy, normally setting to 50kpc.
   :param emi_form:  ['exp','sech'] the distribution form of emissivity, normally choosen 'exponantial'.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool R0_R1_equal: fixed True
   :param bool using_raw_diffuse: if True, using the raw input data without smoothing.
   :param bool only_fit_Anu: if True, only fitting the first params of emissivity. normally setting to False.
   
.. py:function:: delta_m()

   return the residual between interpolated map and integrated map of emissivity in healpix mode with Nside was set, the fitting params of emissivity.

   :return: delt_m, params will be output simultaneously.
   :rtype: np.array, np.array

   
   
   


   
   
