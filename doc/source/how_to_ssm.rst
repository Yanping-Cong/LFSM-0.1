LFSM.fitting_params.how_to_smooth_SSM.smooth
=======================================================
class::

   class LFSM.fitting_params.how_to_smooth_SSM.smooth(object):

method:

.. py:function:: __init__(self, Nside, v, index_type, I_E_form, using_raw_diffuse)
   
   initial parameter function

   :param int Nside: The Nside value one choose in healpix mode, must be 2^N.
   :param v float: The frequency of output sky map
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool using_raw_diffuse: if True, using the raw input data without smoothing.
   
.. py:function:: add_5()

   return the (l, b, pixel_value) array for fitting the params of emissivity.

   :return: (longtitude, latitude, pixel_value) in healpix mode with Nside was set
   :rtype: np.array

   
   
   


   
   
