LFSM.spectral_index.spectral_index_analysis_for_direction_dependence
============================================================================ 
class::

   class LFSM.spectral_index.spectral_index_analysis_for_direction_dependence(object):

method:

.. py:function:: __init__(self)
   
   initial parameter function
   :return: none
.. py:function:: combined_index(downgrade_to)

   return the combined spectral index using error function.

   :param int downgrade_to: The Nside value one choose in healpix mode, must be 2^N.
   :return: the spectral index value
   :rtype: float
   :raises TypeError: if the downgrade_to is not 2^N.

   
   
   


   
   
