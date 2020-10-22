LFSM.I_E_term.I_E_equation.I_E
=======================================================
class::

   class LFSM.I_E_term.I_E_equation.I_E(object):

method:

.. py:function:: __init__(self, v, I_E_form, beta_1, v_1)
   
   initial parameter function

   :param v float: The frequency of output sky map
   :param str I_E_form: ('seiffert','seiffert_freq_depend'), the form of extragalactic component except for CMB.
   :param float beta_1: for frequency dependence index situation, in paper taking beta_1 = 0.7
   :param float v_1: for frequency dependence index situation, in paper taking the v_1 = 1.0 MHz
   
.. py:function:: I_E(freq)

   return the isotropic extragalactic emission value in frequency v.

   :param float freq: The frequency of isotropic extragalactic emission.
   :return: the value of isotropic extragalactic emission in freq　
   :rtype: float
   :raises TypeError: if the freq is not a float.

   
   
   


   
   
