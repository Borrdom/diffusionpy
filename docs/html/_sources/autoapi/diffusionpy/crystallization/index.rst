:py:mod:`diffusionpy.crystallization`
=====================================

.. py:module:: diffusionpy.crystallization


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.crystallization.CNT
   diffusionpy.crystallization.crystallization_mode



.. py:function:: CNT(t, alpha, A, B, n, T, lnS)

   Calculates the crystallization kinetics based on the classical nucleation  coupled with a simple crystal growth model
   :param t: time /s
   :type t: array_like
   :param alpha: crystal fraction           /-
   :type alpha: array_like
   :param A: crystallization rate constant m^2/s
   :type A: float
   :param B: interfacial tension parameter             /N/m
   :type B: float
   :param n: growth order                         /-
   :type n: float
   :param T: Temperature                 /K
   :type T: float
   :param lnS: log of supersaturation  /-
   :type lnS: array_like

   :returns: recrystallization rate       /-

             growth rate    /-
   :rtype: ndarray


.. py:function:: crystallization_mode(ode, A, B, n, T, saftpar, wv_fun=None)

   alter the ode function in diffusionpy.Diffusion_MS, to also solve the crystallization
   :param wvinit: vector of the mass fractions of the mobile components
   :type wvinit: array_like
   :param ode: ode fuinction which is modified by the function for the rest see CNT
   :type ode: array_like

   :returns: new modified ode function with the same format as the input ode function
   :rtype: array_like


