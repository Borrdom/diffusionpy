:py:mod:`diffusionpy.crystallization`
=====================================

.. py:module:: diffusionpy.crystallization


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.crystallization.Bound
   diffusionpy.crystallization.Kris



Attributes
~~~~~~~~~~

.. autoapisummary::

   diffusionpy.crystallization.crystallize


.. py:function:: Bound(t, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, lngi_fun=None, wv_fun=None)

   calculate the time dependent surface concentration during crystallization

   :param t: vector of time
   :type t: array_like
   :param mobiles: index array indicating the mobile components
   :type mobiles: array_like
   :param immobiles: index array indicating the immobile component
   :type immobiles: array_like
   :param crystallize: index array indicating the crystallizing components
   :type crystallize: array_like
   :param wi0: initial mass fractions
   :type wi0: array_like
   :param wi8: mass fractions at time equals infinity
   :type wi8: array_like
   :param rho0i: pure component densities
   :type rho0i: array_like
   :param Mi: molar mass of components
   :type Mi: array_like
   :param DAPI: crystallizing components diffusion coefficient in the vector
   :type DAPI: array_like
   :param sigma: interfacial tension of crystal component and mixture
   :type sigma: array_like
   :param kt: crystal growth rate constant
   :type kt: array_like
   :param g: crsystal growth exponent
   :type g: array_like
   :param deltaHSL: melting enthalpy of crystallizing components
   :type deltaHSL: array_like
   :param TSL: melting temperature of crystallizing components
   :type TSL: array_like
   :param cpSL: differfence in liquid/solid heat capacity of crystallizing components
   :type cpSL: array_like
   :param lngi_fun: function of logarithmic activity coefficients
   :type lngi_fun: array_like, optional
   :param wv_fun: function or vector how the concentration of the volatile components changes with the
   :type wv_fun: array_like, optional
   :param concentration of the cryystallizing components:

   :returns: vector of mass fractions at the surface as a function of time
   :rtype: array_like


.. py:function:: Kris(t, alpha, r, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, lngi_fun=None, wv_fun=None)

   function calculating the crystallization


.. py:data:: crystallize

   

