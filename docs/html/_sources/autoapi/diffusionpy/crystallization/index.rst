:py:mod:`diffusionpy.crystallization`
=====================================

.. py:module:: diffusionpy.crystallization


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.crystallization.crystallization_mode
   diffusionpy.crystallization.CNT
   diffusionpy.crystallization.time_dep_surface_cryst
   diffusionpy.crystallization.CNT_surf



.. py:function:: crystallization_mode(wvinit, ode, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, deltaHSL, TSL, cpSL, DAPI, sigma, kt, g, lngi_tz)

   alter the ode function in diffusionpy.Diffusion_MS, to also solve the crystallization

   :param wvinit: vector of the mass fractions of the mobile components
   :type wvinit: array_like
   :param ode: ode fuinction which is modified by the function
   :type ode: array_like
   :param deltaHSL: Melting Enthalpy
   :type deltaHSL: array_like
   :param TSL: Melting temperature
   :type TSL: array_like
   :param DAPI:
   :type DAPI: array_like
   :param sigma: interfacial tension of crystal to its surrounding
   :type sigma: array_like
   :param kt: rate of crystal growth kinetics
   :type kt: array_like
   :param g: order of crystal growth kinetics
   :type g: array_like

   :returns: new modified ode function with the same format as the input ode function
   :rtype: array_like


.. py:function:: CNT(t, alpha, r, mobiles, immobiles, crystallizes, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, lngi, wv)

   Calculates the crystallization kinetics based on the classical nucleation theory for nucleation and a simple crystal growth model


.. py:function:: time_dep_surface_cryst(t, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, lngi_fun=None, wv_fun=None)

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


.. py:function:: CNT_surf(t, alpha, r, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, lngi_fun=None, wv_fun=None)

   Calculates the crystallization kinetics based on the classical nucleation theory for nucleation and a simple crystal growth model


