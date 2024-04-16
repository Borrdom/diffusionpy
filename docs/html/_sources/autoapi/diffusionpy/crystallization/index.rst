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
   diffusionpy.crystallization.time_dep_surface_cryst
   diffusionpy.crystallization.Diffusion_MS_cryst
   diffusionpy.crystallization.cryst_iter



.. py:function:: CNT(t, alpha, r, mobiles, immobiles, crystallizes, wi0, wi8, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, tnuc, temp, lngi, wv)

   Calculates the crystallization kinetics based on the classical nucleation  coupled with a simple crystal growth model
   :param t: time /s
   :type t: array_like
   :param alpha: crystal fraction /-
   :type alpha: array_like
   :param mobiles: stores index of mobile components
   :type mobiles: array_like
   :param immobiles: stores index of immobile components
   :type immobiles: array_like
   :param crystallizes: stores index of crystallizing components
   :type crystallizes: array_like
   :param wi0: Mass fractions at t=0               /-
   :type wi0: array_like
   :param wi8: Mass fraction at t=infinity         /-
   :type wi8: array_like
   :param rho0i: pure component densitys           /kg/m^3
   :type rho0i: array_like
   :param Mi: Molar mass of components nc         /g/mol
   :type Mi: array_like
   :param DAPI: diffusion coefficient of the API to the crystal m^2/s
   :type DAPI: float
   :param sigma: interfacial tension             /N/m
   :type sigma: float
   :param kt: growth rate constant                /m/s
   :type kt: float
   :param g: growth order                         /-
   :type g: float
   :param deltaHSL: melting enthalpy           /J/mol
   :type deltaHSL: float
   :param TSL: melting temperature           /K
   :type TSL: float
   :param cpSL: heat capacity difference solid and liquid           /J/mol/K
   :type cpSL: float
   :param tnuc: nucleation onset           /s
   :type tnuc: float
   :param temp: temperature                 /K
   :type temp: float
   :param lngi: log of acticity coefficients for supersaturation calculations /-
   :type lngi: array_like
   :param wv: mobile component weight fraction           /-
   :type wv: float

   :returns: recrystallization rate       /-

             growth rate    /-
   :rtype: ndarray


.. py:function:: crystallization_mode(wvinit, ode, mobiles, immobiles, crystallize, wi0, wi8, rho0i, Mi, deltaHSL, TSL, cpSL, tnuc, temp, DAPI, sigma, kt, g, lngi_tz)

   alter the ode function in diffusionpy.Diffusion_MS, to also solve the crystallization

   :param wvinit: vector of the mass fractions of the mobile components
   :type wvinit: array_like
   :param ode: ode fuinction which is modified by the function
   :type ode: array_like
   :param for the rest see CNT:

   :returns: new modified ode function with the same format as the input ode function
   :rtype: array_like


.. py:function:: time_dep_surface_cryst(t, mobile, wi0, wi8, crystallize, rho0i, Mi, DAPI, sigma, kt, g, deltaHSL, TSL, cpSL, tnuc=0.0, temp=298.15, lngi=None, wv_fun=None)

   calculate the time dependent surface concentration during crystallization

   :param t: vector of time
   :type t: array_like
   :param mobile: boolean array indicating the mobile components
   :type mobile: array_like
   :param wi0: initial mass fractions
   :type wi0: array_like
   :param wi8: mass fractions at time equals infinity
   :type wi8: array_like
   :param crystallize: index array indicating the crystallizing components
   :type crystallize: array_like
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


.. py:function:: Diffusion_MS_cryst(t, L, Dvec, wi0, wi8, Mi, mobile, crystpar, lngi_fun, **kwargs)


.. py:function:: cryst_iter(t, mobile, wi0, wi8, crystpar, Mi, lngi_fun, wv_fun)


