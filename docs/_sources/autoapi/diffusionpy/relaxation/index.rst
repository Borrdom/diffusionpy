:py:mod:`diffusionpy.relaxation`
================================

.. py:module:: diffusionpy.relaxation


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.relaxation.MEOS_mode
   diffusionpy.relaxation.MDF
   diffusionpy.relaxation.stress
   diffusionpy.relaxation.boundary
   diffusionpy.relaxation.initialboundary



.. py:function:: MEOS_mode(rhovinit, ode, EJ, etaJ, exponent, M2, v2)

   alter the ode function in diffusionpy.Diffusion_MS, to also solve the relaxation

   :param rhovinit: vector of the partial densities of the mobile components
   :type rhovinit: array_like
   :param ode: ode fuinction which is modified by the function
   :type ode: array_like
   :param EJ: Elasticity Moduli
   :type EJ: array_like
   :param etaJ: Damper constants
   :type etaJ: array_like
   :param exponent: Plasticization factors
   :type exponent: array_like
   :param M2: Molar masses of mobile compnents
   :type M2: array_like
   :param v2: specific volumes of mobile compnents
   :type v2: array_like

   :returns: new modified ode function with the same format as the input ode function
   :rtype: array_like


.. py:function:: MDF(sigmaJ, EJ, RV)

   the mechanical driving force for the stress gradient


.. py:function:: stress(etaWL, EJ, sigmaJ, drhodtNF, v2)

   calculate the change in the stresses of the maxwell elements


.. py:function:: boundary(rhov, etaWL, EJ, sigmaJ, RV, THFaktor, v2)

   calculate the mass fraction vector of the surface elment as afunction of time


.. py:function:: initialboundary(RV, EJ, v2, THFaktors, rho, mobiles, wi0, wi8)

   calculate the boudnary of the surface concnetration and surface stresses as a function time


