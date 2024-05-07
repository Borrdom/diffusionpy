:py:mod:`diffusionpy.relaxation`
================================

.. py:module:: diffusionpy.relaxation


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.relaxation.relaxation_mode
   diffusionpy.relaxation.MDF
   diffusionpy.relaxation.stress



.. py:function:: relaxation_mode(ode, EJ, etaJ, exponent, T=298.15, Mv=18.015, vv=0.001)

   alter the ode function in diffusionpy.Diffusion_MS, to also solve the relaxation of a component
   :param ode: ode function which is modified by this function
   :type ode: array_like
   :param EJ: Elasticity Moduli
   :type EJ: array_like
   :param etaJ: Damper constants
   :type etaJ: array_like
   :param exponent: Plasticization factors
   :type exponent: array_like
   :param T: temperature /K
   :type T: float, optional
   :param Mv: solvent molar mass g/mol
   :type Mv: float, optional
   :param vv: solvent specific volume kg/m3
   :type vv: float, optional

   :returns: new modified ode function with the same format as the input ode function
   :rtype: array_like


.. py:function:: MDF(sigmaJ, EJ, RV)

   the mechanical driving force for the stress gradient


.. py:function:: stress(etaWL, EJ, sigmaJ, drhodtNF, vv)

   calculate the change in the stresses of the maxwell elements


