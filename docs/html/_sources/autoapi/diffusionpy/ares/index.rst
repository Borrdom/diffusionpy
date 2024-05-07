:py:mod:`diffusionpy.ares`
==========================

.. py:module:: diffusionpy.ares


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.ares.ares



.. py:function:: ares(T, eta, xi, mi, si, ui, eAi, kAi, NAi, kij, kijA)

   calculate the reduced residual helmholtz energy, the chemical potential and the real gas factor

   :param T: temperature
   :type T: float
   :param xi: mole/mass fraction. Becomes the mass fraction when the molar mass Mi is not None
   :type xi: array_like
   :param mi: segment number
   :type mi: array_like
   :param si: segment diameter
   :type si: array_like
   :param ui: dispersion energy
   :type ui: array_like
   :param eAi: association energy
   :type eAi: array_like
   :param kAi: association volume
   :type kAi: array_like
   :param NAi: association sites (only symmetric)
   :type NAi: array_like
   :param kij: Matrix of binary interaction parameters for dispersion.
   :type kij: array_like
   :param kijA: Matrix of binary interaction parameters for association.
   :type kijA: array_like

   :returns: reduced residual helmholtz energy
             mures (array_like): reduced residual chemical potential
             Zres (aaray_like): real gas factor
   :rtype: ares (array_like)


