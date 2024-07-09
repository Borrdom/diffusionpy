:py:mod:`diffusionpy.PCSAFT`
============================

.. py:module:: diffusionpy.PCSAFT


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.PCSAFT.eta_iter
   diffusionpy.PCSAFT.vpure
   diffusionpy.PCSAFT.SAFTSAC
   diffusionpy.PCSAFT.lngi
   diffusionpy.PCSAFT.dlnai_dlnxi
   diffusionpy.PCSAFT.dlnai_dlnxi_loop
   diffusionpy.PCSAFT.NETVLE
   diffusionpy.PCSAFT.supersaturation



.. py:function:: eta_iter(p, T, xi, **kwargs)

   solve the density mich yiels a given pressure p


.. py:function:: vpure(p, T, mi, si, ui, eAi, kAi, NAi, **kwargs)

   solve the density mich yields a given pressure p


.. py:function:: SAFTSAC(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij=np.zeros(10), kijA=np.zeros(10))

   Calculate the log of the activity coefficients via the SAFT-SAC approximation
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
   :param vpure: pure component molar volumes
   :type vpure: array_like
   :param Mi: Molar mass. Calculates properties on a mass basis when given. Defaults to None.
   :type Mi: array_like, optional
   :param kij: Matrix of binary interaction parameters for dispersion . Defaults to np.asarray([[0.]]).
   :type kij: array_like, optional
   :param kijA: Matrix of binary interaction parameters for association Defaults to np.asarray([[0.]]).
   :type kijA: array_like, optional

   :returns: vector of activity coefficients
   :rtype: array_like


.. py:function:: lngi(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij=np.zeros(10), kijA=np.zeros(10), **kwargs)


.. py:function:: dlnai_dlnxi(T, xi, **kwargs)

   Generate the derivatives of the mole fraction with concentration

   :param T: temperature
   :type T: float
   :param xi: mole/mass fraction. Becomes the mass fraction when the molar mass Mi is not None
   :type xi: array_like
   :param par: dictionary containg pc-saft parameters
   :type par: dic

   :returns: martrix of derivatives of the mole fraction with concentration
   :rtype: array_like


.. py:function:: dlnai_dlnxi_loop(T, xi, **kwargs)


.. py:function:: NETVLE(T, wi, v0p, ksw, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij=np.zeros(10), kijA=np.zeros(10), n=2)


.. py:function:: supersaturation(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, deltaHSL, TSL, cpSL, kij=np.zeros(10), kijA=np.zeros(10))


