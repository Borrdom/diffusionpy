:py:mod:`diffusionpy.PyCSAFT_nue`
=================================

.. py:module:: diffusionpy.PyCSAFT_nue


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.PyCSAFT_nue.ares
   diffusionpy.PyCSAFT_nue.eta_iter
   diffusionpy.PyCSAFT_nue.vpure
   diffusionpy.PyCSAFT_nue.lngi
   diffusionpy.PyCSAFT_nue.lnphi_TP
   diffusionpy.PyCSAFT_nue.dlnai_dlnxi
   diffusionpy.PyCSAFT_nue.dlnai_dlnxi_loop
   diffusionpy.PyCSAFT_nue.initialize



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


.. py:function:: eta_iter(p, T, xi, mi, si, ui, eAi, kAi, NAi, kij=np.asarray([[0.0]]), kijA=np.asarray([[0.0]]))

   solve the density mich yiels a given pressure p


.. py:function:: vpure(p, T, mi, si, ui, eAi, kAi, NAi, **kwargs)

   solve the density mich yiels a given pressure p


.. py:function:: lngi(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij, kijA)

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


.. py:function:: lnphi_TP(p, T, xi, mi, si, ui, eAi, kAi, NAi, Mi=None, kij=np.asarray([[0.0]]), kijA=np.asarray([[0.0]]), **kwargs)

   calculate the log of the fugacity coeffficients


.. py:function:: dlnai_dlnxi(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij, kijA)

   Generate the derivatives of the mole fraction with concentration

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
   :param idx: index which components mass balance is considered. If None mass balance is ignored. Defaults to None.
   :type idx: int, optional

   :returns: martrix of derivatives of the mole fraction with concentration
   :rtype: array_like


.. py:function:: dlnai_dlnxi_loop(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij, kijA)

   Generate the derivatives of the mole fraction with concentration

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
   :param idx: index which components mass balance is considered. If None mass balance is ignored. Defaults to None.
   :type idx: int, optional

   :returns: martrix of derivatives of the mole fraction with concentration
   :rtype: array_like


.. py:function:: initialize()


