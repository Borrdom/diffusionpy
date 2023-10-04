:py:mod:`diffusionpy`
=====================

.. py:module:: diffusionpy


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   DasDennisSpacing/index.rst
   FEM_collocation/index.rst
   PyCSAFT_nue/index.rst
   Stefan_Maxwell_segmental/index.rst
   crank_and_other/index.rst
   crystallization/index.rst
   liquidseperation/index.rst
   plots/index.rst
   read_componentdatabase/index.rst
   relaxation/index.rst
   surface_activity/index.rst
   xloil_functions/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   diffusionpy.origin_like



Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.Diffusion_MS
   diffusionpy.D_Matrix
   diffusionpy.Diffusion_MS_iter
   diffusionpy.time_dep_surface
   diffusionpy.Diffusion_MS_xloil
   diffusionpy.reduce_points
   diffusionpy.crank_xl
   diffusionpy.interp1d
   diffusionpy.BHX_xloil
   diffusionpy.DasDennis_xloil
   diffusionpy.Diffusion_MS_iter_xloil
   diffusionpy.lngi
   diffusionpy.eta_iter
   diffusionpy.ares
   diffusionpy.lnphi_TP
   diffusionpy.vpure
   diffusionpy.dlnai_dlnxi
   diffusionpy.dlnai_dlnxi_loop
   diffusionpy.circular
   diffusionpy.time_dep_surface_cryst



.. py:function:: Diffusion_MS(tint, L, Dvec, wi0, wi8, Mi, mobile, full_output=False, dlnai_dlnwi=None, **kwargs)

   Method that computes the multi-component diffusion kinetics

   :param t: time
   :type t: array_like
   :param L: dry thickness /m
   :type L: float
   :param Dvec: Vector of diffusion coefficients. See diffusionpy.D_Matrix                       /m^2/s
   :type Dvec: array_like
   :param wi0: Mass fractions at t=0               /-
   :type wi0: array_like
   :param wi8: Mass fraction at t=infinity         /-
   :type wi8: array_like
   :param Mi: Molar mass of components nc         /g/mol
   :type Mi: array_like
   :param mobile: boolean vector indicating the mobility of a component
   :type mobile: array_like
   :param dlnai_dlnwi: estimate for DlnaiDlnx at t          /-
   :type dlnai_dlnwi: array_like
   :param Keyword Arguments: wiB (array_like): Hello

                             rho0iB (array_like): Hello

   :returns: if ``full_output=False``:

             Matrix ``wt`` of mass fractions at t /-


             if ``full_output=True``:

             Matrix of mass fractions at t       /-

             Matrix of mass fractions at t,z     /-
   :rtype: ndarray

   .. rubric:: Examples

   >>> t=np.linspace(0,300,100)
   >>> Dvec=np.asarray([1E-13])
   >>> wi0=np.asarray([0.3,0.7])
   >>> wi8=np.asarray([0.7,0.3])
   >>> L=1E-6
   >>> Mi=np.asarray([18.,72.])
   >>> mobile=np.asarray([True,False])
   >>> wt=Diffusion_MS(tint,L,Dvec,wi0,wi8,Mi,mobile)

   .. seealso:: diffusionpy.D_Matrix


.. py:function:: D_Matrix(Dvec, nc)

   Creates a symmetric square Matrix ``Dmat`` of dimension ``nc``
   using the elements in the vector ``Dvec``.
   It is assumed that the elements of ``Dvec`` fit precisly
   in the upper and lower triangular entries of ``Dmat``.
   :param Dvec: Must have the length of  ``(nc-1)*nc/2`` to fit in the diagionals of the result matrix
   :type Dvec: array_like
   :param nc: Dimension of ``Dmat``.
   :type nc: int

   :returns: square matrix ``Dmat`` of shape ``(nc,nc)``
   :rtype: ndarray

   :raises Exception: Wrong length of ``Dvec``. Provide array with ``(nc-1)*nc/2`` entries

   .. rubric:: Examples

   >>> Dvec=np.array([1E-13,2E-13,3E-13])
   >>> nc=3
   >>> Dmat=D_Matrix(Dvec,nc)
   >>> Dmat
   array([[0.e+00, 1.e-13, 2.e-13],
          [1.e-13, 0.e+00, 3.e-13],
          [2.e-13, 3.e-13, 0.e+00]])


.. py:function:: Diffusion_MS_iter(t, L, Dvec, wi0, wi8, Mi, mobile, full_output=False, dlnai_dlnwi_fun=None, **kwargs)

   iterates dlnai_dlnwi as a function of the concentration wi
   .. seealso:: diffusionpy.Diffusion_MS


.. py:function:: time_dep_surface(t, wi0, wi8, mobile, taui, lngi_t=None)

   _summary_

   :param t: time vector
   :type t: array_like
   :param wi0: Mass fractions at t=0
   :type wi0: array_like
   :param wi8: Mass fraction at t=infinity
   :type wi8: array_like
   :param mobile: boolean vector indicating the mobility of a component
   :type mobile: array_like
   :param taui: time constant of the surface concentration function
   :type taui: array_like
   :param lngi_fun: activity coefficient function.
   :type lngi_fun: array_like, optional

   :returns: mass fraction at the surface as a vector of time
   :rtype: array_like


.. py:function:: Diffusion_MS_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, full_output: bool = False)


.. py:function:: reduce_points(x, n: int)


.. py:function:: crank_xl(t, L0, Ds, ws0, ws8)


.. py:function:: interp1d(x, xp, fp)


.. py:function:: BHX_xloil(t: xlo.Array(float, dims=1), kf: float, kr: float, ws0: float, ws8: float, mfinfty: float, mrinfty: float)


.. py:function:: DasDennis_xloil(p: int, dim: int)


.. py:function:: Diffusion_MS_iter_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, T: float = 298.15, p: float = 100000.0, pure: xlo.Array(object, dims=2) = np.asarray([[]]), kij: xlo.Array(object, dims=2) = np.asarray([[]]), maxit: int = 10, full_output: bool = False)


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


.. py:function:: eta_iter(p, T, xi, mi, si, ui, eAi, kAi, NAi, kij=np.asarray([[0.0]]), kijA=np.asarray([[0.0]]))

   solve the density mich yiels a given pressure p


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


.. py:function:: lnphi_TP(p, T, xi, mi, si, ui, eAi, kAi, NAi, Mi=None, kij=np.asarray([[0.0]]), kijA=np.asarray([[0.0]]), **kwargs)

   calculate the log of the fugacity coeffficients


.. py:function:: vpure(p, T, mi, si, ui, eAi, kAi, NAi, **kwargs)

   solve the density mich yiels a given pressure p


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


.. py:function:: circular(t, zvec, wtz, Lt=None, instances=6, comp=0, cmap='Blues', vmin=None, vmax=None, label=None, tinterp=None)


.. py:class:: origin_like


   .. py:method:: subplots()


   .. py:method:: set_xlabel(xlabel1, xunit1=None)


   .. py:method:: set_ylabel(ylabel1, yunit1=None)


   .. py:method:: plot(x, y, Formatstring, label=None, order=1, yerr=None, z=None)


   .. py:method:: set_ticks(x0=None, x1=None, y0=None, y1=None)


   .. py:method:: ternary()


   .. py:method:: set_labels(label='mass fractions / -', title='T = 298.15 K \np = 1 bar', xlabel='solvent', ylabel='polymer', zlabel='API')


   .. py:method:: filled_line(x, y, z, Formatstring, legend)


   .. py:method:: conodes(RBx, RBy, RBz, LBx, LBy, LBz, Formatstring, legend)



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


