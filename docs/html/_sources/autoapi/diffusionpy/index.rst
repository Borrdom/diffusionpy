:py:mod:`diffusionpy`
=====================

.. py:module:: diffusionpy


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   Extract_DVS/index.rst
   PCSAFT/index.rst
   ares/index.rst
   crystallization/index.rst
   diffusion/index.rst
   distillation/index.rst
   relaxation/index.rst
   surface/index.rst
   utilities/index.rst
   xloil_functions/index.rst


Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.Diffusion_MS
   diffusionpy.Gammaij
   diffusionpy.time_dep_surface
   diffusionpy.lngi
   diffusionpy.vpure
   diffusionpy.NETVLE
   diffusionpy.supersaturation
   diffusionpy.DasDennis
   diffusionpy.Diffusion_MS_iter_xloil
   diffusionpy.Diffusion_MS_xloil



.. py:function:: Diffusion_MS(tint, L, Dvec, wi0, wi8, mobile, T=298.15, p=100000.0, saftpar=None, **kwargs)

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
   :param mobile: boolean vector indicating the mobility of a component
   :type mobile: array_like
   :param dlnai_dlnwi: estimate for DlnaiDlnx at t          /-
   :type dlnai_dlnwi: array_like
   :param Keyword Arguments: wiB (array_like): weight fraction at z=L

   :returns: Matrix of mass fractions at t       /-

             Matrix of mass fractions at t,z     /-

   .. rubric:: Examples

   >>> t=np.linspace(0,300,100)
   >>> Dvec=np.asarray([1E-13])
   >>> wi0=np.asarray([0.3,0.7])
   >>> wi8=np.asarray([0.7,0.3])
   >>> L=1E-6
   >>> Mi=np.asarray([18.,72.])
   >>> mobile=np.asarray([True,False])
   >>> wt=Diffusion_MS(tint,L,Dvec,wi0,wi8,mobile)

   .. seealso:: diffusionpy.D_Matrix


.. py:function:: Gammaij(T, wi, par)


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


.. py:function:: lngi(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij=np.zeros(10), kijA=np.zeros(10), **kwargs)


.. py:function:: vpure(p, T, mi, si, ui, eAi, kAi, NAi, **kwargs)

   solve the density mich yields a given pressure p


.. py:function:: NETVLE(T, wi, v0p, ksw, mi, si, ui, eAi, kAi, NAi, vpure, Mi, kij=np.zeros(10), kijA=np.zeros(10), n=2)


.. py:function:: supersaturation(T, xi, mi, si, ui, eAi, kAi, NAi, vpure, Mi, deltaHSL, TSL, cpSL, kij=np.zeros(10), kijA=np.zeros(10))


.. py:function:: DasDennis(p, dim)

   create a equidistance n dimensional spacing which satisfies the mass balance constraint

   .. rubric:: Examples

   >>> p=30
   >>> dim=3
   >>> spacevec=DasDennis(p, dim)
   >>> pd.DataFrame(spacevec.T).to_excel("test.xlsx")


.. py:function:: Diffusion_MS_iter_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, T: float = 298.15, p: float = 100000.0, pure: xlo.Array(object, dims=2) = np.asarray([[]]), kij: xlo.Array(object, dims=2) = np.asarray([[]]), maxit: int = 10, full_output: bool = False)


.. py:function:: Diffusion_MS_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, full_output: bool = False)


