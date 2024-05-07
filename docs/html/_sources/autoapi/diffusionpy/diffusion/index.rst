:py:mod:`diffusionpy.diffusion`
===============================

.. py:module:: diffusionpy.diffusion


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.diffusion.Diffusion_MS
   diffusionpy.diffusion.D_Matrix
   diffusionpy.diffusion.massbalancecorrection
   diffusionpy.diffusion.Gammaij
   diffusionpy.diffusion.DIdeal2DReal
   diffusionpy.diffusion.wegstein



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


.. py:function:: massbalancecorrection(dlnai_dlnwi, wi0, mobile)


.. py:function:: Gammaij(T, wi, par)


.. py:function:: DIdeal2DReal(Dvec, wave, wi0, dlnai_dlnwi, mobile, realtoideal=False)


.. py:function:: wegstein(fun, x)

   Solving via wegsteins method


