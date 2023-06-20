:py:mod:`diffusionpy.Stefan_Maxwell_segmental`
==============================================

.. py:module:: diffusionpy.Stefan_Maxwell_segmental


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.Stefan_Maxwell_segmental.drhodt
   diffusionpy.Stefan_Maxwell_segmental.Diffusion_MS
   diffusionpy.Stefan_Maxwell_segmental.D_Matrix
   diffusionpy.Stefan_Maxwell_segmental.Diffusion_MS_iter
   diffusionpy.Stefan_Maxwell_segmental.convert



.. py:function:: drhodt(t, rhov, THFaktor, mobiles, immobiles, Mi, D, allflux, swelling, rho, wi0, dmuext, rhoiB, drhovdtB)

   change in the partial density with time


.. py:function:: Diffusion_MS(t, L, Dvec, wi0, wi8, Mi, mobile, full_output=False, dlnai_dlnwi=None, swelling=False, **kwargs)

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


.. py:function:: Diffusion_MS_iter(t, L, Dvec, wi0, wi8, Mi, mobile, full_output=False, swelling=False, dlnai_dlnwi_fun=None, **kwargs)

   iterates dlnai_dlnwi as a function of the concentration wi
   .. seealso:: diffusionpy.Diffusion_MS


.. py:function:: convert(x, M, axis=0)

   convert fractions. e.g mass fractions into mole fractions


