:py:mod:`diffusionpy.xloil_functions`
=====================================

.. py:module:: diffusionpy.xloil_functions


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   diffusionpy.xloil_functions.xlo



Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.xloil_functions.Diffusion_MS_xloil
   diffusionpy.xloil_functions.Diffusion_MS_iter_xloil
   diffusionpy.xloil_functions.time_dep_surface_xloil
   diffusionpy.xloil_functions.gradient
   diffusionpy.xloil_functions.reduce_points
   diffusionpy.xloil_functions.interp1d



.. py:class:: xlo


   .. py:attribute:: func

      

   .. py:method:: Array(**kwargs)



.. py:function:: Diffusion_MS_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, full_output: bool = False)


.. py:function:: Diffusion_MS_iter_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, T: float = 298.15, p: float = 100000.0, pure: xlo.Array(object, dims=2) = np.asarray([[]]), kij: xlo.Array(object, dims=2) = np.asarray([[]]), maxit: int = 10, full_output: bool = False)


.. py:function:: time_dep_surface_xloil(t: xlo.Array(float, dims=1), wi0: xlo.Array(float, dims=1), wi8: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), taui: xlo.Array(float, dims=1), lngi_t: xlo.Array(float, dims=2) = None)


.. py:function:: gradient(x: xlo.Array(float, dims=1), y: xlo.Array(float, dims=1))


.. py:function:: reduce_points(x, n: int)


.. py:function:: interp1d(x, xp, fp)


