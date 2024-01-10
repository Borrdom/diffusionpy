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
   diffusionpy.xloil_functions.origin_plot
   diffusionpy.xloil_functions.crank_xl
   diffusionpy.xloil_functions.interp1d
   diffusionpy.xloil_functions.BHX_xloil
   diffusionpy.xloil_functions.DasDennis_xloil
   diffusionpy.xloil_functions.add_custom
   diffusionpy.xloil_functions.get_path
   diffusionpy.xloil_functions.get_par_xloil
   diffusionpy.xloil_functions.create_header
   diffusionpy.xloil_functions.PC_SAFT_NpT2



Attributes
~~~~~~~~~~

.. autoapisummary::

   diffusionpy.xloil_functions._excelgui


.. py:class:: xlo


   .. py:attribute:: func

      

   .. py:method:: Array(**kwargs)



.. py:function:: Diffusion_MS_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, full_output: bool = False)


.. py:function:: Diffusion_MS_iter_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), swelling: bool = False, witB: xlo.Array(float, dims=2) = None, T: float = 298.15, p: float = 100000.0, pure: xlo.Array(object, dims=2) = np.asarray([[]]), kij: xlo.Array(object, dims=2) = np.asarray([[]]), maxit: int = 10, full_output: bool = False)


.. py:function:: time_dep_surface_xloil(t: xlo.Array(float, dims=1), wi0: xlo.Array(float, dims=1), wi8: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), taui: xlo.Array(float, dims=1), lngi_t: xlo.Array(float, dims=2) = None)


.. py:function:: gradient(x: xlo.Array(float, dims=1), y: xlo.Array(float, dims=1))


.. py:function:: reduce_points(x, n: int)


.. py:function:: origin_plot(x: xlo.Array(float, dims=1), y: xlo.Array(float, dims=1), xexp: xlo.Array(float, dims=1), yexp: xlo.Array(float, dims=1), xmin, xmax, ymin, ymax, xlabel, ylabel, xunit, yunit, spacing: int, fmt)


.. py:function:: crank_xl(t, L0, Ds, ws0, ws8)


.. py:function:: interp1d(x, xp, fp)


.. py:function:: BHX_xloil(t: xlo.Array(float, dims=1), kf: float, kr: float, ws0: float, ws8: float, mfinfty: float, mrinfty: float)


.. py:function:: DasDennis_xloil(p: int, dim: int)


.. py:function:: add_custom(a, b)


.. py:function:: get_path(ctrl)


.. py:function:: get_par_xloil(subst_input, path)


.. py:function:: create_header(ctrl)


.. py:data:: _excelgui

   

.. py:function:: PC_SAFT_NpT2(pure, kij, header, inputs)


