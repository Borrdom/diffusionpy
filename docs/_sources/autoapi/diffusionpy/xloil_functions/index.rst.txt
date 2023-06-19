:py:mod:`diffusionpy.xloil_functions`
=====================================

.. py:module:: diffusionpy.xloil_functions


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.xloil_functions.Diffusion_MS_xloil
   diffusionpy.xloil_functions.Diffusion_MS_iter_xloil
   diffusionpy.xloil_functions.gradient
   diffusionpy.xloil_functions.reduce_points
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


.. py:function:: Diffusion_MS_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), full_output: bool = False, dlnai_dlnwi: xlo.Array(float, dims=2) = None, swelling: bool = False, taui: xlo.Array(float, dims=1) = None, rho0i: xlo.Array(float, dims=1) = None)


.. py:function:: Diffusion_MS_iter_xloil(t: xlo.Array(float, dims=1), L: float, Dvec: xlo.Array(float, dims=1), w0: xlo.Array(float, dims=1), w8: xlo.Array(float, dims=1), Mi: xlo.Array(float, dims=1), mobile: xlo.Array(bool, dims=1), full_output: bool = False, swelling: bool = False, taui: xlo.Array(float, dims=1) = None, rho0i: xlo.Array(float, dims=1) = None, pure: xlo.Array(object, dims=2) = np.asarray([[]]), kij: xlo.Array(object, dims=2) = np.asarray([[]]))


.. py:function:: gradient(x: xlo.Array(float, dims=1), y: xlo.Array(float, dims=1))


.. py:function:: reduce_points(x, n: int)


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


