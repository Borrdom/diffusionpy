:py:mod:`diffusionpy.crank_and_other`
=====================================

.. py:module:: diffusionpy.crank_and_other


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.crank_and_other.crank
   diffusionpy.crank_and_other.Relaxation
   diffusionpy.crank_and_other.BHModel
   diffusionpy.crank_and_other.BHX



.. py:function:: crank(t, kf, mfinfty)

   cranks equation for diffusion


.. py:function:: Relaxation(t, kr, nr, mrinfty)

   "used by the Hophenberg modification


.. py:function:: BHModel(t, kf, kr, nr, mfinfty, mrinfty)

   "Berens and Hopfenberg modification of crank equation for relaxation


.. py:function:: BHX(t, kf, kr, nr, X0, Xinfty, mfinfty, mrinfty)

   "used by the Hophenberg modification


