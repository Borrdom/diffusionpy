:py:mod:`diffusionpy.surface`
=============================

.. py:module:: diffusionpy.surface


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.surface.time_dep_surface



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


