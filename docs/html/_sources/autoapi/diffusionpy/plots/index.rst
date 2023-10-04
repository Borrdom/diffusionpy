:py:mod:`diffusionpy.plots`
===========================

.. py:module:: diffusionpy.plots


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   diffusionpy.plots.origin_like



Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.plots.circular
   diffusionpy.plots.basic_colors



.. py:function:: circular(t, zvec, wtz, Lt=None, instances=6, comp=0, cmap='Blues', vmin=None, vmax=None, label=None, tinterp=None)


.. py:function:: basic_colors(Formatstring)


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



