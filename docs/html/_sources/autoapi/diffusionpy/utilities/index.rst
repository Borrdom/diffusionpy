:py:mod:`diffusionpy.utilities`
===============================

.. py:module:: diffusionpy.utilities


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   diffusionpy.utilities.ternary



Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.utilities.DasDennis
   diffusionpy.utilities.get_line_from_tern



.. py:function:: DasDennis(p, dim)

   create a equidistance n dimensional spacing which satisfies the mass balance constraint

   .. rubric:: Examples

   >>> p=30
   >>> dim=3
   >>> spacevec=DasDennis(p, dim)
   >>> pd.DataFrame(spacevec.T).to_excel("test.xlsx")


.. py:class:: ternary(*args, **kwargs)


   Bases: :py:obj:`mpltern.ternary.TernaryAxes`

   A ternary graph projection, where the input dimensions are *t*, *l*, *r*.
   The plot starts from the top and goes anti-clockwise.

   .. py:attribute:: name
      :value: 'tern'

      

   .. py:method:: set_labels(label='mass fractions / -', title='T = 298.15 K \np = 1 bar', xlabel='solvent', ylabel='polymer', zlabel='API')


   .. py:method:: filled_line(x, y, z, Formatstring, legend)


   .. py:method:: conodes(RBx, RBy, RBz, LBx, LBy, LBz, Formatstring, legend)



.. py:function:: get_line_from_tern(cs)


