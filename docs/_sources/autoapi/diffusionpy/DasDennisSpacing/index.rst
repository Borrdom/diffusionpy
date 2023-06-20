:py:mod:`diffusionpy.DasDennisSpacing`
======================================

.. py:module:: diffusionpy.DasDennisSpacing


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.DasDennisSpacing.DasDennis



.. py:function:: DasDennis(p, dim)

   create a equidistance n dimensional spacing which satisfies the mass balance constraint

   .. rubric:: Examples

   >>> p=30
   >>> dim=3
   >>> spacevec=DasDennis(p, dim)
   >>> pd.DataFrame(spacevec.T).to_excel("test.xlsx")


