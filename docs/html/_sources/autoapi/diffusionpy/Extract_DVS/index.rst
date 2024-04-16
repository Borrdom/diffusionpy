:py:mod:`diffusionpy.Extract_DVS`
=================================

.. py:module:: diffusionpy.Extract_DVS


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   diffusionpy.Extract_DVS.Measurement
   diffusionpy.Extract_DVS.Checkbar
   diffusionpy.Extract_DVS.window3



Functions
~~~~~~~~~

.. autoapisummary::

   diffusionpy.Extract_DVS.averagelist
   diffusionpy.Extract_DVS.averagelisttimeseries
   diffusionpy.Extract_DVS.averagedict
   diffusionpy.Extract_DVS.unnesting



Attributes
~~~~~~~~~~

.. autoapisummary::

   diffusionpy.Extract_DVS.cwd
   diffusionpy.Extract_DVS.filenames


.. py:data:: cwd

   

.. py:function:: averagelist(lis)


.. py:function:: averagelisttimeseries(t, x)


.. py:function:: averagedict(x)


.. py:function:: unnesting(df, explode, axis)


.. py:class:: Measurement(DATA_FILE, desorption=False, actual=False, tend=-1, arc=True, nt=15)


   The class Measurement inherits methods to read and process the DVS excel files from the Apparatus of Surface Measurement Systems.
   The endpoints of each RH step and their kinetics are extracted and the crank equation is fitted to the kinetics.
   The extracted data is added via two new sheets to all selected .xls files. The reading functions are specialized for the format of the .xls files.
   There is little robustness so large changes cause the methids not to work .
   The class accepts more than one excel file and then performs averaging of the kinetics. For that purpose, it is assumed that the all excel sheets are
   replicates of each other, meaning that their RH steps must match

   .. py:method:: FasterXlrdRead(filename)
      :staticmethod:


   .. py:method:: read_excel_file(filename)


   .. py:method:: FitDiffusionData()


   .. py:method:: Cranc()


   .. py:method:: interparc(x, y, N)
      :staticmethod:

      function that interpolates between a given data series to provide datapoints that are equidistant along the arc of the data.
      Hence the name interp(olate)arc. This is quite handy for kinetic data as the most change in concentration is at earlier times and the least
      change is observed at later times. As a result, usually more data points are in the later stages where nothing happens, since measurements are
      usually performed at equidistant time point



.. py:class:: Checkbar(parent=None, picks=[], commands=lambda: None, side='top', anchor='w')


   Bases: :py:obj:`tkinter.Frame`

   Frame widget which may contain other widgets and can have a 3D border.

   .. py:method:: state()



.. py:class:: window3(Film1)


   .. py:method:: enter()


   .. py:method:: indexjump(indJ)


   .. py:method:: indexFile(indF)


   .. py:method:: PlotRHSteps()


   .. py:method:: Desorption()


   .. py:method:: main()

      The main function initializes the GUI-window.



.. py:data:: filenames

   

