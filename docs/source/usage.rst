Usage
==========

This page demonstrates how to use **pysplot** for basic plotting tasks. If you haven't already installed the package, see the :doc:`getting_started` section.

Importing pysplot
----------------------------------

To get started, first import the library:

.. code-block:: python

    import pysplot

Basic Example
---------------

The most common use of ``pysplot`` involves initializing spatial and science timeseries data as data objects. These data objects will make it easier to standardize time cadence, units, and coordinate frames. Then, these data objects can be merged into one combined data object for easy plotting.

.. code-block:: python

    from pysplot.io.data import SpatialData, ScienceData, SpatialTimeData
    from pysplot.plots.plottypes import orbit_plot
    import matplotlib.pyplot as plt

    # input data as dictionaries
    spatial_data_dictionary = {
        'x': ['2020-01-01 01:00', '2020-01-01 02:00', '2020-01-01 03:00'], 
        'y': [[0,1,2],[3,4,5],[6,7,8]]
    }
    science_data_dictionary = {
        'x': ['2020-01-01 02:03', '2020-01-01 02:30', '2020-01-01 02:45'], 
        'y': [[1.08, 14.14, -16.24], [1.39, 14.08, -17.80], [1.86, 13.56, -18.79]]
    }

    # Initiate data objects
    spatial_data = SpatialData(spatial_data_dictionary, spatial_columns_names=['x','y','z'])
    science_data = ScienceData(science_data_dictionary, science_columns_names=['Bx', 'By', 'Bz'])
    combined_data = SpatialTimeData(spatial_data, science=science_data)

    # Make plots
    fig, ax = plt.subplots()
    orbit_plot(combined_data.data, ax, 'x', 'y')


For a full example of how to use ``pysplot``, refer to our `Getting Started notebook <../examples/getting_started.ipynb>`.


Using WHPI tools
------------------------------

``Pysplot`` also includes tools specific to the Whole Heliosphere and Planetary Interactions effort, specifically downloading and plotting McIntosh Carrington Synoptic Maps from a publicly available NOAA repository. See `here <https://whpi.hao.ucar.edu/whpi_mcintosh_maps.php>`_ for more details. for more details.

To use, make sure your install of ``pysplot`` includes the ``whpi`` requirements:

.. code-block:: bash

    pip install pysplot["whpi"]

.. code-block:: python

    from pysplot.util.whpi import get_chmap_image
    import matplotlib.pyplot as plt

    chmap = get_chmap_image('2019-07-26')
    fig, ax = plt.subplots()
    ax.imshow(chmap, aspect='auto')
    ax.set_axis_off()
    plt.show()