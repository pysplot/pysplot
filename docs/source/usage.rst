Usage
=====

This page demonstrates how to use **pysplot** for basic plotting tasks. If you haven't already installed the package, see the :doc:`installation` section.

Importing pysplot
-----------------

To get started, first import the library:

.. code-block:: python

    import pysplot

Basic Line Plot
---------------

The most common use of pysplot is creating line plots from data:

.. code-block:: python

    import pysplot

    x = [0, 1, 2, 3, 4]
    y = [0, 1, 4, 9, 16]

    pysplot.plot(x, y)
    pysplot.show()

This will generate a simple line plot with default settings.

Adding Labels and Title
-----------------------

You can easily add axis labels and a title:

.. code-block:: python

    pysplot.xlabel("X Axis")
    pysplot.ylabel("Y Axis")
    pysplot.title("Sample Plot")

Plotting Multiple Lines
-----------------------

pysplot also supports plotting multiple datasets:

.. code-block:: python

    y2 = [0, 1, 8, 27, 64]

    pysplot.plot(x, y, label="x squared")
    pysplot.plot(x, y2, label="x cubed")
    pysplot.legend()
    pysplot.show()

Saving a Plot to File
---------------------

You can save the current plot to a file using:

.. code-block:: python

    pysplot.save("my_plot.png")

This will save the plot in PNG format to the current working directory.

More Options
------------

pysplot provides various customization options like:

- Changing line styles and colors
- Adding gridlines
- Configuring plot dimensions
- Using subplots

Explore these features in the :doc:`examples` and :doc:`api` sections.