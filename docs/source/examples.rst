Examples
========

This section provides advanced examples to help you get the most out of **pysplot**. For basic usage, see the :doc:`usage` page.

Multiple Line Plot
------------------

You can plot multiple lines on the same graph with labels:

.. code-block:: python

    import pysplot

    x = [0, 1, 2, 3, 4]
    y1 = [0, 1, 4, 9, 16]
    y2 = [0, 1, 8, 27, 64]

    pysplot.plot(x, y1, label='x squared', color='blue')
    pysplot.plot(x, y2, label='x cubed', color='red')
    pysplot.xlabel("X Axis")
    pysplot.ylabel("Y Axis")
    pysplot.title("Comparing Functions")
    pysplot.legend()
    pysplot.grid(True)
    pysplot.show()

Subplots Example
----------------

Use subplots to show multiple plots in one figure:

.. code-block:: python

    import pysplot

    x = [0, 1, 2, 3, 4]
    y1 = [i ** 2 for i in x]
    y2 = [i ** 3 for i in x]

    fig, axs = pysplot.subplots(1, 2, figsize=(10, 4))

    axs[0].plot(x, y1, label='x squared')
    axs[0].set_title("x^2")

    axs[1].plot(x, y2, label='x cubed', color='orange')
    axs[1].set_title("x^3")

    for ax in axs:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.grid(True)

    pysplot.tight_layout()
    pysplot.show()

Customizing Styles
------------------

pysplot allows you to apply custom styles for cleaner visuals:

.. code-block:: python

    import pysplot

    pysplot.style('ggplot')  # Try 'seaborn', 'classic', etc.

    x = range(10)
    y = [i ** 2 for i in x]

    pysplot.plot(x, y, marker='o', linestyle='--')
    pysplot.title("Styled Plot")
    pysplot.show()

Saving Plots
------------

Save a plot to an image file instead of displaying it:

.. code-block:: python

    pysplot.plot([1, 2, 3], [4, 5, 6])
    pysplot.title("Save Example")
    pysplot.save("output_plot.png")