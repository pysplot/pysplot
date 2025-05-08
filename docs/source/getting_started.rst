Getting Started
========================

The **pysplot** package can be easily installed using `pip`. You can install the latest stable release from PyPI, or install the development version directly from GitHub.

Installation
-----------------

To install the latest released version of ``pysplot`` from `PyPI <https://pypi.org/project/pysplot/>`_, run:

.. code-block:: bash

    pip install pysplot

This will install pysplot along with its core dependencies.

To include additional dependencies when running through our examples, run:

.. code-block:: bash

    pip install pysplot["examples"]

or when using the WHPI tool package, run:

.. code-block:: bash

    pip install pysplot["whpi"]

Requirements
============

``pysplot`` requires the following packages:

- **Python** ≥ 3.9  
- **SunPy** ≥ 5.1.5 — for heliophysics-specific tools  
- **Pandas** ≥ 2.2.3 — for data manipulation  
- **Matplotlib** ≥ 3.9.0 — for plotting and visualization

Additional dependencies are installed automatically when using the example notebooks:

- **cdasws** ≥ 1.7.42 — for streaming data from CDAWeb  
- **cdflib** ≥ 0.4.9 — for reading CDF files  
- **xarray** ≥ 0.20.1 — for handling multi-dimensional arrays

or when using the WHPI tooling:

- **imageio** ≥ 2.9.0 — for image manipulation