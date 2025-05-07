Installation
============

The **pysplot** package can be easily installed using `pip`. You can install the latest stable release from PyPI, or install the development version directly from GitHub.

Install from PyPI
-----------------

To install the latest released version of pysplot from [PyPI](https://pypi.org/project/pysplot/), run:

.. code-block:: bash

    pip install pysplot

This will install pysplot along with its core dependencies.

Install from GitHub (Development Version)
-----------------------------------------

To install the latest development version directly from the GitHub repository:

.. code-block:: bash

    pip install git+https://github.com/pysplot/pysplot.git

This is useful if you want to try out new features or contribute to the project.

Clone and Install Locally
--------------------------

Alternatively, you can clone the repository and install it manually:

.. code-block:: bash

    git clone https://github.com/pysplot/pysplot.git
    cd pysplot
    pip install .

Dependencies
------------

pysplot requires the following core Python packages:

- `matplotlib`
- `numpy`

These will be installed automatically if you're using `pip`.

Optional Development Tools
--------------------------

If you're contributing to the project or working on documentation/testing, you may want to install additional dependencies:

.. code-block:: bash

    pip install -r requirements-dev.txt

This will include linters, testing frameworks, and documentation tools like Sphinx.

Python Version
--------------

pysplot supports **Python 3.8 and above**.

For more details, visit the project repository:

https://github.com/pysplot/pysplot