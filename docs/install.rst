Installation instructions
=========================

Sajou only works with Python 3 and will not be made compatible with Python 2.
The code is hosted in `Github <https://github.com/cristobaltapia/sajou/>`_, where you can report bugs and request features.

Using ``pip``
-------------

On any operating system::

    pip install git+https://github.com/cristobaltapia/sajou

Dependencies
------------

Sajou depends on the following libraries:

* `Numpy <http://www.numpy.org/>`_
* `Scipy <https://www.scipy.org/>`_
* `Matplotlib <http://matplotlib.org/>`_: Used for the visualization of the geometry and results
* `Pandas <http://pandas.pydata.org/>`_
* `Seaborn <http://seaborn.pydata.org/>`_: Used for the predefined color palettes of the plots

Optional dependencies
*********************

* `Sphinx <http://www.sphinx-doc.org/en/stable/>`_: Required to build the documentation
* `Napoleon <https://github.com/rtfd/sphinx_rtd_theme>`_: Used to interpret the *numpy* style of the docstrings
* `Read the Docs Sphinx Theme <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/>`_: Used as the documentation theme in the html version
