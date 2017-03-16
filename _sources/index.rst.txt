.. Sajou documentation master file, created by
   sphinx-quickstart on Sun Feb 26 15:52:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Sajou
=====

Introduction
------------
Sajou is a simple Structural Analysis module for Python, that enables the static computation of 2D structures.
It implements the finite element method with a *civil engineering* oriented approach, that allows the obtention of relevant information,
like moments and axial forces of beams, in an easy manner.

Although currently only 2D models are implemented, the code has been designed with 3D support in mind, and will be implemented at
some point in the future.
In a similar way, the element library should be expanded in the future, implementing diverse and useful elements like plates or shells.

.. warning:: Sajou is currently in a very early stage of development. As such, expect major bugs, radical
    changes in the API, inexplicably burn down of your computer and sudden deaths of kittens. Use under your own
    responsability.

Contents
--------

.. toctree::
    :maxdepth: 1

    install
    quickstart
    tutorial/index
    doc/index
    apiref/index


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


