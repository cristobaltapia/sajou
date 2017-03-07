Quick start
===========

Loading the library
-------------------

To use **pybar** all you have to do is import it, as usual::

    import pybar as pb

After this, you are ready to start building your model.

Building the model
------------------

To build the model a :py:class:`pybar.model.Model` has to be created::

    # Initialize a Model instance of a 2D model
    m = pb.model(name='Model 1', dimensionality='2D')

The geometry of the problem can then be defined, by means of :py:class:`pybar.nodes.Node`, which is conveniently wrapped in the method :py:meth:`pybar.model.Model2D.Beam` of the class :py:class:`pybar.model.Model`::

    # Add nodes in the specified positions
    n1 = m.node(0,0)
    n2 = m.node(0,10)
    n3 = m.node(10,10)
    n4 = m.node(10,0)

Beam elements (:py:class:`pybar.elements.beam2d.Beam2D`) are created using a method of the class :py:class:`pybar.model.Model`::

    # Add Beam2D elements
    b1 = m.beam(n1, n2)
    b2 = m.beam(n2, n3)
    b3 = m.beam(n3, n4)

