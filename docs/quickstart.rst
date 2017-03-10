Quick start
===========

Loading the library
-------------------

To use **pybar** simply import the library as you would usually do::

    import pybar as pb

That's it! After this, you are ready to start building your model.

Building the model
------------------

A simple frame structure as described in the figure below will be calculated

.. plot:: pybar_examples/quickstart_geom.py

Geometry
********

To build the model a :py:class:`pybar.model.Model` has to be created::

    # Initialize a Model instance of a 2D model
    m = pb.model(name='Model 1', dimensionality='2D')

The geometry of the problem can then be defined, by means of :py:class:`pybar.nodes.Node`, which is conveniently wrapped in the method :py:meth:`pybar.model.Model2D.beam` of the class :py:class:`pybar.model.Model`::

    # Add nodes in the specified positions
    n1 = m.node(0, 0)
    n2 = m.node(0, 2000)
    n3 = m.node(2000, 2000)
    n4 = m.node(2000, 0)

Beam elements (:py:class:`pybar.elements.beam2d.Beam2D`) are created using a method of the class :py:class:`pybar.model.Model`::

    # Add Beam2D elements
    b1 = m.beam(n1, n2)
    b2 = m.beam(n2, n3)
    b3 = m.beam(n3, n4)

Material and cross-section
**************************

The material is defined by means of the :py:meth:`pybar.model.Model.material` method, which creates an instance of the class :py:class:`pybar.materials.Material`.
This is then assigned to a Section instance (:py:class:`pybar.sections.BeamSection`) using the :py:meth:`pybar.model.Model.beam_section` method::

    # Create a material instance
    steel = m.material(name='steel', data=(200e3, ), type='isotropic')
    # Create a beam section
    section_1 = m.beam_section(name='section 1', material=steel,
                               data=(32e2, 2356e4), type='general')

The above created **Section** now needs to be assigned to a **Beam** instance (:py:class:`pybar.elements.beam2d.Beam2D`)::


    # Assign the section to the beam elements
    b1.assign_section(section_1)
    b2.assign_section(section_1)
    b3.assign_section(section_1)


Applying loads and border conditions
************************************


PyBar supports the application of both point loads as well as distributed loads. For this, the methods :py:meth:`pybar.model.Model.load` and :py:meth:`py:pybar.elements.beam2d.distributed_load` are used.
The border conditions (BCs) are defined with the method :py:meth:`pybar.elements.beam2d.distributed_load`::

    # Add border conditions
    m.bc(node=n1, v1=0., v2=0., r3=0.)
    m.bc(node=n4, v1=0., v2=0.)

    # Add point Load
    m.load(node=n2, f2=-20e3)
    m.load(node=n3, f1=10e3)

    # Add distributed load
    b1.distributed_load(p1=-2, direction='y', coord_system='local')

.. plot:: pybar_examples/quickstart_loads.py

Visualizing the model
*********************

A visual inspection of the model is crucial to easily spot problems in the model.
To see the current state of the model a :py:class:`pybar.plot.Display` instance has to be instantiated and a `Link Matplotlib <http://www.matplotlib.org>`_ axis has to be passed (this might change in the future)::

    # create matplotlib figure and axes
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(20.,14.))
    ax = fig.add_subplot(111)

    # Instatiate a Display object
    disp = pb.Display(theme='dark')

    # plot the current state of the model
    ax = disp.plot_geometry(m, ax)


