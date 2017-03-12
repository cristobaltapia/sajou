Quick start
===========

.. currentmodule:: sajou

Loading the library
-------------------

To use **sajou** simply import the library as you would usually do::

    import sajou as sj

That's it! After this, you are ready to start building your model.

Building the model
------------------

A simple frame structure as described in the figure below will be calculated.

.. plot:: sajou_examples/quickstart_geom.py

Geometry
********

To build the model a :py:class:`model.Model` has to be created::

    # Initialize a Model instance of a 2D model
    m = sj.model(name='Model 1', dimensionality='2D')

The geometry of the problem can then be defined by means of :py:class:`nodes.Node`, which is conveniently wrapped in the method :py:meth:`model.Model2D.beam` of the class :py:class:`model.Model`::

    # Add nodes in the specified positions
    n1 = m.node(0, 0)
    n2 = m.node(0, 2000)
    n3 = m.node(2000, 2000)
    n4 = m.node(2000, 0)

Beam elements (:py:class:`elements.beam2d.Beam2D`) are created using a method of the class :py:class:`model.Model`::

    # Add Beam2D elements
    b1 = m.beam(n1, n2)
    b2 = m.beam(n2, n3)
    b3 = m.beam(n3, n4)

Material and cross-section
**************************

For this example, the material consist of a steel, with an modulus of elasticity (MOE) equals to 200GPa.
The cross-section of the beam is defined in a *general* way by giving both the area and moment of inertia
as a tuple (see :meth:`model.Model.beam_section` to understand how to pass different parameters).

+----------+---------+-------+
| Property | value   | units |
+==========+=========+=======+
| MOE      | 200     | GPa   |
+----------+---------+-------+
| Area     | 3.2e3   | mm^2  |
+----------+---------+-------+
| Inertia  | 2.356e7 | mm^4  |
+----------+---------+-------+

The material is defined by means of the :meth:`model.Model.material` method, which creates an instance of the class :class:`materials.Material`.
This is then assigned to a Section instance (:class:`sections.BeamSection`), using the :meth:`model.Model.beam_section` method and giving the
parameter ``type='general'``.::

    # Create a material instance
    steel = m.material(name='steel', data=(200e3, ), type='isotropic')
    # Create a beam section
    section_1 = m.beam_section(name='section 1', material=steel,
                               data=(32e2, 2.356e7), type='general')

The above created **Section** now needs to be assigned to a **Beam** instance (:py:class:`elements.beam2d.Beam2D`)::


    # Assign the section to the beam elements
    b1.assign_section(section_1)
    b2.assign_section(section_1)
    b3.assign_section(section_1)


Applying loads and border conditions
************************************


PyBar supports the application of both point loads as well as distributed loads. For this, the methods :meth:`model.Model.load` and :meth:`elements.beam2d.distributed_load` are used.
The border conditions (BCs) are defined with the method :meth:`elements.beam2d.distributed_load`::

    # Add border conditions
    m.bc(node=n1, v1=0., v2=0., r3=0.)
    m.bc(node=n4, v1=0., v2=0.)

    # Add point Load
    m.load(node=n2, f2=-5e3)
    m.load(node=n3, f1=-10e3)

    # Add distributed load
    b1.distributed_load(p1=-2, direction='y', coord_system='local')


Visualizing the model
*********************

A visual inspection of the model is crucial to easily spot problems in the model.
To see the current state of the model a :class:`plot.Display` instance has to be instantiated and a `Matplotlib <http://www.matplotlib.org>`_ axis has to be passed (this might change in the future)::

    # create matplotlib figure and axes
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6.,5.3.))
    ax = fig.add_subplot(111)

    # Instatiate a Display object
    disp = sj.Display(theme='dark')

    # plot the current state of the model
    ax = disp.plot_geometry(model=m, ax=ax)
    plt.show()

A figure similar to the one shown below should be created.

.. plot:: sajou_examples/quickstart_loads.py

Solving the system
******************

For this example, the static solver (:class:`solvers.StaticSolver`) is used::

    # instance of StaticSolver
    from sj.solvers import StaticSolver
    # Define output variables
    output = ['nodal displacements', 'internal forces', 'end forces']
    # Create the StaticSolver instance
    solver = StaticSolver(model=m, output=output)
    # Solve the system
    res = solver.solve()

After this, a :class:`solvers.Result` object is created (stored as `res` in the example), which contains the required results of the system.

Postprocessing
**************

The previously obtained :class:`solvers.Result` object is then used for the post-processing of the model.
This is done by means of a :class:`postprocessing.Postprocess` object, which can be used to obtain values
of section forces in a specified element or to plot the results using the above mentioned :class:`plot.Display` class.

Let us begin by extracting values o the moment, shear and axial force along a specified beam (say beam No. 2)::

    # Postprocess the results
    post = Postprocess(result=res)
    # get the values at the center of the beam
    m_0 = post.calc_moment_at(pos=0.5, element=b2, unit_length=True)
    s_0 = post.calc_shear_at(pos=0.5, element=b2, unit_length=True)
    a_0 = post.calc_axial_at(pos=0.5, element=b2, unit_length=True)

As can be seen in the code above, the option ``unit_length`` is set to ``True``, which indicates that the values given
for the ``pos`` parameter must be in the range ``[0, 1]``.

Finally, nice plots can be obtained for the different section forces (moment, shear and axial force)::

    # create the matplotlib figures
    fig1 = plt.figure(figsize=(6.5, 5.5))
    fig2 = plt.figure(figsize=(6.5, 5.5))
    fig3 = plt.figure(figsize=(6.5, 5.5))
    fig4 = plt.figure(figsize=(6.5, 5.5))

    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax3 = fig3.add_subplot(111)
    ax4 = fig4.add_subplot(111)

    # plot the moment along the frame elements
    ax1 = disp.plot_internal_forces(ax=ax1, result=res, component='moment')
    # plot the shear force the frame elements
    ax2 = disp.plot_internal_forces(ax=ax2, result=res, component='shear')
    # plot the axial force along the frame elements
    ax3 = disp.plot_internal_forces(ax=ax3, result=res, component='axial')
    # plot the deformed shape of the structure
    ax4 = disp.plot_deformed_geometry(ax=ax4, result=res, show_undeformed=True,
                                    scale=500)

.. plot:: sajou_examples/quickstart_post.py


