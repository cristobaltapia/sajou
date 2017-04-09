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

To build the model a :class:`.Model` has to be created::

    # Initialize a Model instance of a 2D model
    m = sj.Model(name='Model 1', dimensionality='2D')

The geometry of the problem can then be defined by means of :class:`.Node`, which is conveniently wrapped in the method :meth:`.Model2D.node` of the class :class:`.Model`:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 16-21

Beam elements (:class:`.Beam2D`) are created using a method of the class :class:`model.Model`:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 23-27

Material and cross-section
**************************

For this example, the material consist of a timber glulam, with an modulus of elasticity (MOE) equals to 12 GPa.
The cross-section of the beam is defined as a rectangular section.

.. seealso:: :class:`sections.BeamSection` to understand how to pass different parameters

+----------+---------+-------+
| Property | value   | units |
+==========+=========+=======+
| MOE      | 12      | GPa   |
+----------+---------+-------+
| width    | 100     | mm    |
+----------+---------+-------+
| depth    | 300     | mm    |
+----------+---------+-------+

The material is defined by means of the :meth:`.Model.material` method, which creates an instance of the class :class:`.Material`.
This is then assigned to a :class:`.BeamSection` instance, using the :meth:`.Model.beam_section` method and giving the
parameter ``type='rectangular'``:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 29-34

The above created :class:`.BeamSection` now needs to be assigned to a :class:`.Beam2D` instance:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 36-40

Applying loads and border conditions
************************************

Sajou supports the application of both concentrated loads as well as distributed loads. For this, the methods :meth:`.Model.load` and :meth:`.Beam2D.distributed_load` are used.
The border conditions (BCs) are defined with the method :meth:`.Model.bc`:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 42-52

.. seealso:: Concentrated and distributed moments are also supported. See the methods :meth:`.Model.load` and :meth:`.Beam2D.distributed_moment`

End release (adding a hinge)
****************************

It is also possible to add hinges at a given node, by means of the :meth:`.Beam2D.release_end` method.
This method adds an additional degree of freedom at the respective node, effectively uncoupling the rotation from the rest of the system:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 54-55

Visualizing the model
*********************

A visual inspection of the model is crucial to easily spot problems in the model.
To see the current state of the model a :class:`.Display` instance has to be instantiated and a `Matplotlib <http://www.matplotlib.org>`_ axis has to be passed (this might change in the future):

.. literalinclude:: sajou_examples/quickstart_loads.py
   :lines: 8,55-58,69

A figure similar to the one shown below should be created.

.. plot:: sajou_examples/quickstart_loads.py

Solving the system
******************

For this example, the implemented static solver (:class:`.StaticSolver`) is used:

.. literalinclude:: sajou_examples/quickstart_post.py
   :lines: 12,57-62

After this, a :class:`.Result` object is created (stored as ``res`` in the example), which contains the required results of the system.

Postprocessing
**************

The previously obtained :class:`.Result` object is then used in the post-processing of the model.
This is done by means of a :class:`.Postprocess` object, which defines several methods to obtain values
of section forces in a specified element and to plot the results using the above mentioned :class:`.Display` class.

Let us begin by extracting values o the moment, shear and axial force along a specified beam (say beam No. 2)::

    # Postprocess the results
    post = sj.Postprocess(result=res)
    # get the values at the center of the beam
    m_0 = post.calc_moment_at(pos=0.5, element=b2, unit_length=True)
    s_0 = post.calc_shear_at(pos=0.5, element=b2, unit_length=True)
    a_0 = post.calc_axial_at(pos=0.5, element=b2, unit_length=True)

As can be seen in the code above, the option ``unit_length`` is set to ``True``, which indicates that the values given
for the ``pos`` parameter must be in the range ``[0, 1]``.

.. note:: The ``pos`` parameter of the :meth:`.calc_moment_at` also accept an array-like argument, so that the
    result can be obtained at different points over the beam element at once.
    This also holds for the :meth:`.calc_shear_at` and :meth:`.calc_axial_at` methods.

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


