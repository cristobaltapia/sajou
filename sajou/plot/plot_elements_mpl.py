#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines the functions used to plot each element type in the matplotlib backend
"""


def plot_element(element, ax, elem_options):
    """Plot the given element

    This is a general function that calls the specific plotting function for
    the element passed.

    Parameters
    ----------

    element: Element
        the element that will be plotted

    Returns
    -------

    matplotlib Axis:
        axis containing the plotted element

    """
    if element._kind == 'Beam2D':
        return plot_beam2d(element, ax, elem_options)
    if element._kind == 'Spring2D':
        return plot_spring2d(element, ax, elem_options)


def plot_beam2d(element, ax, elem_options):
    """Plots the given Beam2D element

    element: Beam2D instance
        the element to be plotted
    :returns: TODO

    """
    # get the nodes
    n1 = element._node1
    n2 = element._node2
    # plot the beam element
    ax.plot([n1.x, n2.x], [n1.y, n2.y], **elem_options)

    return ax


def plot_deformed_beam2d(element, results, ax, elem_options, scale=1):
    """Plot the deformed shape of the Beam2D element

    Parameters
    ----------

    element: Beam2D instance
        element to be plotted
    results: Results instance
        result obtained after solving the system
    ax: matplotlib Axis
        axis where to plot the element to
    elem_options: dict
        options to pass to the 'plot' function
    scale: float
        scale applied to the displacements

    Returns
    -------

    matplotlib Axis:
        axis where the element was plotted

    """
    # get node coordinates and displacements
    coord_nodes, displ_nodes = self._get_deformed_node_coordinates(
        elem, model, result)
    # calculate the deformed position (including the scaling factor)
    deformed_nodes = coord_nodes + (scale * displ_nodes)
    # Plot deformed element
    ax.plot(deformed_nodes['x'], deformed_nodes['y'], **elem_options)

    return ax


def plot_spring2d(element, ax, elem_options):
    """Plots the given Beam2D element

    element: Beam2D instance
        the element to be plotted
    :returns: TODO

    """
    # get the nodes
    n1 = element._node1
    n2 = element._node2
    # plot the beam element
    ax.plot([n1.x, n2.x], [n1.y, n2.y], **elem_options)

    return ax
