#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This fil defines all the functions needed to plot the created structures and the results
# obtained with the program.
import numpy as np
import matplotlib as pyplot
from matplotlib.patches import Polygon

def plot_geometry(model, ax):
    """Plots the geometry of the model passed

    :model: TODO
    :ax: a matplotlib axis object
    :returns: TODO

    """
    nodes = model.nodes
    for num, node in nodes.items():
        ax.scatter(node.x, node.y, marker='o', color='b')

    for num, elem in model.beams.items():
        n1 = elem._node1
        n2 = elem._node2
        ax.plot([n1.x, n2.x], [n1.y, n2.y], color='b')

    ax.axis('equal')

    return ax

def plot_internal_forces(model, ax, result, component, scale):
    """Plot the diagrams of internal forces

    :res: three options:
        - 'axial'
        - 'moment'
        - 'shear'

    :returns: matplotlib axis

    """
    # First plot the geometry
    ax = plot_geometry(model, ax)
    # Plot the specified diagram
    # Transform the results from the local coordinates to global
    # coordinates
    for num, elem in model.beams.items():
        T = elem.transformation_matrix[0:2,0:2]
        x_axis = result.internal_forces[num]['x']
        axial = result.internal_forces[num][component]*scale
        diag = np.vstack([x_axis, axial])
        diag = np.dot(diag.T, T)

        x = np.hstack([elem._node1.x, diag.T[0] + elem._node1.x, elem._node2.x])
        d = np.hstack([elem._node1.y, diag.T[1] + elem._node1.y, elem._node2.y])
        polygon = Polygon(np.vstack([x,d]).T, True, facecolor='red', alpha=0.5)
        #ax.plot(x, d)
        ax.add_patch(polygon)

    return ax
