#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This fil defines all the functions needed to plot the created structures and the results
# obtained with the program.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.axes import Axes

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
    # auxiliary axes
    fig_aux = plt.figure()
    ax_aux = fig_aux.add_subplot(111)
    # Plot the specified diagram
    # Transform the results from the local coordinates to global
    # coordinates
    for num, elem in model.beams.items():
        # Get the transformation matrix
        T = elem.transformation_matrix[0:2,0:2]
        # Get the specified member forces
        x_axis = result.internal_forces[num]['x']
        d = result.internal_forces[num][component]*scale

        # Positive values
        d_pos = ax_aux.fill_between(x_axis, d, 0, where=d>0, interpolate=True)
        d_pos = d_pos.get_paths()
        # Negative values
        d_neg = ax_aux.fill_between(x_axis, d, 0, where=d<0, interpolate=True)
        d_neg = d_neg.get_paths()

        # Plot the patches
        # positive part
        for curr_pol in d_pos:
            rot_pos = np.dot(curr_pol.vertices, T)
            rot_pos[:,0] += elem._node1.x
            rot_pos[:,1] += elem._node1.y
            poly_pos = Polygon(rot_pos, True, facecolor='blue', alpha=0.5)
            ax.add_patch(poly_pos)

        # negative part
        for curr_pol in d_neg:
            rot_neg = np.dot(curr_pol.vertices, T)
            rot_neg[:,0] += elem._node1.x
            rot_neg[:,1] += elem._node1.y
            poly_neg = Polygon(rot_neg, True, facecolor='red', alpha=0.5)
            ax.add_patch(poly_neg)

    plt.close(fig_aux)

    return ax
