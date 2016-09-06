#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This fil defines all the functions needed to plot the created structures and the results
# obtained with the program.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.axes import Axes
from matplotlib.path import Path

class Display(object):

    """Display class"""

    def __init__(self, model, width, height, theme='dark'):
        """TODO: to be defined1.

        :width: width of the window in pixels
        :height: height of the window in pixels
        :theme: TODO

        """
        self._model = model
        self._width = width
        self._height = height
        self._theme = theme

        # Dsiplay configurations
        self.display_config = {
                'forces': True,
                'reactions': True,
                'supports': True,
                'nodes': True,
                'elements': True,
                }
        # Color configurations
        if theme == 'dark':
            self.color_config = {
                    'force color': 'green',
                    'background': 'black',
                    'grid color': 'white',
                    'element color': 'yellow',
                    'support color': 'cyan',
                    'node color': 'yellow',
                    'member force border': 'white',
                    }
        else:
            self.color_config = {
                    'force color': 'green',
                    'background': 'white',
                    'grid color': 'black',
                    'element color': 'blue',
                    'support color': 'red',
                    'node color': 'blue',
                    'member force border': 'black',
                    }

    def plot_geometry(self, ax):
        """Plots the geometry of the model passed

        :model: TODO
        :ax: a matplotlib axis object
        :returns: TODO

        """
        color_n = self.color_config['node color']
        color_e = self.color_config['element color']
        color_b = self.color_config['background']
        color_g = self.color_config['grid color']
        color_f = self.color_config['force color']
        color_s = self.color_config['support color']

        # set background color
        ax.set_axis_bgcolor(color_b)
        ax.grid(True, color=color_g)

        model = self._model
        nodes = model.nodes
        for num, node in nodes.items():
            ax.scatter(node.x, node.y, marker='o', color=color_n)

        for num, elem in model.beams.items():
            n1 = elem._node1
            n2 = elem._node2
            ax.plot([n1.x, n2.x], [n1.y, n2.y], color=color_e)

        # Plot forces if requiered
        if self.display_config['forces']==True:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i._Loads.items():
                    if dof==0:
                        if val < 0:
                            halign = 'left'
                        else:
                            halign = 'right'

                        ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(node_i.x,node_i.y),
                                xytext=(-np.sign(val)*50,0), color=color_f, ha=halign,
                                va='center',
                                textcoords='offset points',
                                arrowprops = dict(arrowstyle='->', color=color_f, lw=2.5 ))
                    elif dof==1:
                        ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(node_i.x,node_i.y),
                                xytext=(0,-np.sign(val)*50), color=color_f, ha='center',
                                va='center',
                                textcoords='offset points',
                                arrowprops = dict(arrowstyle='->', color=color_f, lw=2.5 ))

        # Plot supports
        if self.display_config['supports']==True:
            for ix, node_i in model.nodes.items():
                if len(node_i._BC) > 0:
                    ax = self.plot_support(ax, dof=node_i._BC.keys(), at=node_i)

        #x_min = min([n.x for ix, n in nodes.items()])
        #x_max = max([n.x for ix, n in nodes.items()])
        #y_min = min([n.y for ix, n in nodes.items()])
        #y_max = max([n.y for ix, n in nodes.items()])
        #x_range = x_max - x_min

        #ax.set_xlim(xmin=x_min-0.2*x_range, xmax=x_max+0.2*x_range)
        #ax.set_ylim(ymin=y_min-0.2*x_range, ymax=y_max+0.2*x_range)

        ax.axis('equal')

        return ax

    def plot_internal_forces(self, ax, result, component, scale):
        """Plot the diagrams of internal forces

        :res: three options:
            - 'axial'
            - 'moment'
            - 'shear'

        :returns: matplotlib axis

        """
        color_mf = self.color_config['member force border']
        model = self._model
        # First plot the geometry
        ax = self.plot_geometry(ax)
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
                poly_pos = Polygon(rot_pos, True, facecolor='blue',
                        edgecolor=color_mf, alpha=0.5)
                ax.add_patch(poly_pos)

            # negative part
            for curr_pol in d_neg:
                rot_neg = np.dot(curr_pol.vertices, T)
                rot_neg[:,0] += elem._node1.x
                rot_neg[:,1] += elem._node1.y
                poly_neg = Polygon(rot_neg, True, facecolor='red',
                        edgecolor=color_mf, alpha=0.5)
                ax.add_patch(poly_neg)

        plt.close(fig_aux)

        return ax

    def plot_support(self, ax, dof, at):
        """Plot the respective support

        :ax: matplotlib axis
        :dof: DOF that are restrained
        :at: node at which the contraints are applied
        :returns: matplotlib Axis object

        """
        color_s = self.color_config['support color']
        s_size = 25
        # transform to list
        dof = list(dof)
        #
        marker_options = dict(
                markerfacecolor='none',
                markeredgecolor=color_s,
                ms=s_size,
                )
        ############################################################
        # define markers for the supports
        ############################################################
        # rolling x
        x = [0, 1,  -1, 0, -1.2, 1.2]
        y = [0, -1, -1, 0, -1.5, -1.5]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.MOVETO, Path.LINETO]
        roll_x = Path(xy, codes)
        # rolling y
        x = [0, -1, -1, 0, -1.5, -1.5]
        y = [0, 1,  -1, 0, 1.2, -1.2]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.MOVETO, Path.LINETO]
        roll_y = Path(xy, codes)
        # pinned
        x = [0, 1,  -1, 0, -1.2, 1.2, -1.2, -0.8, -0.8, -0.4, -0.4, 0,  0,    0.4, 0.4,  0.8, 0.8,  1.2]
        y = [0, -1, -1, 0, -1,   -1,   -1.5, -1,   -1.5, -1,   -1.5, -1, -1.5, -1,  -1.5, -1,  -1.5, -1]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.MOVETO, Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO ]
        pinned = Path(xy, codes)
        # encastrated
        x = [-1, 1, 1,  -1, -1, -1.2, -0.8, -0.8, -0.4, -0.4, 0,  0,    0.4, 0.4,  0.8, 0.8,  1.2]
        y = [0,  0, -1, -1, 0,  -1.5, -1,   -1.5, -1,   -1.5, -1, -1.5, -1,  -1.5, -1,  -1.5, -1]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
                Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO ]
        encas = Path(xy, codes)
        # no rotation and no displacement in y
        x = [-1, 1, 1,  -1, -1, -1.2, 1.2]
        y = [0,  0, -1, -1, 0,  -1.5, -1.5]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.LINETO, Path.MOVETO, Path.LINETO]
        disp_x = Path(xy, codes)
        # no rotation and no displacement in x
        x = [0, 0,  -1, -1, 0, -1.5, -1.5]
        y = [1, -1, -1, 1,  1, 1.2,  -1.2]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                Path.LINETO, Path.MOVETO, Path.LINETO]
        disp_y = Path(xy, codes)
        # only rotation constrained
        x = [-0.5, 0.5, -0.5, 0.5]
        y = [-0.5, 0.5, 0.5, -0.5]
        xy = list(zip(x, y))
        codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
        rot_z = Path(xy, codes)

        if len(dof) == 1:
            # rolling support free in 'y'
            if dof[0] == 0:
                ax.plot([at.x],[at.y], marker=roll_y, **marker_options)
            # rolling support free in 'x'
            elif dof[0] == 1:
                ax.plot([at.x],[at.y], marker=roll_x, **marker_options)
            # only rotation constrained
            elif dof[0] == 2:
                ax.plot([at.x],[at.y], marker=rot_z, **marker_options)
        #
        if len(dof) == 2:
            # pinned support
            if np.all([0, 1] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=pinned, **marker_options)
            # 
            elif np.all([0, 2] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=disp_y, **marker_options)
            #
            elif np.all([1, 2] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=disp_x, **marker_options)
        #
        if len(dof) == 3:
            # Encastrated
            ax.plot([at.x],[at.y], marker=encas, **marker_options)

        return ax
