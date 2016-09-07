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
            self.draw_config = {
                    'force': {'color':'green'},
                    'background': {'color':'black'},
                    'grid': {'color':'white'},
                    'element': {'color':'yellow', 'linewidth':1},
                    'support': {'markeredgecolor':'cyan',
                                'markerfacecolor':'None',
                                'ms':25},
                    'node': {'color':'yellow'},
                    'member force positive': {'edgecolor':'white',
                                              'facecolor':'blue',
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'white',
                                              'facecolor':'red',
                                              'alpha':0.5},
                    }
        elif theme == 'publication':
            self.draw_config = {
                    'force': {'color':'black'},
                    'background': {'color':'white'},
                    'grid': {'color':'black'},
                    'element': {'color':'black', 'linewidth':2},
                    'support': {'markeredgecolor':'black',
                                'markerfacecolor':'None',
                                'markeredgewidth':2,
                                'ms':25},
                    'node': {'color':'black'},
                    'member force positive': {'edgecolor':'black',
                                              'facecolor':'None',
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'black',
                                              'facecolor':'None',
                                              'alpha':0.5},
                    }
        else:
            self.draw_config = {
                    'force': {'color':'red'},
                    'background': {'color':'white'},
                    'grid': {'color':'black'},
                    'element': {'color':'orange', 'linewidth':1},
                    'support': {'markeredgecolor':'cyan',
                                'markerfacecolor':'None',
                                'ms':25},
                    'node': {'color':'orange'},
                    'member force positive': {'edgecolor':'black',
                                              'facecolor':'blue',
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'black',
                                              'facecolor':'red',
                                              'alpha':0.5},
                    }

    def plot_geometry(self, ax):
        """Plots the geometry of the model passed

        :model: TODO
        :ax: a matplotlib axis object
        :returns: TODO

        """
        node_options = self.draw_config['node']
        elem_options = self.draw_config['element']
        background_options = self.draw_config['background']
        grid_options = self.draw_config['grid']
        force_options = self.draw_config['force']
        support_options = self.draw_config['support']

        # set background color
        ax.set_axis_bgcolor(background_options['color'])
        ax.grid(True, **grid_options)

        model = self._model
        nodes = model.nodes
        if self.display_config['nodes'] == True:
            for num, node in nodes.items():
                ax.scatter(node.x, node.y, marker='o', **node_options)

        for num, elem in model.beams.items():
            n1 = elem._node1
            n2 = elem._node2
            ax.plot([n1.x, n2.x], [n1.y, n2.y], **elem_options)


        # Plot forces if requiered
        if self.display_config['forces']==True:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i._Loads.items():
                    ax = self.plot_nodal_force(ax, dof, at=node_i, val=val)

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

    def plot_internal_forces(self, ax, result, component, scale=1):
        """Plot the diagrams of internal forces

        :ax: matplotlib axis
        :result: Result object
        :component: internal force
            - 'axial'
            - 'shear'
            - 'moment'
        :scale: scale for the member forces (the member forces are already
        automatically scaled to fit in the display. This is a scale factor
        that multiplies that automatically calulated factor)

        :returns: matplotlib axis

        """
        member_force_options_pos = self.draw_config['member force positive']
        member_force_options_neg = self.draw_config['member force negative']
        model = self._model
        # First plot the geometry
        ax = self.plot_geometry(ax)
        # auxiliary axes
        fig_aux = plt.figure()
        ax_aux = fig_aux.add_subplot(111)
        # Determine automatic scaling factor
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        min_range = min([x_range, y_range])
        scale_auto = 0.2*min_range/result._max_member_force[component]
        # Plot the specified diagram
        # Transform the results from the local coordinates to global
        # coordinates
        for num, elem in model.beams.items():
            # Get the transformation matrix
            T = elem.transformation_matrix[0:2,0:2]
            # Get the specified member forces
            x_axis = result.internal_forces[num]['x']
            d = result.internal_forces[num][component]*scale_auto*scale

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
                poly_pos = Polygon(rot_pos, True, **member_force_options_pos)
                ax.add_patch(poly_pos)

            # negative part
            for curr_pol in d_neg:
                rot_neg = np.dot(curr_pol.vertices, T)
                rot_neg[:,0] += elem._node1.x
                rot_neg[:,1] += elem._node1.y
                poly_neg = Polygon(rot_neg, True, **member_force_options_neg)
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
        support_options = self.draw_config['support']
        s_size = 25
        # transform to list
        dof = list(dof)
        #

        if len(dof) == 1:
            # rolling support free in 'y'
            if dof[0] == 0:
                ax.plot([at.x],[at.y], marker=roll_y, **support_options)
            # rolling support free in 'x'
            elif dof[0] == 1:
                ax.plot([at.x],[at.y], marker=roll_x, **support_options)
            # only rotation constrained
            elif dof[0] == 2:
                ax.plot([at.x],[at.y], marker=rot_z, **support_options)
        #
        if len(dof) == 2:
            # pinned support
            if np.all([0, 1] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=pinned, **support_options)
            # 
            elif np.all([0, 2] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=disp_y, **support_options)
            #
            elif np.all([1, 2] == np.sort(dof)):
                ax.plot([at.x],[at.y], marker=disp_x, **support_options)
        #
        if len(dof) == 3:
            # Encastrated
            ax.plot([at.x],[at.y], marker=encas, **support_options)

        return ax

    def plot_nodal_force(self, ax, dof, at, val):
        """TODO: Docstring for plot_forces.

        :ax: matplotlib axis
        :dof: Degree of freedom to which the force is being applied
        :at: node at which the load is applied
        :val: value of the force
        :returns: TODO

        """
        # Get the draw options for the forces
        force_options = self.draw_config['force']

        # Force in x direction
        if dof==0:
            if val < 0:
                halign = 'left'
            else:
                halign = 'right'

            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                xytext=(-np.sign(val)*50,0), color=force_options['color'], ha=halign,
                va='center',
                textcoords='offset points',
                arrowprops = dict(arrowstyle='->', color=force_options['color'], lw=2.5 ))
        # Force in y direction
        elif dof==1:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                xytext=(0,-np.sign(val)*50), color=force_options['color'], ha='center',
                va='center',
                textcoords='offset points',
                arrowprops = dict(arrowstyle='->', color=force_options['color'], lw=2.5 ))

        return ax

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
############################################################
