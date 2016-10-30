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
import seaborn.apionly as sns

class Display(object):

    """Display class"""

    def __init__(self, width, height, theme='dark', backend='matplotlib'):
        """TODO: to be defined1.

        :width: width of the window in pixels
        :height: height of the window in pixels
        :theme: TODO
        :backend: backend used to display the results
            - 'matplotlib'
            - 'python-occ' (not implemented yet)
            - 'mayavi' (not implemented yet)

        """
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
            palette = sns.color_palette('bright')
            self.draw_config = {
                    'force': {'color':palette[1]},
                    'reaction': {'color':'yellow'},
                    'background': {'color':"#12191f"},
                    'grid': {'color':'white'},
                    'element': {'color':palette[4], 'linewidth':1},
                    'support': {'markeredgecolor':palette[5],
                                'markerfacecolor':'None',
                                'ms':25},
                    'node': {'color':palette[4]},
                    'member force positive': {'edgecolor':'white',
                                              'facecolor':palette[0],
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'white',
                                              'facecolor':palette[2],
                                              'alpha':0.5},
                    }
        elif theme == 'publication':
            self.draw_config = {
                    'force': {'color':'black'},
                    'reaction': {'color':'black'},
                    'background': {'color':'white'},
                    'grid': {'color':'black'},
                    'element': {'color':'black', 'linewidth':2},
                    'support': {'markeredgecolor':'black',
                                'markerfacecolor':'None',
                                'markeredgewidth':2,
                                'ms':35},
                    'node': {'color':'black'},
                    'member force positive': {'edgecolor':'black',
                                              'facecolor':'None',
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'black',
                                              'facecolor':'None',
                                              'alpha':0.5},
                    }
        elif theme == 'light':
            palette = sns.color_palette('dark')
            self.draw_config = {
                    'force': {'color':palette[1]},
                    'reaction': {'color':palette[3]},
                    'background': {'color':"#dee5ec"},
                    'grid': {'color':'black'},
                    'element': {'color':palette[4], 'linewidth':1},
                    'support': {'markeredgecolor':palette[5],
                                'markerfacecolor':'None',
                                'ms':25},
                    'node': {'color':palette[4]},
                    'member force positive': {'edgecolor':'black',
                                              'facecolor':palette[0],
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'black',
                                              'facecolor':palette[2],
                                              'alpha':0.5},
                    }
        else:
            palette = sns.color_palette('dark')
            self.draw_config = {
                    'force': {'color':palette[1]},
                    'reaction': {'color':palette[3]},
                    'background': {'color':"#dee5ec"},
                    'grid': {'color':'black'},
                    'element': {'color':palette[4], 'linewidth':1},
                    'support': {'markeredgecolor':palette[5],
                                'markerfacecolor':'None',
                                'ms':25},
                    'node': {'color':palette[4]},
                    'member force positive': {'edgecolor':'black',
                                              'facecolor':palette[0],
                                              'alpha':0.5},
                    'member force negative': {'edgecolor':'black',
                                              'facecolor':palette[2],
                                              'alpha':0.5},
                    }

    def plot_geometry(self, model, ax, ls='-'):
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

        for num, elem in model.beams.items():
            n1 = elem._node1
            n2 = elem._node2
            ax.plot([n1.x, n2.x], [n1.y, n2.y], **elem_options, ls=ls)

        nodes = model.nodes
        if self.display_config['nodes'] == True:
            for num, node in nodes.items():
                ax.scatter(node.x, node.y, marker='o', **node_options)

        # Plot forces if requiered
        if self.display_config['forces']==True:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i._Loads.items():
                    ax = self.plot_nodal_force(ax, dof, at=node_i, val=val)
            # Plot element Loads
            ax = self.plot_element_loads(ax, model)

        # Plot supports
        if self.display_config['supports']==True:
            for ix, node_i in model.nodes.items():
                if len(node_i._BC) > 0:
                    ax = self.plot_support(ax, dof=node_i._BC.keys(), at=node_i)

        ax.axis('equal')

        return ax

    def plot_deformed_geometry(self, result, ax, show_undeformed=False, scale=1.):
        """Plot the system in its deformed configuration.

        :result:
        :ax: matplotlib axis instance
        :scale: scale used to plot the deformations
        :returns: matplotlib axis

        """
        from .model import get_dataframe_of_node_coords
        #
        node_options = self.draw_config['node']
        elem_options = self.draw_config['element']
        background_options = self.draw_config['background']
        grid_options = self.draw_config['grid']
        force_options = self.draw_config['force']
        support_options = self.draw_config['support']
        #
        model = result._model

        # set background color
        ax.set_axis_bgcolor(background_options['color'])
        ax.grid(True, **grid_options)

        for num, elem in model.beams.items():
            n1 = elem._node1.number
            n2 = elem._node2.number
            # Get original coordinates of nodes
            coords_nodes = get_dataframe_of_node_coords(model=model, nodes=[n1, n2])
            # Get displacements of the nodes
            displ_nodes = result.data['nodal displacements'].loc[[n1, n2]]
            # Calculate position of nodes after load application
            deformed_nodes = coords_nodes + displ_nodes * scale
            # Plot deformed element
            ax.plot(deformed_nodes['x'], deformed_nodes['y'], **elem_options)

        nodes = model.nodes
        if self.display_config['nodes'] == True:
            # Get original coordinates of nodes
            coords_nodes = get_dataframe_of_node_coords(model=model, nodes='all')
            # Get displacements of the nodes
            displ_nodes = result.data['nodal displacements']
            # Calculate position of nodes after load application
            deformed_nodes = coords_nodes + displ_nodes * scale
            # plot nodes
            ax.scatter(deformed_nodes['x'], deformed_nodes['y'], marker='o', **node_options)
        
        if show_undeformed == True:
            ax = self.plot_geometry(model=result._model, ax=ax, ls='--')

        ax.axis('equal')

        return ax

    def plot_element_loads(self, ax, model):
        """Plot the element loads.
        Element loads like uniformly distributed loads are drawn accordinlgy.

        :ax: TODO
        :result: TODO
        :returns: TODO

        """
        elems_with_loads = {ix: elem for ix, elem in model.beams.items() 
                if elem._loads != []}

        # Get maximum distributed load
        if len(elems_with_loads) > 0:
            pe = []
            for ix, elem in elems_with_loads.items():
                for load in elem._loads:
                    try:
                        p_max = np.max([np.abs(load._p1), np.abs(load._p2)])
                        pe.append(p_max)
                    except:
                        pass

            pmax = np.max(pe)
        # Determine automatic scaling factor
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        scale = max([x_range, y_range])
        # A loop is done for each element that contains an element load
        for ix, elem in elems_with_loads.items():
            # Plot each element load
            for load_i in elem._loads:
                # Different plot function is applied for the different
                # element loadins
                if load_i._type == 'Distributed Load':
                    ax = self.plot_distributed_load_element(ax, model, elem, load_i,
                            scale=scale, max_p=pmax)
                # TODO: the other types

        return ax

    def plot_distributed_load_element(self, ax, model, element, load, scale=1, max_p=1):
        """Plot a distributed load at the corresponding element.

        :ax: TODO
        :model: TODO
        :element:
        :load: TODO
        :returns: TODO

        """
        # Get the draw options for the forces
        force_options = self.draw_config['force']
        #
        size = scale * 0.1
        # If it is a uniformly distributed load:
        if load._direction == 'z' and load._coord_system == 'local':
            #
            n1 = element._node1
            n2 = element._node2
            p1 = load._p1
            p2 = load._p2
            # Number of arrows
            n_arrows = (element._length*2)//size
            # Generate points on the element line
            nx_e = np.linspace(n1.x, n2.x, n_arrows)
            ny_e = np.linspace(n1.y, n2.y, n_arrows)
            # Offset from the position of the nodes...
            # FIXME: ... according to the system coordinate chosen
            T = element.transformation_matrix
            sign_p = np.sign(load._p1)
            offset = T.T.dot(np.array([0, -p1, 0, 0, -p2, 0]) * scale / max_p * 0.2)
            n1_o = np.array(n1[:]) + offset[:3]
            n2_o = np.array(n2[:]) + offset[3:]
            # Generate points for the arrows with offset
            nx = np.linspace(n1_o[0], n2_o[0], n_arrows)
            ny = np.linspace(n1_o[1], n2_o[1], n_arrows)
            ax.plot(nx, ny, color=force_options['color'])
            # Plot arrows
            for arr_i in range(len(nx)):
                # annotate() is used instead of arrow() because the style
                # of the arrows is better
                ax.annotate('',
                        xy=(nx_e[arr_i], ny_e[arr_i]),
                        xytext=(nx[arr_i], ny[arr_i]),
                        arrowprops=dict(arrowstyle='->', color=force_options['color']))

        elif load._direction == 'x' and load._coord_system == 'local':
            #
            n1 = element._node1
            n2 = element._node2
            # Number of arrows
            n_arrows = (element._length*2)//size
            # Generate points on the element line
            nx_e = np.linspace(n1.x, n2.x, n_arrows)
            ny_e = np.linspace(n1.y, n2.y, n_arrows)
            # Offset from the position of the nodes...
            # FIXME: ... according to the system coordinate chosen
            T = element.transformation_matrix
            sign_p = np.sign(load._p1)
            offset = T.T.dot(np.array([-1,0,0, -1,0,0]) * sign_p * size*0.4)
            n1_o = np.array(n1[:]) + offset[:3]
            n2_o = np.array(n2[:]) + offset[3:]
            # Generate points for the arrows with offset
            nx = np.linspace(n1_o[0], n2_o[0], n_arrows)
            ny = np.linspace(n1_o[1], n2_o[1], n_arrows)
            # Plot arrows
            for arr_i in range(len(nx)):
                # annotate() is used instead of arrow() because the style
                # of the arrows is better
                ax.annotate('',
                            xy=(nx_e[arr_i], ny_e[arr_i]),
                            xytext=(nx[arr_i], ny[arr_i]),
                            arrowprops=dict(arrowstyle='->', color=force_options['color']))

        return ax

    def plot_internal_forces(self, ax, result, component, scale=1):
        """Plot the diagrams of internal forces.

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
        model = result._model
        # First plot the geometry
        ax = self.plot_geometry(model, ax)

        # auxiliary axes
        fig_aux = plt.figure()
        ax_aux = fig_aux.add_subplot(111)
        # Determine automatic scaling factor
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        max_range = max([x_range, y_range])
        max_component = result.metadata['internal forces']['system abs max'][component]
        scale_auto = 0.2 * max_range / max_component
        # Plot the specified diagram
        # Transform the results from the local coordinates to global
        # coordinates
        for num, elem in model.beams.items():
            # Get the transformation matrix
            T = elem.transformation_matrix[0:2,0:2]
            # Get the positions at which the internal forces were
            # caclculated
            x_axis = result.data['internal forces'][num][component]['x']
            # Get the specified member forces
            d = result.data['internal forces'][num][component]['data'] * scale_auto * scale

            # Results are mutiplied by -1 to represent them in the convention
            # where positive moments are drawn 'below' the beam element.
            d = -d
            # Positive values
            d_pos = ax_aux.fill_between(x_axis, d, 0, where=d<0, interpolate=True)
            d_pos = d_pos.get_paths()
            # Negative values
            d_neg = ax_aux.fill_between(x_axis, d, 0, where=d>0, interpolate=True)
            d_neg = d_neg.get_paths()

            # Plot the patches
            # positive part
            for curr_pol in d_pos:
                rot_pos = np.dot(curr_pol.vertices, T.todense())
                rot_pos[:,0] += elem._node1.x
                rot_pos[:,1] += elem._node1.y
                poly_pos = Polygon(rot_pos, True, **member_force_options_pos)
                ax.add_patch(poly_pos)

            # negative part
            for curr_pol in d_neg:
                rot_neg = np.dot(curr_pol.vertices, T.todense())
                rot_neg[:,0] += elem._node1.x
                rot_neg[:,1] += elem._node1.y
                poly_neg = Polygon(rot_neg, True, **member_force_options_neg)
                ax.add_patch(poly_neg)

        # Plot reaction forces
        if self.display_config['reactions']==True:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i.reactions.items():
                    ax = self.plot_nodal_reaction(ax, dof, at=node_i, val=val)

        # close the auxiliary figure
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
        # Moment
        elif dof==2:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                xytext=(20,20), color=force_options['color'], ha='left',
                va='center',
                textcoords='offset points')
            if np.sign(val) >= 0.:
                ax.plot([at.x],[at.y], marker=marker_moment_pos, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)
            else:
                ax.plot([at.x],[at.y], marker=marker_moment_neg, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)

        return ax

    def plot_nodal_reaction(self, ax, dof, at, val):
        """Plot the nodal reactions of the system

        :ax: matplotlib axis
        :dof: Degree of freedom to which the force is being applied
        :at: node at which the load is applied
        :val: value of the reaction
        :returns: TODO

        """
        # Get the draw options for the forces
        force_options = self.draw_config['reaction']

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
        # Moment
        elif dof==2:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                xytext=(20,20), color=force_options['color'], ha='left',
                va='center',
                textcoords='offset points')
            if np.sign(val) >= 0.:
                ax.plot([at.x],[at.y], marker=marker_moment_pos, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)
            else:
                ax.plot([at.x],[at.y], marker=marker_moment_neg, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)

        return ax

class DisplaySym(Display):

    """Symbolic display"""

    def __init__(self, width, height, theme='dark', backend='matplotlib'):
        """TODO: to be defined1. """
        Display.__init__(self, width, height, theme)

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
        model = result._model
        # First plot the geometry
        ax = self.plot_geometry(model, ax)

        # auxiliary axes
        fig_aux = plt.figure()
        ax_aux = fig_aux.add_subplot(111)
        # Determine automatic scaling factor
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        max_range = max([x_range, y_range])
        scale_auto = 0.2*max_range / result._max_member_force[component]

        # Plot the specified diagram
        # Transform the results from the local coordinates to global
        # coordinates
        for num, elem in model.beams.items():
            # Get the transformation matrix
            T = elem.transformation_matrix[0:2,0:2]
            # Get the specified member forces
            x_axis = result.internal_forces[num]['x']
            d = result.internal_forces[num][component]*scale_auto*scale

            # Results are mutiplied by -1 to represent them in the convention
            # where positive moments are drawn 'below' the beam element.
            d = -d
            # Positive values
            d_pos = ax_aux.fill_between(x_axis, d, 0, where=d<0, interpolate=True)
            d_pos = d_pos.get_paths()
            # Negative values
            d_neg = ax_aux.fill_between(x_axis, d, 0, where=d>0, interpolate=True)
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

        # Plot reaction forces
        if self.display_config['reactions']==True:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i.reactions.items():
                    ax = self.plot_nodal_reaction(ax, dof, at=node_i, val=val)

        # close the auxiliary figure
        plt.close(fig_aux)

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
            halign = 'left'
            # Handles the symbolic case
            ax.annotate('{f}'.format(f=val), xy=(at.x, at.y),
                xytext=(-50,0), color=force_options['color'], ha=halign,
                va='center',
                textcoords='offset points',
                arrowprops = dict(arrowstyle='->', color=force_options['color'], lw=2.5 ))
        # Force in y direction
        elif dof==1:
            ax.annotate('{f}'.format(f=-val), xy=(at.x, at.y),
                xytext=(0,50), color=force_options['color'], ha='center',
                va='center',
                textcoords='offset points',
                arrowprops = dict(arrowstyle='->', color=force_options['color'], lw=2.5 ))
        # Moment
        elif dof==2:
            ax.annotate('{f}'.format(f=val), xy=(at.x, at.y),
                xytext=(20,20), color=force_options['color'], ha='left',
                va='center',
                textcoords='offset points')
            ax.plot([at.x],[at.y], marker=marker_moment_pos, markeredgecolor=force_options['color'],
                    markerfacecolor='None', ms=40, markeredgewidth=2)


        return ax

    def plot_nodal_reaction(self, ax, dof, at, val):
        """Plot the nodal reactions of the system

        :ax: matplotlib axis
        :dof: Degree of freedom to which the force is being applied
        :at: node at which the load is applied
        :val: value of the reaction
        :returns: TODO

        """
        # Get the draw options for the forces
        force_options = self.draw_config['reaction']

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
        # Moment
        elif dof==2:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                xytext=(20,20), color=force_options['color'], ha='left',
                va='center',
                textcoords='offset points')
            if np.sign(val) >= 0.:
                ax.plot([at.x],[at.y], marker=marker_moment_pos, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)
            else:
                ax.plot([at.x],[at.y], marker=marker_moment_neg, markeredgecolor=force_options['color'],
                        markerfacecolor='None', ms=40, markeredgewidth=2)

        return ax

############################################################
# define markers for the moment application
############################################################
#
marker_moment_pos = Path.arc(theta1=-45., theta2=135.)
verts_m = marker_moment_pos.vertices
codes_m = marker_moment_pos.codes
arrow_head_s = 0.4
# get angle on the tip (approx.)
ang_aux_pos = np.arctan2((verts_m[-1,1]-verts_m[-3,1]), (verts_m[-1,0]-verts_m[-3,0]))
ang_aux_neg = np.arctan2((verts_m[0,1]-verts_m[2,1]), (verts_m[0,0]-verts_m[2,0]))
# Positive Moment
ang_arrow = np.deg2rad(25.)
ang_p1 = ang_aux_pos - np.pi + ang_arrow
ang_p2 = ang_aux_pos - np.pi - ang_arrow
# Position of the tip of the arrow
tip_pos = verts_m[-1,:]
# Position of the first line of the arrow
p1 = tip_pos + arrow_head_s*np.array([np.cos(ang_p1), np.sin(ang_p1)])
# Position of the second line of the arrow
p2 = tip_pos + arrow_head_s*np.array([np.cos(ang_p2), np.sin(ang_p2)])
# Add to the vertices and code
aux = np.array([p1, tip_pos, p2])
verts_mp = np.vstack((verts_m, aux))
aux_code = np.array([1,2,2])
codes_mp = np.hstack((codes_m, aux_code))
# Create path
marker_moment_pos = Path(verts_mp, codes_mp)

# Negative Moment
ang_p1 = ang_aux_neg - np.pi + ang_arrow
ang_p2 = ang_aux_neg - np.pi - ang_arrow
# Position of the tip of the arrow
tip_neg = verts_m[0,:]
# Position of the first line of the arrow
p1 = tip_neg + arrow_head_s*np.array([np.cos(ang_p1), np.sin(ang_p1)])
# Position of the second line of the arrow
p2 = tip_neg + arrow_head_s*np.array([np.cos(ang_p2), np.sin(ang_p2)])
# Add to the vertices and code
aux = np.array([p1, tip_neg, p2])
verts_mn = np.vstack((verts_m, aux))
aux_code = np.array([1,2,2])
codes_mn = np.hstack((codes_m, aux_code))
marker_moment_neg = Path(verts_mn, codes_mn)
############################################################

############################################################
# define markers for the supports
############################################################
# rolling x
x = [0,  1, -1, 0, -1.2,  1.2]
y = [0, -1, -1, 0, -1.5, -1.5]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.MOVETO, Path.LINETO]
roll_x = Path(xy, codes)
# rolling y
x = [0, -1, -1, 0, -1.5, -1.5]
y = [0,  1, -1, 0,  1.2, -1.2]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.MOVETO, Path.LINETO]
roll_y = Path(xy, codes)
# pinned
x = [0,  1, -1, 0, -1.2, 1.2, -1.2, -0.8, -0.8, -0.4, -0.4,  0,    0, 0.4,  0.4, 0.8,  0.8, 1.2]
y = [0, -1, -1, 0,   -1,  -1, -1.5,   -1, -1.5,   -1, -1.5, -1, -1.5,  -1, -1.5,  -1, -1.5,  -1]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.MOVETO, Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO ]
pinned = Path(xy, codes)
# encastrated
x = [-1, 1,  1, -1, -1, -1.2, -0.8, -0.8, -0.4, -0.4,  0,    0, 0.4,  0.4, 0.8,  0.8, 1.2]
y = [ 0, 0, -1, -1,  0, -1.5,   -1, -1.5,   -1, -1.5, -1, -1.5,  -1, -1.5,  -1, -1.5,  -1]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
        Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO ]
encas = Path(xy, codes)
# no rotation and no displacement in y
x = [-1, 1,  1, -1, -1, -1.2,  1.2]
y = [ 0, 0, -1, -1,  0, -1.5, -1.5]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO, Path.MOVETO, Path.LINETO]
disp_x = Path(xy, codes)
# no rotation and no displacement in x
x = [0,  0, -1, -1, 0, -1.5, -1.5]
y = [1, -1, -1,  1, 1,  1.2, -1.2]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
        Path.LINETO, Path.MOVETO, Path.LINETO]
disp_y = Path(xy, codes)
# only rotation constrained
x = [-0.5, 0.5, -0.5,  0.5]
y = [-0.5, 0.5,  0.5, -0.5]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
rot_z = Path(xy, codes)
############################################################

