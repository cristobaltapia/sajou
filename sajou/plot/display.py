#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Define all the functions needed to plot the created structures and the results

TODO: plot hinges
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn.apionly as sns
from matplotlib.axes import Axes
from matplotlib.patches import Polygon
from matplotlib.path import Path
from sajou.model import get_dataframe_of_node_coords, get_node_coords

from .lines_mpl import Line2D
from .markers_mpl import MarkerStyle
from .plot_elements_mpl import plot_element


class Display(object):
    """
    Class to present the pre and postprocessing graphically.

    This class is not intended to be used directly, but rather their sublclasses, which
    implement different plotting backends (at the moent only matplotlib).

    Parameters
    ----------

    width: float
        width of the window in pixels
    height: float
        height of the window in pixels
    theme: str
        color theme used for the plots. Options are: 'dark', 'light'
        and 'publication'
    backend: str
        backend used to display the results

        - 'matplotlib'
        - 'python-occ' (not implemented yet)
        - 'mayavi' (not implemented yet)

    """

    def __new__(cls, backend='matplotlib', **kwargs):
        if cls is Display:
            if backend == 'matplotlib':
                return super(Display, cls).__new__(Display_mpl)
        else:
            return super(Display, cls).__new__(cls, backend)

    def __init__(self, backend):
        """Instatiate a Display object. """
        self._backend = backend


class Display_mpl(Display):
    """Matplotlib backend for the visualization

    Parameters
    ----------

    width: float
        width of the window in pixels
    height: float
        height of the window in pixels
    theme: str
        color theme used for the plots. Options are: 'dark', 'light'
        and 'publication'

    """

    def __init__(self, width=10., height=10., theme='dark',
                 backend='matplotlib'):
        Display.__init__(self, backend)

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
                'force': {
                    'color': palette[1],
                    'markeredgewidth': 2,
                },
                'reaction': {
                    'color': 'yellow',
                    'markeredgewidth': 2,
                },
                'background': {
                    'color': "#12191f"
                },
                'grid': {
                    'color': 'white',
                    'ls': ':',
                },
                'element': {
                    'color': palette[4],
                    'linewidth': 1
                },
                'support': {
                    'markeredgecolor': palette[5],
                    'markerfacecolor': 'None',
                },
                'node': {
                    'color': palette[4]
                },
                'member force positive': {
                    'edgecolor': 'white',
                    'facecolor': palette[0],
                    'alpha': 0.5
                },
                'member force negative': {
                    'edgecolor': 'white',
                    'facecolor': palette[2],
                    'alpha': 0.5
                },
            }
        elif theme == 'publication':
            self.draw_config = {
                'force': {
                    'color': 'black',
                    'markeredgewidth': 2,
                },
                'reaction': {
                    'color': 'black',
                    'markeredgewidth': 2,
                },
                'background': {
                    'color': 'white'
                },
                'grid': {
                    'color': 'black',
                    'ls': ':',
                },
                'element': {
                    'color': 'black',
                    'linewidth': 2
                },
                'support': {
                    'markeredgecolor': 'black',
                    'markerfacecolor': 'None',
                    'markeredgewidth': 2,
                },
                'node': {
                    'color': 'black'
                },
                'member force positive': {
                    'edgecolor': 'black',
                    'facecolor': 'None',
                    'alpha': 0.5
                },
                'member force negative': {
                    'edgecolor': 'black',
                    'facecolor': 'None',
                    'alpha': 0.5
                },
            }
        elif theme == 'light':
            palette = sns.color_palette('dark')
            self.draw_config = {
                'force': {
                    'color': palette[1],
                    'markeredgewidth': 2,
                },
                'reaction': {
                    'color': palette[3],
                    'markeredgewidth': 2,
                },
                'background': {
                    'color': "#dee5ec"
                },
                'grid': {
                    'color': 'black',
                    'ls': ':',
                },
                'element': {
                    'color': palette[4],
                    'linewidth': 1
                },
                'support': {
                    'markeredgecolor': palette[5],
                    'markerfacecolor': 'None',
                },
                'node': {
                    'color': palette[4]
                },
                'member force positive': {
                    'edgecolor': 'black',
                    'facecolor': palette[0],
                    'alpha': 0.5
                },
                'member force negative': {
                    'edgecolor': 'black',
                    'facecolor': palette[2],
                    'alpha': 0.5
                },
            }
        else:
            palette = sns.color_palette('dark')
            self.draw_config = {
                'force': {
                    'color': palette[1]
                },
                'reaction': {
                    'color': palette[3]
                },
                'background': {
                    'color': "#dee5ec"
                },
                'grid': {
                    'color': 'black',
                    'ls': ':',
                },
                'element': {
                    'color': palette[4],
                    'linewidth': 1
                },
                'support': {
                    'markeredgecolor': palette[5],
                    'markerfacecolor': 'None',
                    'ms': 25
                },
                'node': {
                    'color': palette[4]
                },
                'member force positive': {
                    'edgecolor': 'black',
                    'facecolor': palette[0],
                    'alpha': 0.5
                },
                'member force negative': {
                    'edgecolor': 'black',
                    'facecolor': palette[2],
                    'alpha': 0.5
                },
            }

    def plot_geometry(self, model, ax, ls='-', **kwargs):
        """Plot the geometry of the model passed.

        Parameters
        ----------

        model: sajou.Model
            Model to plot
        ax: a matplotlib axis object
            Axis where to draw to

        Returns
        -------
        matplotlib axis:
            the same axis object given

        """
        node_options = self.draw_config['node']
        elem_options = self.draw_config['element']
        background_options = self.draw_config['background']
        grid_options = self.draw_config['grid']
        force_options = self.draw_config['force']
        support_options = self.draw_config['support']
        show_loads = kwargs.get('show_loads', self.display_config['forces'])
        # add line style to the 'elem_options' dictionary
        elem_options['ls'] = ls

        # set background color
        ax.set_facecolor(background_options['color'])
        ax.grid(True, **grid_options)

        # call the plot function for each element
        for num, elem in model.elements.items():
            plot_element(elem, ax, elem_options)

        nodes = model.nodes
        if self.display_config['nodes']:
            for num, node in nodes.items():
                m_nodes = Line2D([node.x], [node.y], marker='o', ls='',
                                 fillstyle='full', **node_options)
                ax.add_line(m_nodes)

        # Plot forces if requiered
        if show_loads:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i._loads.items():
                    ax = self.plot_nodal_force(ax, dof, at=node_i, val=val)
            # Plot element Loads
            ax = self.plot_element_loads(ax, model)

        # Plot supports
        if self.display_config['supports']:
            for ix, node_i in model.nodes.items():
                if len(node_i._bc) > 0:
                    ax = self.plot_support(ax, dof=node_i._bc.keys(),
                                           at=node_i)

        ax.axis('equal')

        return ax

    def plot_deformed_geometry(self, result, ax, show_undeformed=False,
                               scale=1.):
        """Plot the system in its deformed configuration.

        Parameters
        ----------

        result: Result object
            result object returned by the Solver
        ax: matplotlib axis
            Axis where to draw to
        scale: float
            scale used to plot the deformations

        Return
        ------

        matplotlib axis:
            the same axis object passed

        """
        # get the style configurations
        node_options = self.draw_config['node']
        elem_options = self.draw_config['element']
        background_options = self.draw_config['background']
        grid_options = self.draw_config['grid']
        force_options = self.draw_config['force']
        support_options = self.draw_config['support']
        #
        model = result._model

        # set background color
        ax.set_facecolor(background_options['color'])
        ax.grid(True, **grid_options)

        # Draw each deformed element
        for num, elem in model.elements.items():
            # get node coordinates and displacements
            displ_nodes = get_element_deformed_node_coords(elem, result, scale)
            # convert to ndarray
            node_def = np.array([(d[0], d[1]) for k, d in displ_nodes.items()])
            # Plot deformed element
            ax.plot(node_def[:, 0], node_def[:, 1], **elem_options)

        nodes = model.nodes
        if self.display_config['nodes']:
            # Get original coordinates of nodes
            coords_nodes = get_deformed_node_coords(nodes='all', result=result,
                                                    scale=scale)
            # convert to ndarray
            node_def = np.array(
                [(d[0], d[1]) for k, d in coords_nodes.items()])
            # plot nodes
            ax.plot(node_def[:, 0], node_def[:, 1], ls='', marker='o', **node_options)

        if show_undeformed:
            ax = self.plot_geometry(model=result._model, ax=ax, ls='--')

        ax.axis('equal')

        return ax

    def plot_deformed_element(self, element, result, ax):
        """Plot the given element in its deformed state

        Parameters
        ----------

        element: Element instance
            The element to be plotted
        result: Result instance
            The result obtained after solving the system
        ax: matplotlib Axis
            the axis where to plot the deformed element to

        Returns
        -------

        matplotlib Axis:
            the axis used to plot the element

        """
        return ax

    def plot_element_loads(self, ax, model):
        """Plot the element loads.

        Element loads like uniformly distributed loads are drawn accordinlgy.

        Parameters
        ----------

        ax: matplotlib axis
            axis where to draw to
        result: Result object
            Result object returned by the Solver

        Returns
        -------

        matplotlib axis:
            the same axis originally passed

        """
        elems_with_loads = {
            ix: elem
            for ix, elem in model.beams.items() if elem._loads != []
        }

        # Get maximum distributed load
        if len(elems_with_loads) > 0:
            pe = []
            for ix, elem in elems_with_loads.items():
                for load in elem._loads:
                    if load._type == 'Distributed Load':
                        try:
                            p_max = np.max(
                                [np.abs(load._p1), np.abs(load._p2)])
                            pe.append(p_max)
                        except:
                            pass

            # If there are only distributed moments this will fail, that
            # is why a 'try' is used (FIXME maybe)
            try:
                pmax = np.max(pe)
            except:
                pmax = 1.

        # FIXME: make dependent of the 'bounding box'
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
                    ax = self._plot_distributed_load_element(
                        ax, model, elem, load_i, scale=scale, max_p=pmax)
                elif load_i._type == 'Distributed Moment':
                    ax = self._plot_distributed_moment_element(
                        ax, model, elem, load_i, scale=scale, max_p=pmax)

        return ax

    def plot_internal_forces(self, ax, result, component, scale=1):
        """Plot the diagrams of internal forces.

        Parameters
        ----------

        ax: matplotlib axis
            the Axis instance where to draw the internal forces
        result: Result object
            obtained after solving the system
        component: str
            component that want to be displayed. Options are

                - 'axial'
                - 'shear'
                - 'moment'

        scale: float
            scale for the member forces (the member forces are already
            automatically scaled to fit in the display. This is a scale factor
            that multiplies that automatically calulated factor)

        Returns
        -------

        matplotlib axis:
            the axis with the drawn internal forces

        """
        member_force_options_pos = self.draw_config['member force positive']
        member_force_options_neg = self.draw_config['member force negative']
        model = result._model

        # First plot the geometry
        ax = self.plot_geometry(model, ax, show_loads=False)

        # FIXME: instead of using an auxiliary figure to create the
        # polygon and then plot, calculate the respective points of the
        # polygon with scipy.interp()
        # auxiliary axes
        fig_aux = plt.figure()
        ax_aux = fig_aux.add_subplot(111)
        # Determine automatic scaling factor
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        max_range = max([x_range, y_range])
        max_component = result.metadata['internal forces']['system abs max'][
            component]
        scale_auto = 0.2 * max_range / max_component
        # Plot the specified diagram
        # Transform the results from the local coordinates to global
        # coordinates
        for num, elem in model.elements.items():
            # Get the transformation matrix
            T = elem.transformation_matrix[0:2, 0:2]
            # Get the positions at which the internal forces were
            # calculated
            x_axis = result.data['internal forces'][num][component]['x']
            # Get the specified member forces
            d = result.data['internal forces'][num][component][
                'data'] * scale_auto * scale

            # Results are multiplied by -1 to represent them in the convention
            # where positive moments are drawn 'below' the beam element.
            d = -d
            # Positive values
            d_pos = ax_aux.fill_between(x_axis, d, 0, where=d < 0,
                                        interpolate=True)
            d_pos = d_pos.get_paths()
            # Negative values
            d_neg = ax_aux.fill_between(x_axis, d, 0, where=d > 0,
                                        interpolate=True)
            d_neg = d_neg.get_paths()

            # Plot the patches
            # positive part
            for curr_polygon in d_pos:
                rot_pos = np.dot(curr_polygon.vertices, T.todense())
                rot_pos[:, 0] += elem._node1.x
                rot_pos[:, 1] += elem._node1.y
                poly_pos = Polygon(rot_pos, True, **member_force_options_pos)
                ax.add_patch(poly_pos)

            # negative part
            for curr_polygon in d_neg:
                rot_neg = np.dot(curr_polygon.vertices, T.todense())
                rot_neg[:, 0] += elem._node1.x
                rot_neg[:, 1] += elem._node1.y
                poly_neg = Polygon(rot_neg, True, **member_force_options_neg)
                ax.add_patch(poly_neg)

        # Plot reaction forces
        if self.display_config['reactions']:
            for ix, node_i in model.nodes.items():
                for dof, val in node_i.reactions.items():
                    ax = self.plot_nodal_reaction(ax, dof, at=node_i, val=val)

        # close the auxiliary figure
        plt.close(fig_aux)

        return ax

    def plot_support(self, ax, dof, at):
        """Plot the respective support

        Parameters
        ----------

        ax: matplotlib Axis
            axis where to plot the support
        dof: ndarray, list[int]
            DOF that are restrained
        at: Node object
            node at which the contraints are applied

        Returns
        -------

        matplotlib Axis
            the same axis passed

        """
        support_options = self.draw_config['support']
        # transform to list
        dof = list(dof)
        #

        if len(dof) == 1:
            # rolling support free in 'y'
            if dof[0] == 0:
                suppt = Line2D([at.x], [at.y], marker='rsy', **support_options)
            # rolling support free in 'x'
            elif dof[0] == 1:
                suppt = Line2D([at.x], [at.y], marker='rsx', **support_options)
            # only rotation constrained
            elif dof[0] == 2:
                suppt = Line2D([at.x], [at.y], marker=rot_z, **support_options)
        #
        if len(dof) == 2:
            # pinned support
            if np.all([0, 1] == np.sort(dof)):
                suppt = Line2D([at.x], [at.y], marker='psx', **support_options)
            #
            elif np.all([0, 2] == np.sort(dof)):
                suppt = Line2D([at.x], [at.y], marker='rey', **support_options)
            #
            elif np.all([1, 2] == np.sort(dof)):
                suppt = Line2D([at.x], [at.y], marker='rex', **support_options)
        #
        if len(dof) == 3:
            # Encastrated
            suppt = Line2D([at.x], [at.y], marker='es', **support_options)

        ax.add_line(suppt)

        return ax

    def plot_nodal_force(self, ax, dof, at, val):
        """TODO: Docstring for plot_forces.

        Parameters
        ----------

        ax: matplotlib Axis
            axis to plot to
        dof: int
            Degree of freedom to which the force is being applied
        at: Node object
            node at which the load is applied
        val: float
            value of the force

        Returns
        -------

        matplotlib Axis
            the same axis passed as asrgument

        """
        # Get the draw options for the forces
        force_options = self.draw_config['force']

        # Force in x direction
        if dof == 0:
            if val < 0:
                halign = 'left'
            else:
                halign = 'right'

            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(-np.sign(val) * 50, 0),
                        color=force_options['color'], ha=halign, va='center',
                        textcoords='offset points', arrowprops=dict(
                            arrowstyle='->',
                            linewidth=force_options['markeredgewidth'],
                            color=force_options['color']))
        # Force in y direction
        elif dof == 1:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(0, -np.sign(val) * 50),
                        color=force_options['color'], ha='center', va='center',
                        textcoords='offset points', arrowprops=dict(
                            arrowstyle='->',
                            linewidth=force_options['markeredgewidth'],
                            color=force_options['color']))
        # Moment
        elif dof == 2:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(20, 20), color=force_options['color'],
                        ha='left', va='center', textcoords='offset points')
            if np.sign(val) >= 0.:
                mom = Line2D([at.x], [at.y], marker='ap',
                             markerfacecolor='None', **force_options)
                ax.add_line(mom)
            else:
                mom = Line2D([at.x], [at.y], marker='an',
                             markerfacecolor='None', **force_options)
                ax.add_line(mom)

        return ax

    def plot_nodal_reaction(self, ax, dof, at, val):
        """Plot the nodal reactions of the system

        Parameters
        ----------

        ax: matplotlib Axis
            axis to plot to
        dof: int
            Degree of freedom to which the force is being applied
        at: Node object
            node at which the load is applied
        val: float
            value of the reaction

        Returns
        -------

        matplotlib Axis
            the same Axis passed as argument

        """
        # Get the draw options for the forces
        force_options = self.draw_config['reaction']

        # Force in x direction
        if dof == 0:
            if val < 0:
                halign = 'left'
            else:
                halign = 'right'

            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(-np.sign(val) * 50, 0),
                        color=force_options['color'], ha=halign, va='center',
                        textcoords='offset points', arrowprops=dict(
                            arrowstyle='->', color=force_options['color'],
                            lw=force_options['markeredgewidth']))
        # Force in y direction
        elif dof == 1:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(0, -np.sign(val) * 50),
                        color=force_options['color'], ha='center', va='center',
                        textcoords='offset points', arrowprops=dict(
                            arrowstyle='->', color=force_options['color'],
                            lw=force_options['markeredgewidth']))
        # Moment
        elif dof == 2:
            ax.annotate('{f:.2E}'.format(f=abs(val)), xy=(at.x, at.y),
                        xytext=(20, 20), color=force_options['color'],
                        ha='left', va='center', textcoords='offset points')
            if np.sign(val) >= 0.:
                mom = Line2D([at.x], [at.y], marker='ap',
                             markeredgecolor=force_options['color'],
                             markerfacecolor='none', **force_options)
                ax.add_line(mom)
            else:
                mom = Line2D([at.x], [at.y], marker='an',
                             markeredgecolor=force_options['color'],
                             markerfacecolor='none', **force_options)
                ax.add_line(mom)

        return ax

    def _plot_distributed_load_element(self, ax, model, element, load, scale=1,
                                       max_p=1):
        """Plot a distributed load at the corresponding element.

        Parameters
        ----------

        ax: matplotlib axis
            axis where to draw to
        model: Model object
            Model containing the element
        element: Element object
            element for which the distributed loads will be drawn
        load: Load object
            The Load object that wants to be drawn

        Returns
        -------

        matplotlib axis:
            the same axis passed

        """
        # Get the draw options for the forces
        force_options = self.draw_config['force']
        #
        size = scale * 0.08
        #
        n1 = element._node1
        n2 = element._node2
        p1 = load._p1
        p2 = load._p2
        # Number of arrows
        n_arrows = (element._length * 2) // size
        # Generate points on the element line
        nx_e = np.linspace(n1.x, n2.x, n_arrows)
        ny_e = np.linspace(n1.y, n2.y, n_arrows)

        if load._direction == 'y':
            # Generate points on the element line
            nx_e = np.linspace(n1.x, n2.x, n_arrows)
            ny_e = np.linspace(n1.y, n2.y, n_arrows)
            # Offset from the position of the nodes...
            T = element.transformation_matrix
            # Differentiate between local and global coordinates
            # - For local coordinates
            if load._coord_system == 'local':
                offset = T.T.dot(
                    np.array([0, -p1, 0, 0, -p2, 0]) * scale / max_p * 0.2)
            # - For global coordinates
            elif load._coord_system == 'global':
                offset = np.array([0, -p1, 0, 0, -p2, 0]) * scale / max_p * 0.2

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
                ax.annotate('', xy=(nx_e[arr_i], ny_e[arr_i]),
                            xytext=(nx[arr_i], ny[arr_i]), arrowprops=dict(
                                arrowstyle='->', color=force_options['color']))

        elif load._direction == 'x':
            # Offset from the position of the nodes...
            # FIXME: ... according to the system coordinate chosen
            T = element.transformation_matrix
            # Differentiate between local and global coordinates
            # - For local coordinates
            if load._coord_system == 'local':
                # Generate points on the element line
                nx_e = range_with_ratio(n1.x, n2.x, n_arrows, p1 / p2)
                ny_e = range_with_ratio(n1.y, n2.y, n_arrows, p1 / p2)
                offset = T.T.dot(
                    np.array([-p1, 0, 0, -p2, 0, 0]) * size / max_p * 0.8)
            # - For global coordinates
            elif load._coord_system == 'global':
                # Generate points on the element line
                nx_e = np.linspace(n1.x, n2.x, n_arrows)
                ny_e = np.linspace(n1.y, n2.y, n_arrows)
                offset = np.array([-p1, 0, 0, -p2, 0, 0]) * scale / max_p * 0.2

            n1_o = np.array(n1[:]) + offset[:3]
            n2_o = np.array(n2[:]) + offset[3:]
            # Generate points for the arrows with offset
            nx = range_with_ratio(n1_o[0], n2_o[0], n_arrows, p1 / p2)
            ny = range_with_ratio(n1_o[1], n2_o[1], n_arrows, p1 / p2)
            if load._coord_system == 'global':
                ax.plot(nx, ny, color=force_options['color'])

            # Plot arrows
            for arr_i in range(len(nx)):
                # annotate() is used instead of arrow() because the style
                # of the arrows is better
                if load._coord_system == 'local':
                    ax.annotate('', xy=(nx_e[arr_i], ny_e[arr_i]),
                                xytext=(nx[arr_i], ny[arr_i]),
                                arrowprops=dict(arrowstyle='->',
                                                color=force_options['color']))
                elif load._coord_system == 'global':
                    ax.annotate('', xy=(nx_e[arr_i], ny_e[arr_i]),
                                xytext=(nx[arr_i], ny[arr_i]),
                                arrowprops=dict(arrowstyle='->',
                                                color=force_options['color']))

        return ax

    def _plot_distributed_moment_element(self, ax, model, element, load,
                                         scale=1, max_p=1):
        """Plot a distributed load at the corresponding element.

        Parameters
        ----------

        ax: matplotlib axis
            axis where to draw to
        model: Model object
            Model used
        element: Element object
            Element for which the distributed moment is drawn
        load: Load object
            Load that wants to be drawn

        Returns
        -------

        matplotlib axis:
            the same axis passed

        """
        # Get the draw options for the forces
        force_options = self.draw_config['force']
        #
        size = scale * 0.08
        #
        n1 = element._node1
        n2 = element._node2
        m1 = load._m1
        m2 = load._m2
        # Number of arrows
        n_arrows = (element._length * 2) // size
        # Generate points on the element line
        nx_e = np.linspace(n1.x, n2.x, n_arrows)
        ny_e = np.linspace(n1.y, n2.y, n_arrows)

        if load._direction == 'z' and load._coord_system == 'local':
            # Offset from the position of the nodes...
            # FIXME: ... according to the system coordinate chosen
            T = element.transformation_matrix
            # Generate points for the arrows with offset
            nx_e = range_with_ratio(n1.x, n2.x, n_arrows, m1 / m2)
            ny_e = range_with_ratio(n1.y, n2.y, n_arrows, m1 / m2)

            # Plot arrows
            ax.plot(
                nx_e,
                ny_e,
                marker='ap',
                ls='None',
                ms=20,
                markeredgewidth=1.5,
                markeredgecolor=force_options['color'],
                markerfacecolor='None', )

        return ax


def get_element_deformed_node_coords(element, result, scale=1):
    """Get the deformed coordinates of the nodes of the given element

    Parameters
    ----------

    element: ELement instance
        Element from which the deformed position of its nodes will be computed
    result: Result object
        results obtained from the solver
    sacele: float
        scale applied to the deformations (default to 1)

    Returns
    -------

    dict:
        node deformed coordinates

    """
    # get the nodes
    nodes_e = element._nodal_connectivity
    # List of the nodes in the element
    list_nodes = [node for k, node in nodes_e.items()]
    # Get original coordinates of nodes
    coord_nodes = get_node_coords(model=result._model, nodes=list_nodes)
    # Get displacements of the nodes
    node_displ = result.data['nodal displacements']

    # FIXME: make this compatible for 3D
    displ_nodes = dict()
    for n in list_nodes:
        displ_nodes[n.number] = coord_nodes[n.number] + node_displ[
            n.number][:2] * scale

    return displ_nodes


def get_deformed_node_coords(nodes, result, scale=1):
    """Get the deformed coordinates of the given nodes

    Parameters
    ----------

    nodes: List[Node], str
        Node from which the deformed position will be computed
    result: Result object
        results obtained from the solver
    sacele: float
        scale applied to the deformations (default to 1)

    Returns
    -------

    dict:
        node deformed coordinates

    """
    model = result._model
    if nodes == 'all':
        list_nodes = [node for k, node in model.nodes.items()]
    else:
        list_nodes = nodes

    # Get original coordinates of nodes
    coord_nodes = get_node_coords(model=result._model, nodes=list_nodes)
    # Get displacements of the nodes
    node_displ = result.data['nodal displacements']

    # FIXME: make this compatible for 3D
    displ_nodes = dict()
    for n in list_nodes:
        displ_nodes[n.number] = coord_nodes[n.number] + node_displ[
            n.number][:2] * scale

    return displ_nodes


def range_with_ratio(x1, x2, n, a):
    """Create an array with numbers increasing their relative distance
    according to the ratio 'a'.

    The range is divided such that:

    .. math::

        \ell_{tot} = \\ell_1 + \\ell_2 + ... + \\ell_n

    where:

    .. math::

        \\ell_i &= \\ell_i - 1 + b \\\\
        b &= \\frac{\ell_{tot} - (\\frac{\ell_{tot}}{n} \cdot a)}{n}

    Parameters
    ----------

    x1: float
        begin of the range
    x2: float
        end of the range
    n: int
        number of points
    a: float
        ratio between the beginning and end

    Returns
    -------

    ndarray
        a range with increasingly separated values

    """
    # base array
    base = np.ones(n)
    # add to obtain required ratio
    base = base * np.linspace(0, a, n)
    # Put in the desired range
    new_range = (base / np.max(base) * (x2 - x1)) + x1

    return new_range


# only rotation constrained
x = [-0.5, 0.5, -0.5, 0.5]
y = [-0.5, 0.5, 0.5, -0.5]
xy = list(zip(x, y))
codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
rot_z = Path(xy, codes)
############################################################
