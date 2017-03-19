#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Postprocessing module.
"""
import numpy as np
import pandas as pd


class Postprocess(object):
    """Class for managing the postprocessing of results."""

    def __init__(self, result):
        """Initialize a Postprocess instance.

        :result: a Result instance

        """
        self._result = result
        self._model = result._model

    def calc_axial_at(self, pos, element, unit_length=False):
        """
        Calculate the axial force of element at position x.

        :param pos: TODO
        :param element: TODO
        :param unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # FIXME: verify first if 'pos' is array
        pos = np.array(pos)
        # Decide which data to use
        if unit_length:
            x_l = pos * element._length
        else:
            x_l = pos

        x_l = np.array([np.ones(pos.shape), x_l, x_l**2, x_l**3]).T

        # Calculate moment
        axial = x_l @ element._poly_sec_force[:, 0]

        return axial

    def calc_shear_at(self, pos, element, unit_length=False):
        """
        Calculate the shear force of element at position x.

        :param pos: TODO
        :param element: TODO
        :param unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # FIXME: verify first if 'pos' is array
        pos = np.array(pos)
        # Decide which data to use
        if unit_length:
            x_l = pos * element._length
        else:
            x_l = pos

        x_l = np.array([np.ones(pos.shape), x_l, x_l**2, x_l**3]).T

        # Calculate moment
        shear = x_l @ element._poly_sec_force[:, 1]

        return shear

    def calc_moment_at(self, pos, element, unit_length=False):
        """
        Calculate the moment of element at position x.

        :param pos: TODO
        :param element: TODO
        :param unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # FIXME: verify first if 'pos' is array
        pos = np.array(pos)
        # Decide which data to use
        if unit_length == True:
            x_l = pos * element._length
        else:
            x_l = pos

        x_l = np.array([np.ones(pos.shape), x_l, x_l**2, x_l**3]).T

        # Calculate moment
        moment = x_l @ element._poly_sec_force[:, 2]

        return moment

    def calc_all_internal_forces(self, n=11):
        """
        Compute the internal forces at every element of the model.

        The number of points on which the internal forces are evaluated in each element
        are set by the variable 'n' (default n=11).

        Parameters
        ----------

        n: int
            number of points at which the internal forces are evaluated. Must grater than
            2 (n>=2)

        Returns
        -------

        ndarray: array with internal forces

        """
        # Initialize dictionary for the results
        internal_forces = dict()
        # Check the number 'n'
        if n < 2:
            # FIXME: throw error
            return 0

        # Initialize list with maximum values of internal forces
        max_internal_force = dict()
        # The same with minimum values
        min_internal_force = dict()
        # A loop for each element of the model is made
        for num_e, curr_element in self._model.beams.items():
            # define the points
            pos = np.linspace(0, 1, n, dtype=np.float64) * curr_element._length
            # call a function to calculate the internal forces in a
            # single element
            moment = self.calc_internal_force_element(curr_element, pos,
                                                      component='moment')
            shear = self.calc_internal_force_element(curr_element, pos,
                                                     component='shear')
            axial = self.calc_internal_force_element(curr_element, pos,
                                                     component='axial')

            # TODO: rest of the internal forces (shear and axial)
            res_moment = {
                'moment': {
                    'data': moment,
                    'x': pos
                },
                'shear': {
                    'data': shear,
                    'x': pos
                },
                'axial': {
                    'data': axial,
                    'x': pos
                },
            }
            #
            internal_forces[num_e] = res_moment
            # Calculate max and min
            max_internal_force[num_e] = {
                'moment': np.max(moment),
                'shear': np.max(shear),
                'axial': np.max(axial)
            }
            min_internal_force[num_e] = {
                'moment': np.min(moment),
                'shear': np.min(shear),
                'axial': np.min(axial)
            }

        # Get the maximum of the system
        max_df = pd.DataFrame(max_internal_force)
        min_df = pd.DataFrame(min_internal_force)
        df_aux = pd.concat([max_df, np.abs(min_df)], axis=1)

        abs_max = df_aux.max(axis=1)

        max_force = max(max_internal_force)
        min_force = min(min_internal_force)
        # Create dictionary with the maximum and minimum data
        min_max_internal_forces = {
            'min': min_internal_force,
            'max': max_internal_force,
            'system abs max': abs_max
        }

        # Add results to Result object
        result = self._result
        result.add_result('internal forces', internal_forces)
        # Add metadata to the results object
        result.add_metadata('internal forces', min_max_internal_forces)

        return internal_forces

    def calc_internal_force_element(self, element, pos, component,
                                    unit_length=False):
        """
        Compute the internal forces of a given element.
        The variable 'pos' defines te position where the internal force is calculated. It
        can also be an array, specifying different positions at which the internal forces
        are needed.

        :param element: a Beam instance
        :param pos: position of evaluation (float or array)
        :param component: TODO
        :param unit_length: bool, sets wheather the local coordinate 'x' of the beam move between [0,1]
               or [0,Le]

        :returns: TODO

        """
        if component == 'moment':
            internal_force = self.calc_moment_at(pos, element, unit_length)
        elif component == 'shear':
            internal_force = self.calc_shear_at(pos, element, unit_length)
        elif component == 'axial':
            internal_force = self.calc_axial_at(pos, element, unit_length)

        return internal_force

    def calc_all_deflections(self, n=11):
        """
        Compute the deflections of every element in the model.
        The number of points at whicn the deflections are calculated is specified in the
        variable 'n' (default n=11).

        :param n: number of points at which the deflections are copmuted
        :returns: TODO

        """
        pass

    def calc_node_new_coord(self):
        """Calculate the new coordinates of the nodes in the deformed state
        :returns: TODO

        """
        # get the nodes
        nodes_coords = self._result['nodal displacements']
        # initilize array for the numbers of the nodes in the element
        arr_nodes = np.zeros(len(nodes_coords))
        for key, node_i in nodes_coords.items():
            arr_nodes[key] = node_i.number

        # Get original coordinates of nodes
        coords_nodes = get_dataframe_of_node_coords(model=result._model,
                                                    nodes=arr_nodes)
        # Get displacements of the nodes
        displ_nodes = result.data['nodal displacements'].loc[arr_nodes]

        return coords_nodes, displ_nodes
