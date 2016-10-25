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
        
    def calc_moment_at(self, pos, element, unit_length=False):
        """Calculate the moment of element at position x.

        :pos: TODO
        :element: TODO
        :unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # Decide which data to use
        if unit_length == True:
            x_l = pos * element._length
        else:
            x_l = pos

        # Initialize moment
        try:
            moment = np.zeros(len(x_l))
        except:
            moment = 0.
        # Calculate for every load applied
        for load in element._loads:
            moment += self.calc_moment_with_member_load(element, load, x_l)

        num = element.number
        # Get the end forces results
        end_forces = self._result.data['end forces'] 
        # Add effect of end forces to the total moment at position 'pos'
        moment += end_forces[num][1]*x_l - end_forces[num][2]

        return moment

    def calc_shear_at(self, pos, element, unit_length=False):
        """Calculate the shear force of element at position x.

        :pos: TODO
        :element: TODO
        :unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # Decide which data to use
        if unit_length == True:
            x_l = pos * element._length
        else:
            x_l = pos

        # Initialize shear force results
        try:
            shear = np.zeros(len(x_l))
        except:
            shear = 0.
        # Calculate for every load applied
        for load in element._loads:
            shear += self.calc_shear_force_with_member_load(element, load, x_l)

        num = element.number
        # Get the end forces results
        end_forces = self._result.data['end forces'] 
        # Add effect of end forces to the total moment at position 'pos'
        shear += end_forces[num][1]

        return shear

    def calc_axial_at(self, pos, element, unit_length=False):
        """Calculate the axial force of element at position x.

        :pos: TODO
        :element: TODO
        :unit_length: boolean. Defines whether the range is [0, 1] or [0, Le] (default)
        :returns: TODO

        """
        # Decide which data to use
        if unit_length == True:
            x_l = pos * element._length
        else:
            x_l = pos

        # Initialize shear force results
        try:
            axial = np.zeros(len(x_l))
        except:
            axial = 0.
        # Calculate for every load applied
        for load in element._loads:
            axial += self.calc_axial_force_with_member_load(element, load, x_l)

        num = element.number
        # Get the end forces results
        end_forces = self._result.data['end forces'] 
        # Add effect of end forces to the total moment at position 'pos'
        axial += -end_forces[num][0]

        return axial

    def calc_axial_force_with_member_load(self, element, load, x):
        """Calculate the axial force in the given element.
        :returns: TODO

        :element:
        :load:
        :x:

        returns:
        """
        if load._direction == 'x' and load._coord_system == 'local':
            p1 = load._p1
            return -x * p1
        # TODO: case with global coord system
        else:
            axial = np.zeros(len(x))
            return axial

    def calc_shear_force_with_member_load(self, element, load, x):
        """Calculate the shear force in the given element.
        :returns: TODO

        :element:
        :load:
        :x:

        returns:
        """
        if load._direction == 'z' and load._coord_system == 'local':
            p1 = load._p1
            return x * p1
        # TODO: case with global coord system
        else:
            axial = np.zeros(len(x))
            return axial

    def calc_moment_with_member_load(self, element, load, x):
        """Calculate the shear force in the given element.
        :returns: TODO

        :element:
        :load:
        :x:

        returns:
        """
        # TODO: linearly varying distributed load
        if load._direction == 'z' and load._coord_system == 'local':
            p1 = load._p1
            return x**2 * 0.5 * p1
        # TODO: case with global coord system
        else:
            axial = np.zeros(len(x))
            return axial

    def calc_all_internal_forces(self, n=11):
        """Compute the internal forces at every element of the model.
        The number of points on which the internal forces are evaluated in each element
        are set by the variable 'n' (default n=11).

        :n: number of points at which the internal forces are evaluated. Must grater than
        2 (n>=2)
        :returns: TODO

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
            pos = np.linspace(0, 1, n) * curr_element._length
            # call a function to calculate the internal forces in a
            # single element
            moment = self.calc_internal_force_element(curr_element, pos, component='moment')
            shear = self.calc_internal_force_element(curr_element, pos, component='shear')
            axial = self.calc_internal_force_element(curr_element, pos, component='axial')

            # TODO: rest of the internal forces (shear and axial)
            res_moment = {'moment':{'data':moment, 'x':pos},
                          'shear': {'data':shear,  'x':pos},
                          'axial': {'data':axial,  'x':pos},
                          }
            #
            internal_forces[num_e] = res_moment
            # Calculate max and min
            max_internal_force[num_e] = {'moment':np.max(moment),
                                         'shear': np.max(shear),
                                         'axial': np.max(axial)
                                         }
            min_internal_force[num_e] = {'moemnt':np.min(moment),
                                         'shear': np.min(shear),
                                         'axial': np.min(axial)
                                         }

        # Get the maximum of the system
        max_df = pd.DataFrame(max_internal_force)
        min_df = pd.DataFrame(min_internal_force)
        df_aux = pd.concat([max_df, -min_df], axis=1)

        abs_max = df_aux.max(axis=1)

        max_force = max(max_internal_force)
        min_force = min(min_internal_force)
        # Create dictionary with the maximum and minimum data
        min_max_internal_forces = {'min': min_internal_force,
                                   'max': max_internal_force,
                                   'system abs max': abs_max}

        # Add results to Result object
        result = self._result
        result.add_result('internal forces', internal_forces)
        # Add metadata to the results object
        result.add_metadata('internal forces', min_max_internal_forces)

        return internal_forces

    def calc_internal_force_element(self, element, pos, component, unit_length=False):
        """Compute the internal forces of a given element.
        The variable 'pos' defines te position where the internal force is calculated. It
        can also be an array, specifying different positions at which the internal forces
        are needed.

        :element: a Beam instance
        :pos: position of evaluation (float or array)
        :component: TODO
        :unit_length: bool, sets wheather the local coordinate 'x' of the beam move between [0,1]
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
        """Compute the deflections of every element in the model.
        The number of points at whicn the deflections are calculated is specified in the
        variable 'n' (default n=11).

        :n: number of points at which the deflections are copmuted
        :returns: TODO

        """
        pass

