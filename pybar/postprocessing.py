#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Postprocessing module.
"""

class Postprocess(object):

    """Class for managing the postprocessing of results."""

    def __init__(self, result):
        """Initialize a Postprocess instance.

        :result: a Result instance

        """
        self._result = result
        self._model = result._model
        
    def _calc_internal_forces(self):
        """Calculate member forces

        :result: TODO
        :returns: TODO

        """
        # First calculate end forces on each member
        end_forces = self._result.end_forces

        max_N = 0.
        max_V = 0.
        max_M = 0.

        for num, elem in self._model.beams.items():
            # Points where data is going to be extracted
            x_axis = np.linspace(0,1,11) * elem._length
            axial = np.zeros(len(x_axis))
            shear = np.zeros(len(x_axis))
            moment = np.zeros(len(x_axis))
            # Take account for element loads (distributed, etc)
            # FIXME: still need to include distributed loads, etc
            if len(elem._loads) > 0:
                for load in elem._loads:
                    axial += self.calc_axial_force_with_member_load(elem, load, x_axis)
                    shear += self.calc_shear_force_with_member_load(elem, load, x_axis)
                    moment += self.calc_moment_with_member_load(elem, load, x_axis)

            # Calculate axial force on the member
            axial -= np.ones(len(x_axis)) * end_forces[num][0]
            # Calculate shear force on the member
            shear += np.ones(len(x_axis)) * end_forces[num][1]
            # Calculate moment on the member
            moment +=  end_forces[num][1]*x_axis - end_forces[num][2]
            # Add to the results
            res = {'axial':axial, 'shear':shear, 'moment':moment, 'x':x_axis}
            result.internal_forces[num] = res
            # maximum absolute values
            max_N = max([max(abs(axial)), max_N])
            max_V = max([max(abs(shear)), max_V])
            max_M = max([max(abs(moment)), max_M])

        result._max_member_force['axial'] = max_N
        result._max_member_force['shear'] = max_V
        result._max_member_force['moment'] = max_M

        return result.internal_forces

    def calc_moment_at(self, x, element, unit_length=False):
        """Calculate the moment of element at position x.

        :x: TODO
        :element: TODO
        :unit_length: boolean. Defines whether the range is [0, L] (default) or [0, 1]
        :returns: TODO

        """
        moment = 0.
        for load in element._loads:
            moment += self.calc_moment_with_member_load(element, load, x)

        num = element.number
        moment +=  self._result.end_forces[num][1]*x - self._result.end_forces[num][2]

        return moment

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

