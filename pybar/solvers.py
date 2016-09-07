#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import solve as np_solve
#import scipy as sp
""" This module contains the different solvers used by the program.
"""

class Solver(object):

    """Parent class for the solvers used"""

    def __init__(self, model, **kwargs):
        """Initialize solver

        :model: TODO
        :**kwargs: Additional options:

            - :output: list of keywords specifying the different results that need to be
              computed after the solution is found.
              Accepted values are:
                - 'reactions'
                - 'internal forces'
                - 'stresses'
                - 'strains'
                - 'energy'


        """
        self._model = model
        # The variable 'output' stores the different results that need
        # to be computed after the solution is found.
        self._output_request = kwargs.get('output', None)


class StaticSolver(Solver):

    """Linear Static solver used to solve the most typical problems"""

    def __init__(self, model, **kwargs):
        """Initialize the solver instance

        Optional parameters:
        nlg: boolean, consider geometrical nonlinearities (Default False)

        """
        Solver.__init__(self, model, **kwargs)

    def solve(self):
        """Solves the problem
        :returns: numpy array

        """
        # Assemble the stiffness matrix of the system
        K = self._model._assemble_global_K()
        # Assemble the vector of applied loads
        P = self._model._generate_loading_vector()
        # Apply border conditions
        V = self._model._generate_displacement_vector()
        #FIXME: this can be done when establishing the incidence table
        # dof corresponding to border conditions
        bc_ind = self._model._dof_restraied
        K_red = np.delete(K, bc_ind, axis=0)
        K_red = np.delete(K_red, bc_ind, axis=1)

        # reduce the loading vector, by taking out the unknown loading
        # reactions
        P_red = np.delete(P, bc_ind, axis=0)
        # Solve the system
        V_res = np_solve(K_red, P_red)

        # Add the known displacements to the results obtained
        free_ind = np.arange(len(P))
        free_ind = np.delete(free_ind, bc_ind, axis=0)
        V[free_ind] = V_res

        P_react = np.dot(K, V)[bc_ind]

        result = Result()
        result._V = V
        result._R = P_react

        self.postprocess(result)

        return result

    def postprocess(self, result):
        """Calculates the specified results given in the 'output' variable

        :result: TODO
        :returns: TODO

        """
        for curr_output in self._output_request:
            if curr_output == 'internal forces':
                self._calc_internal_forces(result)
            elif curr_output == 'stresses':
                self._calc_stresses(result)

        return result

    def _calc_internal_forces(self, result):
        """Calculate member forces

        :result: TODO
        :returns: TODO

        """
        # First calculate end forces on each member
        end_forces = self._calc_end_forces(result)

        max_N = 0.
        max_V = 0.
        max_M = 0.

        # FIXME: still need to include distributed loads, etc
        for num, elem in self._model.beams.items():
            x_axis = np.linspace(0,1,11) * elem._length
            # Calculate axial force on the member (a function is returned)
            axial = np.ones(len(x_axis)) * -end_forces[num][0]
            # Calculate shear force on the member
            shear = np.ones(len(x_axis)) * -end_forces[num][1]
            # Calculate moment on the member
            moment = -shear*x_axis - end_forces[num][2]
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

    def _calc_end_forces(self, result):
        """Calculates the internal forces of beam elements

        :result: TODO
        :returns: TODO

        """
        # calculate for each element
        for num, elem in self._model.beams.items():
            # Get the transformation matrix for the element
            T = elem.transformation_matrix
            # Get the stiffness matrix of the element
            Ke = elem._Ke
            # Get the displacements of the corresponding DOFs
            dof_pn = elem._dof_per_node
            v_i = np.zeros(dof_pn*elem._n_nodes)
            for n_node, node in enumerate(elem._nodes):
                i1 = n_node*dof_pn
                i2 = dof_pn*(n_node+1)
                j1 = node.number*dof_pn
                j2 = dof_pn*(1+node.number)
                v_i[i1:i2] = result._V[j1:j2]

            P_i = np.dot(Ke, v_i)

            P_i_local = np.dot(T, P_i)

            result.end_forces[num] = P_i_local

        return result.end_forces


class Result(object):

    """An object to store the data obtained from the solving process"""

    def __init__(self):
        """TODO: to be defined1. """
        self._V = None # displacement results
        self._R = None # nodal reactions
        self.end_forces = dict() # end forces of elements
        self.internal_forces = dict() # internal forces of elements
        # maximum absolut value of the different member forces in the
        # model
        self._max_member_force = dict()



