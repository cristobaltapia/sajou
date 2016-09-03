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
        :**kwargs: TODO

        """
        self._model = model


class StaticSolver(Solver):

    """Linear Static solver used to solve the most typical problems"""

    def __init__(self, model, **kwargs):
        """Initialize the solver instance """
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
        print(K_red)
        # reduce the loading vector, by taking out the unknown loading
        # reactions
        P_red = np.delete(P, bc_ind, axis=0)
        # Solve the system
        V_res = np_solve(K_red, P_red)

        # Add the known displacements to the results obtained
        free_ind = np.arange(len(P))
        free_ind = np.delete(free_ind, bc_ind, axis=0)
        V[free_ind] = V_res

        P_react = np.dot(K, V)

        result = Result()
        result._V = V
        result._R = P_react

        return result
        
        
class Result(object):

    """An object to store the data obtained from the solving process"""

    def __init__(self):
        """TODO: to be defined1. """
        self._V = None # displacement results
        self._R = None # nodal reactions
        self._elem_force = None # nodal forces on the elements


