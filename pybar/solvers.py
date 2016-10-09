#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.linalg import lu_factor, lu_solve
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
        # FIXME: implement this
        deffault_output = ['displacements', 'forces']
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
        # Assemble the vector of applied element loads
        P_e = self._model._generate_element_loading_vector()
        # Apply Dirichlet boundary conditions
        # Get the dof corresponding to Dirichlet border conditions
        dirich_ind = self._model._dof_dirichlet
        V = self._model._generate_displacement_vector()
        # Calculate forces generated by the BC
        P_dirich = np.dot(K[dirich_ind,:], V)
        P_new = P[:]
        P_new[dirich_ind] = P[dirich_ind] - P_dirich
        # Add the element loads to the global loading vector
        P_new += P_e

        # Matrix to multiply the global stiffness matrix and put zeros on the
        # rows and columns corresponding to the Dirichlet BC.
        #T = np.eye(len(K))
        # Put a zero on the respective Dirichlet DOF
        #for index in dirich_ind:
        #    T[index,index] = 0.

        # Generate the new augumented stiffness matrix
        #K_new = np.dot(T, K)
        K_new = np.copy(K)
        K_new[dirich_ind,:] = 0.
        K_new[:,dirich_ind] = 0.

        for index in dirich_ind:
            K_new[index,index] = 1.
            P_new[index] = V[index]

        # LU decomposition
        LU, piv = lu_factor(K_new)
        # Solve the augumented system
        V_res = lu_solve((LU, piv), P_new, trans=0)

        # Obtain reactions at the DOF constrained with Dirichlet BCs
        # (Take the element loads into account with 'P_e')
        P_react = np.dot(K, V_res)[dirich_ind] - P_e[dirich_ind]

        # Copy the data of the model
        model_data = self._model.export_model_data()

        # Create results object
        result = Result(model=model_data)
        result._V = V_res

        # Make dictionary with nodes and respective node reactions
        for ix, index_r in enumerate(dirich_ind):
            node_i, dof_i = self._model.get_node_and_dof(index_r)
            # Add the reactions to the dictionary of reactions of the
            # corresponding node
            node_i.reactions[dof_i] = P_react[ix]

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
            else:
                print('Post-processing of '+curr_output+' not implemented yet.')

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
            axial += np.ones(len(x_axis)) * end_forces[num][0]
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

    def _calc_end_forces(self, result):
        """Calculate the internal forces of beam elements

        :result: Result object
        :returns: TODO

        """
        # Assemble the vector of applied element loads
        # FIXME: only take the contribution of the analyzed beam element
        # and not the cotiguous too
        P_e = self._model._generate_element_loading_vector()
        # calculate for each element
        for num, elem in self._model.beams.items():
            # Get the transformation matrix for the element
            T = elem.transformation_matrix
            # Get the stiffness matrix of the element in global coordinates
            Ke = elem._Ke
            # Get the displacements of the corresponding DOFs in global coordinates
            dof_pn = elem._dof_per_node
            v_i = np.zeros(dof_pn*elem._n_nodes)
            P_e_i = np.zeros(dof_pn*elem._n_nodes)
            # Assemble matrix with element nodal displacements of the current
            # beam element
            for n_node, node in enumerate(elem._nodes):
                i1 = n_node*dof_pn
                i2 = dof_pn*(n_node+1)
                j1 = node.number*dof_pn
                j2 = dof_pn*(1+node.number)
                v_i[i1:i2] = result._V[j1:j2]
                # FIXME:
                if len(elem._loads) > 0:
                    P_e_i[i1:i2] = P_e[j1:j2]

            # Get the End Forces of the element in global coordinates
            P_i = np.dot(Ke, v_i) - P_e_i
            # Transform End Forces to local coordinates
            P_i_local = np.dot(T, P_i)
            # Add End Forces to the dictionary of End Forces of the result object
            # for the current beam element
            result.end_forces[num] = P_i_local

        return result.end_forces

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
            return x * p1
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
        

class SymbolicSolver(Solver):

    """Symbolic solver, meant to be used with the sympy library. It solves linear
    problems in a symbolic manner.
    """

    def __init__(self, model, list_sym, **kwargs):
        """TODO: to be defined1. """
        Solver.__init__(self, model, **kwargs)
        self.list_sym = list_sym

    def solve(self):
        """Solves the problem
        :returns: numpy array

        """
        from sympy.solvers.solveset import linsolve
        from sympy.matrices import Matrix, eye
        from sympy import Symbol
        list_sym = self.list_sym
        # Assemble the stiffness matrix of the system
        K = self._model._assemble_global_sym_K()
        # Assemble the vector of applied loads
        P = self._model._generate_loading_vector()
        # Apply Dirichlet boundary conditions
        # Get the dof corresponding to Dirichlet border conditions
        V = self._model._generate_displacement_vector()
        dirich_ind = self._model._dof_dirichlet
        # Calculate forces generated by the BC
        P_dirich = K[dirich_ind,:].multiply(V)
        P_new = P
        # FIXME: is there a better way of doing this?
        for ix, index in enumerate(dirich_ind):
            P_new[index] = P[index] - P_dirich[ix]
        # Matrix to multiply the global stiffness matrix and put zeros on the
        # rows and columns corresponding to the Dirichlet BC.
        T = eye(len(P))
        # Put a zero on the respective Dirichlet DOF
        for index in dirich_ind:
            T[index,index] = 0.

        # Generate the new augumented stiffness matrix
        K_new = T.multiply(K)

        for index in dirich_ind:
            K_new[index,index] = 1.
            P_new[index] = V[index]

        # Solve the augumented system
        v_sym = []
        for ix in range(len(P_new)):
            v_aux = Symbol('v{n}'.format(n=ix))
            v_sym.append(v_aux)

        system = K_new, P_new
        sym_res = linsolve(system, list_sym )
        # Convert to sympy matrix
        V_res = Matrix(sym_res.args).T

        # Obtained reactions at the DOF constrained with Dirichlet BCs
        P_react = K.multiply(V_res)#[dirich_ind]

        # Copy the data of the model
        model_data = self._model.export_model_data()

        # Create results object
        result = Result(model=model_data)
        result._V = V_res

        # Make dictionary with nodes and respective node reactions
        for ix, index_r in enumerate(dirich_ind):
            node_i, dof_i = self._model.get_node_and_dof(index_r)
            # Add the reactions to the dictionary of reactions of the
            # corresponding node
            node_i.reactions[dof_i] = P_react[ix]

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
            else:
                print('Post-processing of '+curr_output+' not implemented yet.')

        return result

    def _calc_internal_forces(self, result):
        """Calculate member forces

        :result: TODO
        :returns: TODO

        """
        from sympy.matrices import ones, Matrix
        # First calculate end forces on each member
        end_forces = self._calc_end_forces(result)

        max_N = 0.
        max_V = 0.
        max_M = 0.

        # FIXME: still need to include distributed loads, etc
        for num, elem in self._model.beams.items():
            x_axis = Matrix(np.linspace(0,1,11) * elem._length)
            # Calculate axial force on the member (a function is returned)
            axial = ones(len(x_axis)) * -end_forces[num][0]
            # Calculate shear force on the member
            shear = ones(len(x_axis)) * -end_forces[num][1]
            # Calculate moment on the member
            moment = -shear*x_axis - ones(len(x_axis),1)*end_forces[num][2]
            # Add to the results
            res = {'axial':axial, 'shear':shear, 'moment':moment, 'x':x_axis}
            result.internal_forces[num] = res
            # maximum absolute values
            #max_N = max([max(abs(axial)), max_N])
            #max_V = max([max(abs(shear)), max_V])
            #max_M = max([max(abs(moment)), max_M])

        #result._max_member_force['axial'] = max_N
        #result._max_member_force['shear'] = max_V
        #result._max_member_force['moment'] = max_M

        return result.internal_forces

    def _calc_end_forces(self, result):
        """Calculates the internal forces of beam elements

        :result: Result object
        :returns: TODO

        """
        from sympy.matrices import zeros, Matrix
        # calculate for each element
        for num, elem in self._model.beams.items():
            # Get the transformation matrix for the element
            T = Matrix(elem.transformation_matrix)
            # Get the stiffness matrix of the element in global coordinates
            Ke = elem._Ke
            # Get the displacements of the corresponding DOFs in global coordinates
            dof_pn = elem._dof_per_node
            # Initialize matrix for the results
            v_i = zeros(dof_pn*elem._n_nodes, 1)
            # Assemble matrix with element nodal displacements of the current
            # beam element
            for n_node, node in enumerate(elem._nodes):
                i1 = n_node*dof_pn
                i2 = dof_pn*(n_node+1)
                j1 = node.number*dof_pn
                j2 = dof_pn*(1+node.number)
                v_i[i1:i2,:] = result._V[j1:j2]

            # Get the End Forces of the element in global coordinates
            P_i = Ke.multiply(v_i)
            # Transform End Forces to local coordinates
            P_i_local = T.multiply(P_i)
            # Add End Forces to the dictionary of End Forces of the result object
            # for the current beam element
            result.end_forces[num] = P_i_local

        return result.end_forces

class Result(object):

    """An object to store the data obtained from the solving process"""

    def __init__(self, model):
        """TODO: to be defined1. """
        self._V = None # displacement results
        self.end_forces = dict() # end forces of elements
        self.internal_forces = dict() # internal forces of elements
        # maximum absolut value of the different member forces in the
        # model
        self._max_member_force = dict()
        # Data of the model associated with the results
        self._model = model



