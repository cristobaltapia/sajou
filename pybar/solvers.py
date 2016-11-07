#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import dsolve
import pandas as pd
from scipy.linalg import lu_factor, lu_solve
from .postprocessing import Postprocess
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
        # Assemble the vector of applied loads (FIXME)
        P = self._model._generate_loading_vector()
        # Assemble the vector of applied element loads
        P_e = self._model._generate_element_loading_vector()

        # Apply Dirichlet boundary conditions
        # Get the dof corresponding to Dirichlet border conditions
        dirich_ind = self._model._dof_dirichlet
        V = self._model._generate_displacement_vector()
        # Calculate forces generated by the BC
        P_dirich = np.dot(K.todense()[dirich_ind,:], V)
        # Create new loading vector, including the forces generated by
        # the Dirichlet BC and element forces:
        P_new = P[:]
        # - Dirichlet BC
        P_new[dirich_ind] = P[dirich_ind] - P_dirich
        # - Add the element loads to the global loading vector
        P_new += P_e

        # Generate the new augmented stiffness matrix
        # Takes into account the Dirichlet BC
        K_new = np.copy(K.todense())
        K_new[dirich_ind,:] = 0.
        K_new[:,dirich_ind] = 0.

        for index in dirich_ind:
            K_new[index,index] = 1.
            P_new[index] = V[index]

        # Convert to sparse matrix
        K_new = sparse.csr_matrix(K_new)

        # Solve the augmented system (nodal displacements are computed)
        V_res = dsolve.spsolve(K_new, P_new, use_umfpack=False)

        # Initialize a Result object
        # - Copy the data of the model
        model_data = self._model.export_model_data()
        # - Create Result object
        result = Result(model=model_data)

        # Add nodal displacements to the result object
        nodal_displ = self.process_nodal_displ(V_res)
        result.add_result('nodal displacements', nodal_displ)

        # Add nodal reactions to the results object
        self.calc_nodal_reactions(result, V_res, K, P_e)

        # Add nodal forces to the result object
        self.calc_nodal_forces(result, V_res, K, P_e)

        # Add end forces to the results object
        self.calc_end_forces(result, V_res)

        # Postprocess the results according to the specified in 'output'
        # variable
        result = self.postprocess(result)

        return result

    def postprocess(self, result):
        """Postprocess the specified results given in the 'output' variable

        :result: TODO
        :returns: TODO

        """
        #
        postp = Postprocess(result)
        #
        for curr_output in self._output_request:
            if curr_output == 'internal forces':
                # Calculate forces along the frame elements
                postp.calc_all_internal_forces()
            elif curr_output == 'deflections':
                postp.calc_all_deflections()
            # FIXME: not sure if put this here
            elif curr_output == 'stresses':
                postp.calc_stresses(result)
            else:
                print('Post-processing of '+curr_output+' not implemented yet.')

        return result
    
    def process_nodal_displ(self, nodal_displ):
        """Return nodal displacements as a pandas DataFrame instance.
        Each column represents a DOF and the index correspond to the respective node
        number.

        :nodal_displ: TODO
        :returns: TODO

        """
        n_dimensions = self._model.n_dimensions

        nodes = [n for i, n in self._model.nodes.items()]

        # Initialize numpy array
        ar_nodal_displ = np.zeros((len(nodes), n_dimensions), dtype=np.float64)
        index_nodes = np.zeros(len(nodes), dtype=np.int)

        # Loop for each node of the model
        # TODO: implement 3D case
        for ix, curr_node in enumerate(nodes):
            ix_node = curr_node.number
            aux_arr = np.array([nodal_displ[ix_node*3], nodal_displ[ix_node*3+1]])
            ar_nodal_displ[ix,:] = aux_arr
            index_nodes[ix] = curr_node.number

        index_label = ['x','y']

        # Create the data frame
        df_nodal_displ = pd.DataFrame(data=ar_nodal_displ, index=index_nodes, dtype=np.float64,
                columns=index_label)

        return df_nodal_displ

    def calc_nodal_forces(self, result, nodal_displ, K, elem_load):
        """TODO: Docstring for calc_nodal_forces.

        :result: TODO
        :nodal_displ: TODO
        :K: TODO
        :elem_load: TODO
        :returns: TODO

        """
        nodal_forces = K.dot(nodal_displ) - elem_load
        #
        n_dof = self._model.n_dof_per_node

        nodes = [n for i, n in self._model.nodes.items()]

        # Initialize numpy array
        ar_nodal_forces = np.zeros((len(nodes), n_dof), dtype=np.float64)
        index_nodes = np.zeros(len(nodes), dtype=np.int)

        # Loop for each node of the model
        # TODO: implement 3D case
        for ix, curr_node in enumerate(nodes):
            ix_node = curr_node.number
            aux_arr = np.array([nodal_forces[ix_node*3: ix_node*3+3]])
            ar_nodal_forces[ix,:] = aux_arr
            index_nodes[ix] = curr_node.number

        index_label = ['fx', 'fy', 'mz']

        # Create the data frame
        df_nodal_forces = pd.DataFrame(data=ar_nodal_forces, index=index_nodes, dtype=np.float64,
                columns=index_label)

        result.add_result('nodal forces', df_nodal_forces)

        return df_nodal_forces

    def calc_nodal_reactions(self, result, nodal_displ, K, elem_load):
        """Calculate the nodal reactions of the model.

        :result: TODO
        :nodal_displ: nodal displacements
        :K: stiffness matrix
        :elem_forces: element forces vector
        :returns: TODO

        """
        # Get the positions of the Dirichlet birder conditions
        dirich_ind = result._model._dof_dirichlet
        # (Take the element loads into account with 'P_e')
        nodal_react = K.dot(nodal_displ)[dirich_ind] - elem_load[dirich_ind]
        # Add result to the Result object
        result.add_result('nodal reactions', nodal_react)

        # Make dictionary with nodes and respective node reactions
        for ix, index_r in enumerate(dirich_ind):
            node_i, dof_i = self._model.get_node_and_dof(index_r)
            # Add the reactions to the dictionary of reactions of the
            # corresponding node
            node_i.reactions[dof_i] = nodal_react[ix]

        return nodal_react

    def calc_end_forces(self, result, nodal_displ):
        """Calculate the internal forces of beam elements.

        :result: Result object
        :returns: TODO

        """
        # Initialize dictionary with results of the end forces
        end_forces = dict()
        # calculate for each element
        for num, elem in self._model.beams.items():
            # Get the transformation matrix for the element
            T = elem.transformation_matrix
            # Get the stiffness matrix of the element in global coordinates
            # FIXME!
            Ke = elem._Ke
            #Ke = elem.__assemble_Ke__()
            # Get number of DOF per node
            dof_pn = elem._dof_per_node
            # Get the displacements of the corresponding DOFs in global coordinates
            # - Initialize matrix
            v_i = np.zeros(dof_pn * elem._n_nodes, dtype=np.float64)
            P_e_i = np.zeros(dof_pn * elem._n_nodes, dtype=np.float64)
            # Assemble matrix with element nodal displacements of the current
            # beam element
            for n_node, node in enumerate(elem._nodes):
                # indices for the array v_i
                i1 = dof_pn * n_node
                i2 = dof_pn * (n_node + 1)
                # Indices corresponding to the position of the DOF of
                # the current node analyzed
                j1 = dof_pn * node.number
                j2 = dof_pn * (1 + node.number)
                # Add the results for these DOFs to the v_i array
                v_i[i1:i2] = nodal_displ[j1:j2] # DOF of node selected
            """
            If an element uses release_end option, then the displacement (rotation) needs
            to be calculated sepparately.
            """
            if elem.release_end_1 == True or elem.release_end_2 == True:
                aux = elem._calc_condensed_displacements(v_i)
                #v_i[5] = aux

            if len(elem._loads) > 0:
                # Element load vector FIXME:
                P_e_i = elem._loads[0]._load_vector_global
            else:
                P_e_i = 0.

            # Get the End Forces of the element in global coordinates
            P_i_global = Ke.dot(v_i) - P_e_i
            # Transform End Forces to local coordinates
            P_i_local = T.dot(P_i_global)
            # Add End Forces to the dictionary of End Forces of the result object
            # for the current beam element
            end_forces[num] = P_i_local
            # Add the corresponding nodal forces to the matrix to
            # calculate the sectional forces
            # - Axial force
            elem._poly_sec_force[0,0] += -P_i_local[0]
            # - Shear force
            elem._poly_sec_force[0,1] += P_i_local[1]
            # - Moment
            elem._poly_sec_force[0,2] += -P_i_local[2]
            elem._poly_sec_force[1,2] += P_i_local[1]

        # Add results to the result object
        result.add_result('end forces', end_forces)

        return end_forces


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
        # Analyzed model
        self._model = model
        # Initiate the container for the results
        self.data = dict()
        # Initialize dictionary with information of minimum and maximum
        # values of the results, and other useful data.
        self.metadata = dict()

    def add_result(self, name, results):
        """Add results to the results dictionary

        :name: TODO
        :results: TODO
        :returns: TODO

        """
        self.data[name] = results

        return self

    def add_metadata(self, name, results):
        """Add results to the results dictionary

        :name: TODO
        :results: TODO
        :returns: TODO

        """
        self.metadata[name] = results

        return self


