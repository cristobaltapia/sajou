#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from sajou.postprocessing import Postprocess
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import dsolve

""" This module contains the different solvers used by the module."""


class Solver(object):
    """
    Parent class for the solvers used

    Parameters
    ----------
    model: TODO
    **kwargs: Additional options:

        - output: list of keywords specifying the different results that need to be
            computed after the solution is found.
            Accepted values are:
            - 'reactions'
            - 'internal forces'
            - 'stresses'
            - 'strains'
            - 'energy'

    Attributes
    ----------
    nfat: dict
        Node Freedom Allocation Table
    n_dof_per_node: int
        number of degrees of freedom per node
        spacial dimensions used in the model
    _K: numpy ndarray
        global stiffness matrix
    _P: numpy ndarray
        global load vector
    _V: numpy ndarray
        global displacement vector
    _dof_dirichlet: list
        number of degrees of freedom with Dirichlet border conditions
    _nfmt: dict
        Node Freedom Map Table

        """

    def __init__(self, model, **kwargs):
        self._model = model
        # The variable 'output' stores the different results that need
        # to be computed after the solution is found.
        # FIXME: implement this
        deffault_output = ['displacements', 'forces']
        self._output_request = kwargs.get('output', None)
        self.nfat = dict()
        # Specify dofs that are not active due to border conditions
        self._dof_dirichlet = []

    def get_node_and_dof(self, dof):
        """Return the node and element dof (number of the dof in a specific element)
        corresponding to the global dof given.

        Parameters
        ----------
        dof: int
            global *dof*

        Returns
        -------
        sajou.Node
            Node correpsonding to the *dof* specified
        int
            number of the *dof* of the element corresponding to the global
            *dof* supplied


        """
        # FIXME: node freedom allocation table should be used here
        # Get the node
        nodes = self._model.nodes
        nfmt = self.nfmt

        # search for the node
        for n, node in nodes.items():
            n_dof_node = sum(node.nfs)
            if dof >= nfmt[node.number] and dof < (nfmt[node.number] + n_dof_node):
                node_num = node.number
                break

        # search for the node in the element
        node_sel = nodes[node_num]
        # Get the intra-element dof
        n_dof = dof - nfmt[node_num]

        return node_sel, n_dof

    def _assemble_global_K(self):
        """Assemble the global stiffness matrix using the addition method.

        Returns
        -------
        numpy.array
           Global stiffness matrix of the system

        """
        elements = self._model.elements
        # Generate Node Freedom Al[location Table and total number of
        # active DOFs of the system
        self.__generate_node_freedom_allocation_dict__()
        # Generate Node Freedom Map dictionary
        self.__generate_node_freedom_map_dict__()
        # number of dof per node
        n_dof = self.n_active_dof
        # Initialize the global stiffness matrix
        K = np.zeros([n_dof, n_dof], dtype=np.float64)
        # Fill the global stiffness matrix
        for n_elem, element in elements.items():
            # Get nodes of the respective element
            nodes_elem = element._nodal_connectivity
            # Assemble element stiffness matrix
            element.assemble_Ke()
            # List to store the global system indices
            g_i = []
            # List of indices of used element DOFs
            g_e = []
            # For each node in the current element
            for n_node_e, node in nodes_elem.items():
                # Get Element Freedom Signature
                efs = element.efs[n_node_e]
                # Number of active DOF in the node of the element
                active_dof = np.sum((efs >= 1))
                if active_dof > 0:
                    # Get value of th Node Freedom Assign Table for the
                    # current node
                    nfat_node = self.nfmt[node.number]
                    # Get NFS of the node in the element
                    enfs_node = element.enfmt[n_node_e]
                    # for the total of used DOF in the node
                    # FIXME!!
                    index_base = element.get_node_active_dof(n_node_e)
                    active_nodes = nfat_node + index_base
                    # Extend the list
                    g_i.extend(active_nodes)
                    #
                    index_base_e = element.get_element_active_dof(n_node_e)
                    active_nodes_e = enfs_node + index_base_e
                    g_e.extend(active_nodes_e)

            # Convert list to numpy array in order to broadcast more
            # easily to the global stiffness matrix
            g_i = np.array(g_i)
            g_e = np.array(g_e)
            # Add the contributions to the respective DOFs in global system
            K[g_i[:, None], g_i] += element._Ke[g_e[:, None], g_e]

        # Generate sparse matrix
        K_s = sparse.csr_matrix(K)

        return K_s

    def _generate_loading_vector(self):
        """Generate the global matrix of applied forces P.

        Returns
        -------
        numpy.array
            Loading vector of the system

        """
        nodes = self._model.nodes
        # number of dof per node
        n_dof = self.n_active_dof
        # Get the node freedom allocation map table
        nfmt = self.nfmt
        # Initialize a zero vector of the size of the total number of
        # dof
        P = np.zeros(n_dof, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node in nodes.items():
            # Get the Node Freedom Signature of the current node
            nfs = node.nfs
            #
            index_i = np.array([kx for kx in node._loads.keys()],
                               dtype=np.int) + nfmt[node.number]
            P[index_i] = np.array([kx for kx in node._loads.values()])

        self._P = P

        return P

    def _generate_element_loading_vector(self):
        """Generate the global element vector of forces.

        Returns
        -------
        numpy.array
            Loading vector of the system

        """
        elements = self._model.elements
        # number of dof per node
        n_dof = self.n_active_dof
        # Initialize a zero vector of the size of the total number of
        # DOF
        P = np.zeros(n_dof, dtype=np.float64)

        # Add loads applied to the elements (distributed loads)
        for ix, element in elements.items():
            # Check if the element has element loads defined
            if len(element._loads) > 0:
                # Get nodes of the respective element
                nodes_elem = element._nodal_connectivity
                # List of indices of the global system
                g_i = []
                # List of indices of active element DOFs
                g_e = []
                # For each node in the current element
                for n_node_e, node in nodes_elem.items():
                    # Get Element Freedom Signature
                    efs = element.efs[n_node_e]
                    # Number of active DOF in the node of the element
                    active_dof = np.sum((efs >= 1))

                    if active_dof > 0:
                        # Get value of th Node Freedom Assign Table for the
                        # current node
                        nfat_node = self.nfmt[node.number]
                        # Get node freedom signature of the node in the element
                        enfs_node = element.enfmt[n_node_e]
                        # Get the corresponding active indices of the
                        # node in the element
                        index_base = element.get_node_active_dof(n_node_e)
                        active_nodes = nfat_node + index_base
                        # Extend the list
                        g_i.extend(active_nodes)
                        #
                        index_base_e = element.get_element_active_dof(
                            n_node_e)
                        active_nodes_e = enfs_node + index_base_e
                        g_e.extend(active_nodes_e)

                # Add to the global load vector
                P[g_i] += element._load_vector_e[g_e]

        self._P = P

        return P

    def _generate_displacement_vector(self):
        """Generate the global displacement matrix V, containing the border
        conditions applied to each dof.

        Returns
        -------
        numpy.array
            Dsiplacement vector of the system

        """
        nodes = self._model.nodes
        # number of dof per node
        n_dof = self.n_active_dof
        # Get the node freedom allocation map table
        nfmt = self.nfmt
        # Initialize a zero vector of the size of the total number of
        # dof
        V = np.zeros(n_dof, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node in nodes.items():
            # Get the Node Freedom Signature of the current node
            nfs = node.nfs
            #
            index_i = np.array([kx for kx in node._bc.keys()],
                               dtype=np.int) + nfmt[node.number]
            V[index_i] = np.array([kx for kx in node._bc.values()])
            # Add to the list of restrained DOFs
            self._dof_dirichlet.extend(index_i.tolist())

        self._V = V

        return V

    def __generate_node_freedom_allocation_dict__(self):
        """
        Generate the Node Freedom Allocation Table.

        Generates a dictionary of arrays, containing the Node Freedom
        Signature of each Node of the model. The key values correspond to the
        node numbers.
        Also counts the total number of active DOFs of the system and stores
        it on the variable self.n_active_dof

        Returns
        -------
        dict
            the Node Freedom Allocation Table of the system

        """
        nodes = self._model.nodes
        n_active_dof = 0
        # Loop over each node and add its Node Freedom Signature to the
        # dictionary
        for n_num, node in nodes.items():
            node.__generate_node_freedom_signature__()
            self.nfat[node.number] = node.nfs
            n_active_dof += np.sum(node.nfs)

        self.n_active_dof = n_active_dof

        return self.nfat

    def __generate_node_freedom_map_dict__(self):
        """
        Generate the Node Freedom Map Table of the system.

        The Node Freedom Map Table is a dictionary that contains the index,
        relative to the global system, to which each node's first active DOF
        contributes.
        It is assumed that the Node Freedom Allocation Table has already been
        generated using the function __generate_node_freedom_allocation_dict__().

        Returns
        -------
        numpy.array
            TODO

        """
        nodes = self._model.nodes
        # Obtain the number of active DOFs in each node:
        n_active_dof = {n: sum(node.nfs) for n, node in nodes.items()}
        # get the list of the node numbers
        n_num = [n for n, node in nodes.items()]
        # order the list
        n_sorted = np.sort(n_num)
        # Obtain the cumulative sum
        nfmt = dict()
        cumsum = 0
        for n in n_sorted:
            cumsum += n_active_dof[n]
            nfmt[n] = cumsum - n_active_dof[n]

        self.nfmt = nfmt

        return nfmt


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
        K = self._assemble_global_K()
        # Assemble the vector of applied loads (FIXME)
        P = self._generate_loading_vector()
        # Assemble the vector of applied element loads
        P_e = self._generate_element_loading_vector()

        # Apply Dirichlet boundary conditions
        # Get the dof corresponding to Dirichlet border conditions
        V = self._generate_displacement_vector()
        dirich_ind = self._dof_dirichlet
        # Calculate forces generated by the BC
        P_dirich = K.todense()[dirich_ind, :] @ V
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
        K_new[dirich_ind, :] = 0.
        K_new[:, dirich_ind] = 0.

        for index in dirich_ind:
            K_new[index, index] = 1.
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
        self.calc_nodal_displ(result, V_res, K, P_e)

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
                print('Post-processing of ' + curr_output +
                      ' not implemented yet.')

        return result

    def calc_nodal_displ(self, result, nodal_displ, K, elem_load):
        """Calculate the nodal displacements

        Parameters
        ----------

        result: Results object
            result obtained from the Solver
        nodal_displ: ndarray
            displacements array calculated after the solve process
        K: ndarray
            Stiffness matrix
        elem_load: ndarray
            Elemetn load

        """
        # generate a list of the nodes of the system
        nodes = [n for i, n in self._model.nodes.items()]
        # get the node freedom map table of the model
        nfmt = self.nfmt
        # get the node freedom allocation table
        nfat = self.nfat

        # initialize empty dictionary
        dict_displ = dict()

        # Loop for each node of the model
        for node_i in nodes:
            # the the index of the node
            ix_node = node_i.number
            # get the indices for the node's DOFs on the global system
            ix_i = nfmt[ix_node]
            ix_f = ix_i + sum(node_i.nfs)
            # get the displacements for the specified node
            dict_displ[ix_node] = nodal_displ[ix_i: ix_f]

        result.add_result('nodal displacements', dict_displ)

        return dict_displ

    def calc_nodal_forces(self, result, nodal_displ, K, elem_load):
        """Compute the nodal forces

        Parameters
        ----------
        result: Result object
            result obtained from the solver
        nodal_displ: ndarray
            nodal displacements calculated with the stiffness matrix
        K: ndarray
            Stiffness matrix
        elem_load: TODO

        Returns
        -------
        dict: dictionary with the nodal forces

        """
        # calculate the nodal forces
        nodal_forces = K @ nodal_displ - elem_load

        # get the node freedom map table
        nfmt = self.nfmt
        nodes = self._model.nodes

        dict_nf = dict()
        # Loop for each node of the model
        for num, curr_node in nodes.items():
            # get total DOFs of the node
            ndof = curr_node.n_dof
            # get node number
            ix = curr_node.number
            # get node freedom signature
            nfs = curr_node.nfs
            # initialize nodal force vector
            nod_force = np.zeros(ndof)
            #
            indices = np.arange(ndof)[nfs>0]
            nod_force[indices] = nodal_forces[nfmt[ix]: nfmt[ix] + sum(nfs)]
            dict_nf[ix] = nod_force

        result.add_result('nodal forces', dict_nf)

        return dict_nf

    def calc_nodal_reactions(self, result, nodal_displ, K, elem_load):
        """Calculate the nodal reactions of the model.

        Parameters
        ----------

        result: Result object
            result obtained with Solver
        nodal_displ: ndarray
            nodal displacements
        K: ndarray
            stiffness matrix
        elem_forces: ndarray
            element forces vector

        Returns
        -------

        ndarray: nodal reactions

        """
        # Get the positions of the Dirichlet birder conditions
        dirich_ind = self._dof_dirichlet
        # (Take the element loads into account with 'P_e')
        nodal_react = (K @ nodal_displ)[dirich_ind] - elem_load[dirich_ind]
        # Add result to the Result object
        result.add_result('nodal reactions', nodal_react)

        # Make dictionary with nodes and respective node reactions
        for ix, index_r in enumerate(dirich_ind):
            node_i, dof_i = self.get_node_and_dof(index_r)
            # Add the reactions to the dictionary of reactions of the
            # corresponding node
            node_i.reactions[dof_i] = nodal_react[ix]

        return nodal_react

    def calc_end_forces(self, result, nodal_displ):
        """
        Calculate the internal forces of elements.

        Parameters
        ----------

        result: Result object
            used to store the results to
        nodal_displ: ndarray
            nodal displacements obtained from the solver

        Returns
        -------

        dict: End forces

        """
        elements = self._model.elements
        # Initialize dictionary with results of the end forces
        end_forces = dict()
        # Get the node freedom map table
        nfmt = self.nfmt
        # Get the node freedom allocation of every node
        nfat = self.nfat
        # Calculate end forces for each element
        for num, element in elements.items():
            # Get the transformation matrix for the element
            Te = element.transformation_matrix
            # Get the stiffness matrix of the element in global coordinates
            # FIXME!
            Ke = element._Ke
            # Get the displacements of the corresponding DOFs in global coordinates
            # - Initialize matrix
            v_i = np.zeros(element._k_size, dtype=np.float64)
            P_e_i = np.zeros(element._k_size, dtype=np.float64)
            # Get the node freedom allocation map table of the current element
            enfmt = element.enfmt
            # Initiate element load vector
            # Assemble matrix with element nodal displacements of the current
            # element
            for n_node_e, node in element._nodal_connectivity.items():
                # Indices for the array v_i
                index_base = element.get_element_active_dof(n_node_e)
                i_index = index_base + enfmt[n_node_e]
                # Indices corresponding to the position of the DOF of
                # the current node analyzed (global system)
                index_base_n = element.get_node_active_dof(n_node_e)
                j_index = index_base_n + nfmt[node.number]
                # Add the results for these DOFs to the v_i array
                v_i[i_index] = nodal_displ[j_index]  # DOF of node selected
                # Account for element loads
                if len(element._loads) > 0:
                    # Same indices are used here
                    for load in element._loads:
                        P_e_i[i_index] += load._load_vector_global[i_index]

            # Get the End Forces of the element in global coordinates
            P_i_global = Ke @ v_i - P_e_i
            # Transform End Forces to local coordinates
            P_i_local = Te @ P_i_global
            # Add End Forces to the dictionary of End Forces of the result object
            # for the current element
            end_forces[element.number] = P_i_local
            # Add the corresponding nodal forces to the matrix to
            # calculate the sectional forces
            # - Axial force
            element._poly_sec_force[0, 0] += -P_i_local[0]
            # - Shear force
            element._poly_sec_force[0, 1] += P_i_local[1]
            # - Moment
            element._poly_sec_force[0, 2] += -P_i_local[2]
            element._poly_sec_force[1, 2] += P_i_local[1]

        # Add results to the result object
        result.add_result('end forces', end_forces)

        return end_forces


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
