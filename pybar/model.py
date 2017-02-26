#!/usr/bin/env python
# encoding: utf-8
"""
Defines the model classes for 2D and 3D models.
"""
__docformat__ = 'reStructuredText'

import numpy as np
import pandas as pd
import scipy.sparse as sparse

from .materials import Material
from .nodes import Node2D
from .sections import BeamSection

class Model(object):
    """Defines a model object"""

    def __new__(cls, name, dimensionality):
        if cls is Model:
            if dimensionality == '2D':
                return super(Model, cls).__new__(Model2D)
            if dimensionality == '3D':
                return super(Model, cls).__new__(Model3D)
        else:
            return super(Model, cls).__new__(cls, name, dimensionality)

    def __init__(self, name, dimensionality):
        """Initialization of model instance.

        :name: TODO
        :dimensionality: TODO

        """
        self._name = name
        self._dimensionality = dimensionality
        # Node Freedome Allocation Table:
        #
        self.nfat = dict()
        self.nodes = dict()
        self.beams = dict()
        self.beam_sections = dict()
        self.materials = dict()
        self.n_nodes = 0
        self.n_beams = 0
        self.n_materials = 0
        self._connectivity = None
        self._K = None  # global stiffness matrix
        self._P = None  # load matrix
        self._V = None  # global displacement matrix
        # Number of dimensions of the model
        self.n_dimensions = None
        # Number of dof per node. Initialized in the respective models
        self.n_dof_per_node = None
        # Specify dofs that are not active due to border conditions
        self._dof_dirichlet = []
        # Node Freedom Allocation Table
        self._nfat = dict()
        # Node Freedom Map Table:
        # Stores the index of the first used DOF of a node in the global
        # system.
        self._nfmt = dict()

    def Material(self, name, data, type='isotropic'):
        """Function used to create a Material instance in the model

        :param name: name of the material
        :param data: data for the material
        :param type: type of the material
        :returns: a Material instance

        """
        material = Material(name=name, data=data, type=type)
        # Add the material to the dictionary of materials in the current
        # model
        self.materials[name] = material
        self.n_materials += 1

        return material

    def BeamSection(self, name, material, data, type='rectangular'):
        """Function use to create a BeamSection instance in the model

        :param name: name of the section
        :param material: material for the section
        :param data: data (see BeamSection class definition)
        :param type: type of the section (see BeamSection class definition)
        :returns: a beam section instance

        """
        # The material can be passed both as a string, corresponding to
        # a key of the material dictionary of the model, or as a
        # material instance directly.

        if isinstance(material, str):
            material_section = self.materials[material]
        else:
            material_section = material

        section = BeamSection(
            name=name, material=material_section, data=data, type=type)
        # Add section to the list of beam sections
        self.beam_sections[name] = section

        return section

    def __str__(self):
        """
        Printable string
        """
        return str(
            'Model: Name: {name}, Nodes: {n_nodes}, Beams: {n_beams}'.format(
                name=self._name, n_nodes=self.n_nodes, n_beams=self.n_beams))

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return str(
            'Model: Name: {name}, Nodes: {n_nodes}, Beams: {n_beams}'.format(
                name=self._name, n_nodes=self.n_nodes, n_beams=self.n_beams))

    def _generate_connectivity_table2D(self):
        """Generates the connectivity table for the model
        :returns: numpy array


                element | N1 | N2
                ------------------
        row -->   0     | 0  |  1
                  1     | 1  |  2
                  2     | 0  |  3
                  .       .     .
                  .       .     .
                  .       .     .
                  9     | 3  |  4

        """
        # Connectivity matrix for the Beams
        conn_matrix = np.zeros((len(self.beams), 3), dtype=np.float64)
        #
        count = 0
        for num, curr_line in self.beams.items():
            conn_matrix[count, 0] = curr_line.number
            conn_matrix[count, 1] = curr_line._node1.number
            conn_matrix[count, 2] = curr_line._node2.number
            count += 1

        self._connectivity = conn_matrix

        return conn_matrix

    def _assemble_global_K(self):
        """Assemble the global stiffness matrix using the addition method.
        :returns: numpy array

        """
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
        for n_elem, element in self.beams.items():
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
                active_dof = np.sum(efs)
                if active_dof > 0:
                    # Get value of th Node Freedom Assign Table for the
                    # current node
                    nfat_node = self.nfmt[node.number]
                    # Get NFS of the node in the element
                    enfs_node = element.enfmt[n_node_e]
                    # for the total of used DOF in the node
                    index_base = element.get_index_array_of_node(n_node_e)
                    active_nodes = nfat_node + index_base
                    # Extend the list
                    g_i.extend(active_nodes)
                    #
                    index_base_e = element.get_index_array_of_node(n_node_e)
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
        :returns: numpy array

        """
        # number of dof per node
        n_dof = self.n_active_dof
        # Get the node freedom allocation map table
        nfmt = self.nfmt
        # Initialize a zero vector of the size of the total number of
        # dof
        P = np.zeros(n_dof, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node in self.nodes.items():
            # Get the Node Freedom Signature of the current node
            nfs = node.nfs
            #
            index_i = np.array(
                [kx for kx in node._loads.keys()],
                dtype=np.int) + nfmt[node.number]
            P[index_i] = np.array([kx for kx in node._loads.values()])

        self._P = P

        return P

    def _generate_element_loading_vector(self):
        """Generate the global element vector of forces.
        :returns: numpy array

        """
        # number of dof per node
        n_dof = self.n_active_dof
        # Initialize a zero vector of the size of the total number of
        # DOF
        P = np.zeros(n_dof, dtype=np.float64)

        # Add loads applied to the elements (distributed loads)
        for ix, element in self.beams.items():
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
                    active_dof = np.sum(efs)

                    if active_dof > 0:
                        # Get value of th Node Freedom Assign Table for the
                        # current node
                        nfat_node = self.nfmt[node.number]
                        # Get node freedom signature of the node in the element
                        enfs_node = element.enfmt[n_node_e]
                        # Get the corresponding active indices of the
                        # node in the element
                        index_base = element.get_index_array_of_node(n_node_e)
                        active_nodes = nfat_node + index_base
                        # Extend the list
                        g_i.extend(active_nodes)
                        #
                        index_base_e = element.get_index_array_of_node(
                            n_node_e)
                        active_nodes_e = enfs_node + index_base_e
                        g_e.extend(active_nodes_e)

                # Add to the global load vector
                P[g_i] += element._load_vector_e[g_e]

        self._P = P

        return P

    def _generate_displacement_vector(self):
        """This function generates the global displacement matrix V, containing the border
        conditions applied to each dof.

        :returns: numpy array

        """
        # number of dof per node
        n_dof = self.n_active_dof
        # Get the node freedom allocation map table
        nfmt = self.nfmt
        # Initialize a zero vector of the size of the total number of
        # dof
        V = np.zeros(self.n_nodes * self.n_dof_per_node, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node in self.nodes.items():
            # Get the Node Freedom Signature of the current node
            nfs = node.nfs
            #
            index_i = np.array(
                [kx for kx in node._bc.keys()],
                dtype=np.int) + nfmt[node.number]
            V[index_i] = np.array([kx for kx in node._bc.values()])
            # Add to the list of restrained DOFs
            self._dof_dirichlet.extend(index_i.tolist())

        self._V = V

        return V

    def BC(self, node, type='displacement', coord_system='global', **kwargs):
        """
        Introduces a border condition to the node.

        :param node: a Node instance
        :param type: type of border condition

            - Options:
                ``'displacement'``, ``...``
        :param coord_system: spcifies the coordinate system

        **Optional arguments:** (**kwargs) At least one of the following parameters must be supplied

        :param v1: displacement in the direction 1
        :param v2: displacement in the direction 2
        :param v3: displacement in the direction 3
        :param r1: rotation in the direction 1
        :param r2: rotation in the direction 2
        :param r3: rotation in the direction 3
        :type v1: float
        :type v2: float
        :type v3: float
        :type r1: float
        :type r2: float
        :type r3: float

        :returns: TODO

        """
        # TODO: currently only in global coordintaes. Implement
        # transformation in other coordinate systems.

        # Get the BC applied
        v1 = kwargs.get('v1', None)
        v2 = kwargs.get('v2', None)
        v3 = kwargs.get('v3', None)
        r1 = kwargs.get('r1', None)
        r2 = kwargs.get('r2', None)
        r3 = kwargs.get('r3', None)

        # For the case of the 2D model
        if self.n_dof_per_node == 3:
            list_dof = [v1, v2, r3]
            for dof, curr_bc in enumerate(list_dof):
                if curr_bc is not None:
                    node.set_BC(dof=dof, val=curr_bc)

        # For the case of the 3D model
        elif self.n_dof_per_node == 6:
            list_dof = [v1, v2, v3, r1, r2, r3]
            for dof, curr_bc in enumerate(list_dof):
                if curr_bc is not None:
                    node.set_BC(dof=dof, val=curr_bc)

        return None

    # TODO: there has to give a 'Load' class to handle the different
    # type of loads.
    def Load(self, node, coord_system='global', **kwargs):
        """Introduces a Load in the given direction according to the selected coordinate
        system at the specified node.

        :param node: a Node instance
        :param coordinate: coordinate system
        :param **kwargs: optional arguments. The BC is defined for the different degree of freedom (*dof*) available to the node.

        **Optional arguments:** (**kwargs) At least one of the following parameters must be supplied

        :param f1: force in direction 1
        :param f2: force in direction 2
        :param f3: force in direction 3
        :param m1: moment in direction 1
        :param m2: moment in direction 2
        :param m3: moment in direction 3

        :returns: a Load instance

        """
        # TODO: currently only in global coordintaes. Implement
        # transformation in other coordinate systems.

        # Get the BC applied
        f1 = kwargs.get('f1', None)
        f2 = kwargs.get('f2', None)
        f3 = kwargs.get('f3', None)
        m1 = kwargs.get('m1', None)
        m2 = kwargs.get('m2', None)
        m3 = kwargs.get('m3', None)

        # For the case of the 2D model
        if self.n_dof_per_node == 3:
            list_dof = [f1, f2, m3]
            for dof, curr_force in enumerate(list_dof):
                if curr_force is not None:
                    node.set_Load(dof=dof, val=curr_force)

        # For the case of the 3D model
        elif self.n_dof_per_node == 6:
            list_dof = [f1, f2, f3, m1, m2, m3]
            for dof, curr_force in enumerate(list_dof):
                if curr_force is not None:
                    node.set_Load(dof=dof, val=curr_force)

        return None

    def export_model_data(self):
        """Export all the data of the model. This means the nodes, elements,
        border conditions and forces are exported to a ModelData object.
        :returns: TODO

        """
        model_data = ModelData(self)

        return model_data

    def get_node_and_dof(self, dof):
        """Returns the node and element dof (number of the dof in a specific element)
        corresponding to the global dof given.

        :param dof: global *dof*
        :type dof: int
        :returns: node, int

        """
        # Get the node
        n_node = dof // self.n_dof_per_node
        node = self.nodes[n_node]
        # Get the intra-element dof
        n_dof = dof - (n_node * self.n_dof_per_node)

        return node, n_dof

    def add_hinge(self, node):
        """Add hinge to the specified node. Also supports list of nodes

        :param node: Node instance or list of node instances
        :returns: TODO

        """
        #FIXME: not yet implemented!
        if isinstance(node, list):
            for node_i in node:
                node_i.add_hinge()
        else:
            node.add_hinge()

    def __generate_node_freedom_allocation_dict__(self):
        """
        Generate the Node Freedom Allocation Table.

        Generates a dictionary of arrays, containing the Node Freedom
        Signature of each Node of the model. The key values correspond to the node
        numbers.
        Also counts the total number of active DOFs of the system and stores it on the
        variable self.n_active_dof

        :returns: a dictionary

        """
        n_active_dof = 0
        # Loop over each node and add its Node Freedom Signature to the
        # dictionary
        for n_num, node in self.nodes.items():
            node.__generate_node_freedom_signature__()
            self.nfat[node.number] = node.nfs
            n_active_dof += np.sum(node.nfs)

        self.n_active_dof = n_active_dof

        return self.nfat

    def __generate_node_freedom_map_dict__(self):
        """
        Generate the Node Freedom Map Table of the system.

        The Node Freedom Map Table is a dictionary that contains the index, relative to
        the global system, to which each node's first active DOF contributes.
        It is assumed that the Node Freedom Allocation Table has already been generated
        using the function __generate_node_freedom_allocation_dict__().

        :returns: TODO

        """
        from numpy import cumsum
        # Obtain the number of active DOFs in each node:
        n_active_dof = [sum(node.nfs) for n_node, node in self.nodes.items()]
        # Obtain the cumulative sum
        nfmt = cumsum(n_active_dof, dtype=np.int) - n_active_dof[0]
        # TODO: make this a dictionary

        self.nfmt = nfmt

        return nfmt


class Model2D(Model):
    """
    Subclass of the 'Model' class. It is intended to be used for the 2-dimensional
    models of frame structures.

    Allocation of DOFs in each node:

            [1 2 3] = [ux, uy, rz]

    """

    def __init__(self, name, dimensionality='2D'):
        """TODO: to be defined1. """
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)
        # Numer of dimensions
        self.n_dimensions = 2
        # Number of degrees of freedom per node:
        self.n_dof_per_node = 3

    def Node(self, x, y):
        """
        2D implementation of the Node.

        :param x: x position
        :param y: y position
        :type x: float
        :type y: float

        :returns: instance of Node

        """
        # A coordinate z=0 is passed to initiate the Node Instance
        node = Node2D(x=x, y=y, z=0.0, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def Beam(self, node1, node2):
        """
        Define a line between two nodes.

        :param node1: first node
        :param node2: second node

        """
        from .elements.beam2d import Beam2D
        line = Beam2D(node1=node1, node2=node2, number=self.n_beams)
        self.beams[line.number] = line
        self.n_beams += 1

        return line

    def distributed_load(self, elements, **kwargs):
        """
        Add a distributed load to a list of beam elements.
        A list of elements has to be supplied for the first variable. The rest of the
        variables are exactly the same as in the 'distributed_load' function of the
        corresponding elements.

        :param elements: list of beam elements
        :param p1: value of the force at start node
        :param p2: value of the force at end node
        :param direction: direction (default: *2*)
        :param coord_system: coordinate system (default: global)
        :type p1: float
        :type p2: float
        :returns: TODO

        """
        for curr_elem in elements:
            # Add distributed load
            curr_elem.distributed_load(**kwargs)

        return 1


class Model3D(Model):
    """Subclass of the 'Model' class. It is intended to be used for the 3-dimensional
    models of frame structures."""

    def __init__(self, name, dimensionality='3D'):
        """TODO: to be defined1. """
        dimensionality = '3D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 6  # dof per node

    def Node(self, x, y, z):
        """3D implementation of the Node.

        :*kwargs: TODO
        :returns: instance of Node

        """
        node = Node(x=x, y=y, z=z, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def Beam(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Beam3D(node1=node1, node2=node2, number=self.n_beams)
        self.beams[line.number] = line
        self.n_beams += 1

        return line


class ModelData(object):
    """Object to store the data of a model object. It is used to pass it to the results object"""

    def __init__(self, model):
        """Initializes the ModelData instance

        :model: a Model instance

        """
        from copy import copy

        self._name = model._name
        self._dimensionality = model._dimensionality
        self.nodes = copy(model.nodes)
        self.beams = copy(model.beams)
        self.beam_sections = copy(model.beam_sections)
        self.materials = copy(model.materials)
        self.n_nodes = model.n_nodes
        self.n_beams = model.n_beams
        self.n_dimensions = model.n_dimensions
        self.n_materials = model.n_materials
        self._connectivity = copy(model._connectivity)
        # Number of dof per node. Initialized in the respective models
        self.n_dof_per_node = model.n_dof_per_node
        # Specify dofs that are not active due to border conditions
        self._dof_dirichlet = copy(model._dof_dirichlet)


def get_dataframe_of_node_coords(model, nodes='all'):
    """Return a pandas dataframe with coordinates of selected nodes of the model

    :nodes: list of nodes or 'all'
    :returns: TODO

    """
    dimensions = model.n_dimensions
    #
    if nodes == 'all':
        nodes = [i for i, n in model.nodes.items()]

    ar_coords = np.zeros((len(nodes), dimensions), dtype=np.float)
    index_nodes = np.zeros(len(nodes), dtype=np.int)

    for ix_node, curr_node in enumerate(nodes):
        node_i = model.nodes[curr_node]
        ar_coords[ix_node, :] = node_i.coords
        index_nodes[ix_node] = curr_node

    # Set coordinate labels according to the model
    if dimensions == 2:
        index_label = ['x', 'y']
    else:
        index_label = ['x', 'y', 'z']

    # Append to the Dta Frame
    df_coords = pd.DataFrame(
        data=ar_coords,
        index=index_nodes,
        dtype=np.float64,
        columns=index_label)

    return df_coords
