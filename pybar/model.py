#!/usr/bin/env python
# encoding: utf-8
"""
Defines the model classes for 2D and 3D models.
"""
__docformat__ = 'reStructuredText'

import numpy as np
import pandas as pd
import scipy.sparse as sparse

from pybar.materials import Material
from pybar.nodes import Node2D
from pybar.sections import BeamSection


class Model(object):
    """Defines a model object

    Parameters
    ----------
    name: str
        name of the model
    dimensionality: str
        spacial dimensions used in the model ('2D' or '3D')

    Attributes
    ----------
    nfat: dict
        Node Freedom Allocation Table
    nodes: dict
        dictionary with all the nodes of the system
    beams: dict
        dictionary with all the beams of the system
    beam_sections: dict
        dictionary with the beam sections defined in the model
    materials: dict
        dictionary with the materials defined in the model
    n_nodes: int
        number of nodes of the system
    n_beams: int
        number of beams in the system
    n_materials: int
        number of materials defined
    n_dimensions: int
        number of spacial dimensions of the model
    n_dof_per_node: int
        number of degrees of freedom per node
    _name: str
        name of the model
    _dimensionality: str
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

    def __new__(cls, name, dimensionality):
        if cls is Model:
            if dimensionality == '2D':
                return super(Model, cls).__new__(Model2D)
            if dimensionality == '3D':
                return super(Model, cls).__new__(Model3D)
        else:
            return super(Model, cls).__new__(cls, name, dimensionality)

    def __init__(self, name, dimensionality):
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

    def material(self, name, data, type='isotropic'):
        """Function used to create a Material instance in the model

        Parameters
        ----------
        name: str
            name of the material
        data:
            data for the material
        type: str
            type of the material

        Returns
        -------
        pybar.Material
            a Material instance

        """
        material = Material(name=name, data=data, type=type)
        # Add the material to the dictionary of materials in the current
        # model
        self.materials[name] = material
        self.n_materials += 1

        return material

    def beam_section(self, name, material, data, type='rectangular'):
        """Function use to create a BeamSection instance in the model

        Parameters
        ----------
        name: str
            name of the section
        material: pybar.Material
            material for the section
        data:
            data (see BeamSection class definition)
        type:
            type of the section (see BeamSection class definition)
        Returns
        -------
        returns
            a beam section instance

        """
        # The material can be passed both as a string, corresponding to
        # a key of the material dictionary of the model, or as a
        # material instance directly.

        if isinstance(material, str):
            material_section = self.materials[material]
        else:
            material_section = material

        section = BeamSection(name=name, material=material_section, data=data,
                              type=type)
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

    def _assemble_global_K(self):
        """Assemble the global stiffness matrix using the addition method.

        Returns
        -------
        numpy.array
           Global stiffness matrix of the system

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

        Returns
        -------
        numpy.array
            Loading vector of the system

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
        """Generate the global displacement matrix V, containing the border
        conditions applied to each dof.

        Returns
        -------
        numpy.array
            Dsiplacement vector of the system

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
            index_i = np.array([kx for kx in node._bc.keys()],
                               dtype=np.int) + nfmt[node.number]
            V[index_i] = np.array([kx for kx in node._bc.values()])
            # Add to the list of restrained DOFs
            self._dof_dirichlet.extend(index_i.tolist())

        self._V = V

        return V

    def bc(self, node, type='displacement', coord_system='global', **kwargs):
        """Introduces a border condition to the node.

        Parameters
        ----------
        node: pybar.Node
            Node to which the border condition will be applied
        type: str
            type of border condition

            - Options:
                ``'displacement'``, ``...``

        coord_system:
            spcifies the coordinate system to be used when applying the BC
        **kwargs:
            keyword arguments. At least one of the following parameters must
            be supplied

        Keyword Arguments
        -----------------
        v1: float
            displacement in the direction 1
        v2: float
            displacement in the direction 2
        v3: float
            displacement in the direction 3
        r1: float
            rotation in the direction 1
        r2: float
            rotation in the direction 2
        r3: float
            rotation in the direction 3

        Returns
        -------
        bool
            True if successful

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

        return True

    # TODO: there has to give a 'Load' class to handle the different
    # type of loads.
    def load(self, node, coord_system='global', **kwargs):
        """Introduces a Load in the given direction according to the selected
        coordinate system at the specified node.

        Parameters
        ----------
        node: pybar.Node
            a Node instance
        coordinate:
            coordinate system
        **kwargs:
            keyword arguments. The BC is defined for the different degree of
            freedom (*dof*) available to the node.
            At least one of the following parameters must be supplied:

        Keyword Arguments
        -----------------
        f1: float
            force in direction 1
        f2: float
            force in direction 2
        f3: float
            force in direction 3
        m1: float
            moment in direction 1
        m2: float
            moment in direction 2
        m3: float
            moment in direction 3

        Returns
        -------
        pybar.Load
            the instance of the Load object created

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

        Returns
        -------
        pybar.model.ModelData
            the data of the whole analyzed model

        """
        model_data = ModelData(self)

        return model_data

    def get_node_and_dof(self, dof):
        """Return the node and element dof (number of the dof in a specific element)
        corresponding to the global dof given.

        Parameters
        ----------
        dof: int
            global *dof*

        Returns
        -------
        pybar.Node
            Node correpsonding to the *dof* specified
        int
            number of the *dof* of the element corresponding to the global
            *dof* supplied


        """
        # Get the node
        n_node = dof // self.n_dof_per_node
        node = self.nodes[n_node]
        # Get the intra-element dof
        n_dof = dof - (n_node * self.n_dof_per_node)

        return node, n_dof

    def add_hinge(self, node):
        """Add hinge to the specified node. Also supports list of nodes

        Parameters
        ----------
        node: pybar.Node
            Node instance or list of node instances

        Returns
        -------
        bool
            TODO

        Todo
        ----
        This function still needs work
        """
        #FIXME: not yet implemented!
        if isinstance(node, list):
            for node_i in node:
                node_i.add_hinge()
        else:
            node.add_hinge()

        return True

    def __generate_node_freedom_allocation_dict__(self):
        """
        Generate the Node Freedom Allocation Table.

        Generates a dictionary of arrays, containing the Node Freedom
        Signature of each Node of the model. The key values correspond to the node
        numbers.
        Also counts the total number of active DOFs of the system and stores it on the
        variable self.n_active_dof

        Returns
        -------
        dict
            the Node Freedom Allocation Table of the system

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

        Returns
        -------
        numpy.array
            TODO

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

    Parameters
    ----------
    name: str
        name of the model

    Attributes
    ----------
    n_dimensions: int
        number of spacial dimensions (2 for Model2D)
    n_dof_per_node: int
        number of degrees of freedom per node


    """

    def __init__(self, name, dimensionality='2D'):
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)
        # Numer of dimensions
        self.n_dimensions = 2
        # Number of degrees of freedom per node:
        self.n_dof_per_node = 3

    def node(self, x, y):
        """2D implementation of the Node.

        Parameters
        ----------
        x: float
            x position
        y: float
            y position

        Returns
        -------
        pybar.Node
            the node created

        """
        # A coordinate z=0 is passed to initiate the Node Instance
        node = Node2D(x=x, y=y, z=0.0, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def beam(self, node1, node2):
        """Define a line between two nodes.

        Parameters
        ----------
        node1: pybar.Node
            first node
        node2: pybar.Node
            second node

        Returns
        -------
        pybar.Beam
            the beam element created

        """
        from .elements.beam2d import Beam2D
        line = Beam2D(node1=node1, node2=node2, number=self.n_beams)
        self.beams[line.number] = line
        self.n_beams += 1

        return line

    def distributed_load(self, elements, **kwargs):
        """Add a distributed load to a list of beam elements.

        A list of elements has to be supplied for the first variable. The rest of the
        variables are exactly the same as in the 'distributed_load' function of the
        corresponding elements.

        Parameters
        ----------
        elements: list
            list of beam elements
        p1: float
            value of the force at start node
        p2: float
            value of the force at end node
        direction: int
            direction of the applied load (default: *2*)
        coord_system: str
            coordinate system (default: global)

        Returns
        -------
        bool
            TODO

        """
        for curr_elem in elements:
            # Add distributed load
            curr_elem.distributed_load(**kwargs)

        return True


class Model3D(Model):
    """
    Subclass of the 'Model' class. It is intended to be used for the 3-dimensional
    models of frame structures.

    Allocation of DOFs in each node:

            [1 2 3 4 5 6] = [ux, uy, uz, rx, ry, rz]

    """

    def __init__(self, name, dimensionality='3D'):
        dimensionality = '3D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 6  # dof per node

    def node(self, x, y, z):
        """
        3D implementation of the Node.

        Parameters
        ----------

        x: float
            x-position of the node
        y: float
            y-position of the node
        z: float
            z-position of the node

        Returns
        -------

        Node: instance of Node

        """
        node = Node(x=x, y=y, z=z, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def beam(self, node1, node2):
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
    df_coords = pd.DataFrame(data=ar_coords, index=index_nodes,
                             dtype=np.float64, columns=index_label)

    return df_coords
