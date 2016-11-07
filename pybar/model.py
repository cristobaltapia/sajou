#!/usr/bin/env python
# encoding: utf-8
"""
Defines the model classes for 2D and 3D models.
"""
import numpy as np
import scipy.sparse as sparse
import pandas as pd
from .utils import Local_Csys_two_points
from .solvers import StaticSolver
from .elements import Beam2D, Beam3D
from .nodes import Node2D
from .materials import Material
from .sections import BeamSection

class Model(object):
    """Defines a model object"""

    def __new__(cls, name, dimensionality):
        if cls is Model:
            if dimensionality == '2D': return super(Model, cls).__new__(Model2D)
            if dimensionality == '3D': return super(Model, cls).__new__(Model3D)
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
        self._K = None # global stiffness matrix
        self._P = None # load matrix
        self._V = None # global displacement matrix
        # Number of dimensions of the model
        self.n_dimensions = None
        # Number of dof per node. Initialized in the respective models
        self.n_dof_per_node = None
        # Specify dofs that are not active due to border conditions
        self._dof_dirichlet = []
        # Node Freedom Allocation Table
        self._nfat = dict()

    def Material(self, name, data, type='isotropic'):
        """Function used to create a Material instance in the model

        :name: name of the material
        :data: data for the material
        :type: type of the material
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

        :name: name of the section
        :material: material for the section
        :data: data (see BeamSection class definition)
        :type: type of the section (see BeamSection class definition)
        :returns: a beam section instance

        """
        # The material can be passed both as a string, corresponding to
        # a key of the material dictionary of the model, or as a
        # material instance directly.

        if isinstance(material, str):
            material_section = self.materials[material]
        else:
            material_section = material

        section = BeamSection(name=name, material=material_section, data=data, type=type)
        # Add section to the list of beam sections
        self.beam_sections[name] = section

        return section

    def __str__(self):
        """
        Printable string
        """
        return str('Model: Name: {name}, Nodes: {n_nodes}, Beams: {n_beams}'.format(
            name=self._name, n_nodes=self.n_nodes, n_beams=self.n_beams))

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return str('Model: Name: {name}, Nodes: {n_nodes}, Beams: {n_beams}'.format(
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

    def _generate_connectivity_matrix2D(self):
        """Generates the connectivity matrix 'a' (also known as incidence matrix or locator) such that,

            v = a*V

            where "v" correspond to the vector containing displacements at the same node
            from different elements

            "V" contains the global displacements of all DOF's
            "a" is the incidence matrix

        :returns: numpy array

        """
        # Get the total DOF per node of each element and add it
        n_v = 0
        #
        for keys, elem in self.beams.items():
            n_v += elem._dof_per_node * elem._n_nodes

        # Total number of DOFs
        num_dof = self.n_nodes * self.beams[0]._dof_per_node

        # DOF per node
        n_dof = self.n_dof_per_node

        # create a zero matrix with the adequate size
        connectivity = np.zeros([n_v, num_dof], dtype=np.float64)

        # Assemble the connectivity matrix
        for n_elem, elem in self.beams.items():
            n_nodes = elem._n_nodes
            for ix, node in enumerate(elem._nodes):
                n_node = node.number
                i1 = n_elem*n_dof*n_nodes + ix * n_dof
                i2 = n_elem*n_dof*n_nodes + ix * n_dof + n_dof
                j1 = n_node * n_dof
                j2 = n_node * n_dof + n_dof
                connectivity[i1:i2,j1:j2] = np.eye(n_dof, dtype=np.float64)

        self._connectivity = connectivity

        return connectivity

    def _assemble_global_K(self):
        """Assemble the global stiffness matrix using the addition method.
        :returns: numpy array

        """
        # FIXME: generalize the assembly of the stiffness matrix for
        # other elements, too
        connect_table = self._generate_connectivity_table2D()
        # Prealocate the matrix
        num_dof = self.n_nodes * self.n_dof_per_node
        # number of dof per node
        n_dof = self.n_dof_per_node
        #
        K = np.zeros([num_dof, num_dof], dtype=np.float64)
        # Fill the global stiffness matrix
        for n_elem in self.beams:
            n1 = connect_table[n_elem, 1]
            n2 = connect_table[n_elem, 2]
            self.beams[n_elem].assemble_K()
            #
            i1 = int(n_dof*n1)
            i2 = int(n_dof*(n1 + 1))
            j1 = int(n_dof*n2)
            j2 = int(n_dof*(n2 + 1))
            K[i1:i2,i1:i2] += self.beams[n_elem]._Ke[0:n_dof,0:n_dof]
            K[j1:j2,j1:j2] += self.beams[n_elem]._Ke[n_dof:,n_dof:]
            K[j1:j2,i1:i2] += self.beams[n_elem]._Ke[n_dof:,0:n_dof]
            K[i1:i2,j1:j2] += self.beams[n_elem]._Ke[0:n_dof,n_dof:]

        # Generate sparse matrix
        K_s = sparse.csr_matrix(K)

        return K_s

    def _generate_loading_vector(self):
        """Generate the global matrix of applied forces P.
        :returns: numpy array

        """
        # Initialize a zero vector of the size of the total number of
        # dof
        P = np.zeros(self.n_nodes*self.n_dof_per_node, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node_i in self.nodes.items():
            for dof, val in node_i._Loads.items():
                ind = ix*self.n_dof_per_node + dof
                P[ind] = val

        self._P = P

        return P

    def _generate_element_loading_vector(self):
        """Generate the global element vector of forces.
        :returns: numpy array

        """
        # Initialize a zero vector of the size of the total number of
        # DOF
        P = np.zeros(self.n_nodes*self.n_dof_per_node, dtype=np.float64)

        # Add loads applied to the elements (distributed loads)
        for ix, elem in self.beams.items():
            if len(elem._loads) > 0:
                # Get the correct indices
                # First node:
                n1 = elem._node1.number
                # DOFs coresponding to the node 1
                ind1 = n1 * self.n_dof_per_node
                # Second node:
                n2 = elem._node2.number
                # DOFs coresponding to the node 2
                ind2 = n2 * self.n_dof_per_node
                # Add to the global load vector
                P[[ind1, ind1+1, ind1+2, ind2, ind2+1, ind2+2]] += elem._load_vector_e

        self._P = P

        return P

    def _generate_displacement_vector(self):
        """This function generates the global displacement matrix V, containing the border
        conditions applied to each dof.

        :returns: numpy array

        """
        # Initialize a zero vector of the size of the total number of
        # dof
        V = np.zeros(self.n_nodes*self.n_dof_per_node, dtype=np.float64)
        # Assign the values corresponding to the loads in each dof
        for ix, node_i in self.nodes.items():
            for dof, val in node_i._BC.items():
                # Compute index corresponding to the current dof
                ind = int(ix*self.n_dof_per_node + dof)
                # Add value to the vector
                V[ind] = val
                # Add to the list of restrained DOFs
                self._dof_dirichlet.append(ind)

        self._V = V

        return V

    def BC(self, node, type='displacement', coord_system='global', **kwargs):
        """Introduces a border condition to the node.

        :node: a Node instance
        :type: type of border condition
            - Options: - 'displacement'
                       - ...
        :**kwargs: optional arguments. The BC is defined for the different degree of
            freedom available to the node.

            vi : translation of the i-th direction in the specified dof
            ri : rotation of the i-th direction in the specified dof

            For Beam2D elements:
                - v1, v2, r3
            For Beam3D elements:
                - v1, v2, v3, r1, r2, r3
            ...
        :returns: TODO

        """
        # FIXME: currently only in global coordintaes. Implement
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

    # FIXME: there has to give a 'Load' class to handle the different
    # type of loads.
    def Load(self, node, coord_system='global', **kwargs):
        """Introduces a Load in the given direction according to the selected coordinate
        system at the specified node.

        :node: a Node instance
        :**kwargs: optional arguments. The BC is defined for the different degree of
            freedom available to the node.

            fi : force on the i-th direction in the specified dof
            mi : moment on the i-th direction in the specified dof

            For Beam2D elements:
                - f1, f2, m3
            For Beam3D elements:
                - f1, f2, f3, m1, m2, m3
            ...
        :returns: TODO

        """
        # FIXME: currently only in global coordintaes. Implement
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

        :dof: global DOF (integer)
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

        :node: Node instance or list of node instances
        :returns: TODO

        """
        #FIXME: not yet implemented!
        if isinstance(node, list):
            for node_i in node:
                node_i.add_hinge()
        else:
            node.add_hinge()

    def generate_node_freedom_allocation_dict(self):
        """
        Generate the Node Freedom Allocation Table.

        Generates a dictionary of arrays, containing the Node Freedom
        Signature of each Node of the model. The key values correspond to the node
        numbers.

        :returns: a dictionary

        """
        # Loop over each node and add its Node Freedom Signature to the
        # dictionary
        for n_num, node in self.nodes.items():
            self.nfat[node.number] = node.nfs

        return self.nfat


class Model2D(Model):
    """Subclass of the 'Model' class. It is intended to be used for the 2-dimensional
    models of frame structures."""

    def __init__(self, name, dimensionality='2D'):
        """TODO: to be defined1. """
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 3 # dof per node
        self.n_dimensions = 2

    def Node(self, x, y):
        """2D implementation of the Node.

        :*kwargs: TODO
        :returns: instance of Node

        """
        # A coordinate z=0 is passed to initiate the Node Instance
        node = Node2D(x=x, y=y, z=0.0, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def Beam(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Beam2D(node1=node1, node2=node2, number=self.n_beams)
        self.beams[line.number] = line
        self.n_beams += 1

        return line

    def distributed_load(self, elements, p1, p2=None, direction='z', coord_system='local'):
        """Add a distributed load to a list of beam elements.
        A list of elements has to be supplied for the first variable. The rest of the
        variables are exactly the same as in the 'distributed_load' function of the
        corresponding elements.

        :elements: list of beams elements
        :p1: TODO
        :p2: TODO
        :returns: TODO

        """
        for curr_elem in elements:
            # Add distributed load
            curr_elem.distributed_load(p1, p2, direction, coord_system)

class Model3D(Model):
    """Subclass of the 'Model' class. It is intended to be used for the 3-dimensional
    models of frame structures."""

    def __init__(self, name, dimensionality='3D'):
        """TODO: to be defined1. """
        dimensionality = '3D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 6 # dof per node

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

class SymbolicModel2D(Model):
    """Defines a symbolic model object"""


    def __init__(self, name):
        """TODO: to be defined1. """
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 3 # dof per node

    def _assemble_global_sym_K(self):
        """Assembles the global stiffness matrix using the addition method in
        a symbolical manner
        :returns: numpy array

        """
        from sympy.matrices import zeros
        # FIXME: generalize the assembly of the stiffness matrix for
        # other elements, too
        connect_table = self._generate_connectivity_table2D()
        # Prealocate the matrix
        num_dof = self.n_nodes * self.n_dof_per_node
        # number of dof per node
        n_dof = self.n_dof_per_node
        #
        K = zeros(num_dof, num_dof)
        # Fill the global stiffness matrix
        for n_elem in self.beams:
            n1 = connect_table[n_elem, 1]
            n2 = connect_table[n_elem, 2]
            self.beams[n_elem].assemble_sym_K()
            #
            i1 = int(n_dof*n1)
            i2 = int(n_dof*(n1 + 1))
            j1 = int(n_dof*n2)
            j2 = int(n_dof*(n2 + 1))
            K[i1:i2,i1:i2] += self.beams[n_elem]._Ke[0:n_dof,0:n_dof]
            K[j1:j2,j1:j2] += self.beams[n_elem]._Ke[n_dof:,n_dof:]
            K[j1:j2,i1:i2] += self.beams[n_elem]._Ke[n_dof:,0:n_dof]
            K[i1:i2,j1:j2] += self.beams[n_elem]._Ke[0:n_dof,n_dof:]

        return K

    def _generate_loading_vector(self):
        """This function generates the global matrix of applied forces P
        :returns: numpy array

        """
        from sympy.matrices import zeros
        # Initialize a zero vector of the size of the total number of
        # dof
        P = zeros(self.n_nodes*self.n_dof_per_node, 1)
        # Assign the values corresponding to the loads in each dof
        for ix, node_i in self.nodes.items():
            for dof, val in node_i._Loads.items():
                ind = ix*self.n_dof_per_node + dof
                P[ind] = val

        self._P = P

        return P

    def _generate_displacement_vector(self):
        """This function generates the global displacement matrix V, containing the border
        conditions applied to each dof.

        :returns: numpy array

        """
        from sympy.matrices import zeros
        # Initialize a zero vector of the size of the total number of
        # dof
        V = zeros(self.n_nodes*self.n_dof_per_node, 1)
        # Assign the values corresponding to the loads in each dof
        for ix, node_i in self.nodes.items():
            for dof, val in node_i._BC.items():
                # Compute index corresponding to the current dof
                ind = int(ix*self.n_dof_per_node + dof)
                # Add value to the vector
                V[ind] = val
                # Add to the list of restrained DOFs
                self._dof_dirichlet.append(ind)

        self._V = V

        return V

    def Node(self, x, y):
        """2D implementation of the Node.

        :*kwargs: TODO
        :returns: instance of Node

        """
        # A coordinate z=0 is passed to initiate the Node Instance
        node = Node2D(x=x, y=y, z=0.0, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def Beam(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Beam2D(node1=node1, node2=node2, number=self.n_beams)
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

    ar_coords = np.zeros((len(nodes), dimensions) , dtype=np.float)
    index_nodes = np.zeros(len(nodes), dtype=np.int)

    for ix_node, curr_node in enumerate(nodes):
        node_i = model.nodes[curr_node]
        ar_coords[ix_node,:] = node_i.coords
        index_nodes[ix_node] = curr_node

    # Set coordinate labels according to the model
    if dimensions == 2:
        index_label = ['x', 'y']
    else:
        index_label = ['x', 'y', 'z']

    # Append to the Dta Frame
    df_coords = pd.DataFrame(data=ar_coords, index=index_nodes, dtype=np.float64,
            columns=index_label)

    return df_coords

