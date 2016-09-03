#!/usr/bin/env python
# encoding: utf-8
import numpy as np
from .utils import Local_Csys_two_points
from .stiffness import assemble_Ke_2D
"""
Includes all the necessary objects to create the model.
"""

class Model(object):
    """Defines a model object"""

    def __init__(self, name, dimensionality):
        """TODO: to be defined1.

        :name: TODO
        :dimensionality: TODO

        """
        self._name = name
        self._dimensionality = dimensionality
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
        # Number of dof per node. Initialized in the respective models
        self.n_dof_per_node = None

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
        conn_matrix = np.zeros((len(self.beams), 3))
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
        connectivity = np.zeros([n_v, num_dof])

        # Assemble the connectivity matrix
        for n_elem, elem in self.beams.items():
            n_nodes = elem._n_nodes
            for ix, node in enumerate(elem._nodes):
                n_node = node.number
                i1 = n_elem*n_dof*n_nodes + ix * n_dof
                i2 = n_elem*n_dof*n_nodes + ix * n_dof + n_dof
                j1 = n_node * n_dof
                j2 = n_node * n_dof + n_dof
                connectivity[i1:i2,j1:j2] = np.eye(n_dof)

        self._connectivity = connectivity

        return connectivity

    def _assemble_global_K(self):
        """Assembles the global stiffness matrix using the addition method
        :returns: numpy array

        """
        # FIXME: generalize the assembly of the stiffness matrix for
        # other elements, too
        connect_table = self._generate_connectivity_table2D()
        # Prealocate the matrix
        num_dof = self.n_nodes * self.n_dof_per_node
        n_dof = self.n_dof_per_node
        #
        K = np.zeros([num_dof, num_dof])
        # Fill the global stiffness matrix
        for n_elem in self.beams:
            n1 = connect_table[n_elem, 1]
            n2 = connect_table[n_elem, 2]
            self.beams[n_elem].assembleK()
            #
            i1 = int(n_dof*n1)
            i2 = int(n_dof*(n1 + 1))
            K[i1:i2,i1:i2] = self.beams[n_elem]._Ke[0:n_dof,0:n_dof]
            j1 = int(n_dof*n2)
            j2 = int(n_dof*(n2 + 1))
            K[j1:j2,j1:j2] = self.beams[n_elem]._Ke[n_dof:,n_dof:]

        self._K = K

        return K

    def _generate_load_matrix(self):
        """This function generates the global load matrix P
        :returns: numpy array

        """
        # Initialize a zero vector of the size of the total number of
        # dof
        P = np.zeros(self.n_nodes*self.n_dof_per_node)
        # Assign the values corresponding to the loads in each dof
        for ix, node_i in self.nodes.items():
            for dof, val in node_i._Loads.items():
                ind = ix*self.n_dof_per_node + dof
                P[ind] = val

        self._P = P

        return P

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

class Model2D(Model):
    """Subclass of the 'Model' class. It is intended to be used for the 2-dimensional
    models of frame structures."""

    def __init__(self, name):
        """TODO: to be defined1. """
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)
        self.n_dof_per_node = 3 # dof per node

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

class Node(np.ndarray):
    """3-dimensional implementation of Nodes"""
    def __new__(cls, x, y, z, number):
        # A z-coordinate is incorporated if only two are given
        obj = np.asarray([x, y, z]).view(cls)
        obj.number = number
        obj.x = x
        obj.y = y
        obj.z = z
        # number of degree of freedom of the node
        obj.n_dof = None

        # Border conditions on the DOF:
        # The border conditions are defined in this dictionary.
        # The keys correspond to the dof restrained and the value to the
        # border condition applied. The convention used ist that each
        # key value is associated with the corresponding DOF, i.e. in a
        # 2D model the key value '0' is associated with a translation in
        # the first (x) direction, while a key value '2' is associated
        # with a rotation around the third (z)-axis.
        # Only global coordinates are allowed here. For use local
        # coordinate system or another coordinate system the BC have
        # to be transformed accordingly first.
        # The function set_BC() is used to do this.
        #
        obj._BC = dict()
        # Forces (moments) on the node
        # The forces on the node are defined analogously as for the
        # application of border conditions. The function used to apply
        # the loads is set_Load().
        #
        obj._Loads = dict()

        # dictionary containing the beams that use this node
        obj.beams = dict()

        return obj

    def set_BC(self, dof, val):
        """Adds a BC to the specified dof of the Node.

        :dof: specified degree of freedom
        :val: value given
        :returns: TODO

        """
        self._BC[dof] = val

        return 1

    def set_Load(self, dof, val):
        """Adds a Load to the specified dof of the Node.

        :dof: specified degree of freedom
        :val: value given
        :returns: TODO

        """
        self._Loads[dof] = val

        return 1

    def append_beam(self, beam, node):
        """Appends the information of the beam that uses the node and the corresponding
        denomintaion: node 1 or 2

        :beam: Beam instance
        :node: '1' or '2'
        :returns: nothing FIXME

        """
        self.beams[beam.number] = node

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Node {number}'.format(number=self.number)

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Node {number}: ({x},{y},{z})'.format(number=self.number, x=self.x, y=self.y, z=self.z)

class Node2D(Node):
    """2-dimensional implementation of Nodes"""
    def __init__(self, x, y, z, number):
        """ Initializes a Node2D instance
        """
        #Node.__init__(self, x=x, y=y, z=0.0, number=number)
        # Numbe od DOF per node:
        # Translation in x, translation in y, rotation in z
        self.n_dof = 3

        return None

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Node2D {number}'.format(number=self.number)

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Node2D {number}: ({x},{y},{z})'.format(number=self.number, x=self.x, y=self.y, z=self.z)

class Beam(object):
    """Line objects, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        self._beam_section = None
        # TODO: accept tuples with coordinates also
        self._node1 = node1
        self._node2 = node2
        self._nodes = [self._node1, self._node2]
        # TODO: check that the number is not in use
        self.number = number
        #
        node1.append_beam(self, 1)
        node2.append_beam(self, 2)
        # Calculate the length of the element
        self._length = np.linalg.norm(node2-node1)
        # Local coordinate system
        self._localCSys = Local_Csys_two_points(point1=node1, point2=node2, type='cartesian')
        # Section
        self._beam_section = None
        # Directive cosines
        delta = node2 - node1
        cx = delta[0] / self._length
        cy = delta[1] / self._length
        cz = delta[2] / self._length

        # Transformation matrix
        self.transformation_matrix = self._localCSys.calc_tranformation_matrix(self._length, cx, cy, cz)

        self._Ke = None

    def assign_section(self, beam_section):
        """Assign a beam section instance to the beam

        :beam_section: a BeamSection instance
        :returns: self

        """
        self._beam_section = beam_section

        return self

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam {number}: (N{n1}, N{n2})'.format(number=self.number,
                n1=self._node1.number, n2=self._node2.number)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam {number}'.format(number=self.number)

class Beam2D(Beam):
    """Line objects, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        Beam.__init__(self, node1, node2, number)
        # Number of nodes of the element
        self._n_nodes = 2
        # Number of DOF per node
        self._dof_per_node = 3
        # displacement/rotation of each degree of freedom, for each node
        # - Node 1
        self.dof1 = 0. # trans x
        self.dof2 = 0. # trans y
        self.dof3 = 0. # rot z

        # Node 2
        self.dof4 = 0. # trans x
        self.dof5 = 0. # trans y
        self.dof6 = 0. # rot z

        # Release rotation on the ends of the beam
        self.release_node_1 = False # first node
        self.release_node_2 = False # second node

    def assembleK(self):
        """Assembles the stiffnes matrix for the element
        :returns: stiffness matrix in global coordinates

        """
        Ke = assemble_Ke_2D(self)

        self._Ke = Ke

        return Ke

class Beamt3D(Beam):
    """Line objects, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        Beam.__init__(self, node1, node2, number)
        # displacement/rotation of each degree of freedom, for each node
        # - Node 1
        self.dof1 = 0. # trans x
        self.dof2 = 0. # trans y
        self.dof3 = 0. # trans z
        self.dof4 = 0. # rot x
        self.dof5 = 0. # rot y
        self.dof6 = 0. # rot z

        # Node 2
        self.dof7 = 0.  # trans x
        self.dof8 = 0.  # trans y
        self.dof9 = 0.  # trans z
        self.dof10 = 0. # rot x
        self.dof11 = 0. # rot y
        self.dof12 = 0. # rot z

class Material(object):
    """Material properties"""

    # TODO: create function to print information of the material
    def __init__(self, name, data, type='isotropic'):
        """TODO: to be defined1.

        :name: name of the material
        :table: properties of the material
        :type: type of the material:
            - 'isotropic': data = (E, )
            - 'orthotropic': data = (E_1, E_2, E_3)

        """
        self._name = name
        self._data = data
        self._type = type

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Material: {name}'.format(name=self._name)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Material: {name}'.format(name=self._name)

class BeamSection(object):
    """Defines a beam section"""

    def __init__(self, name, material, data, type='rectangular'):
        """Initiates the section

        :name: name of the section
        :material: material of the section defined as an instance of Material object
        :data: properties of the section:
        :type: defines the type of cross-section
                - rectangular: data=(width, height,)
                - circular:    data=(r, )
                - I-section:  data=(H, h_f, w_web, w_f)
                - general:    data=(A, I_3,)

        """
        self._name = name
        self._material = material
        self._data = data
        self._type = type
        self._area = 0
        self._Iz = 0
        self._Iy = 0
        self._Jx = 0
        self.compute_properties()

    def print_properties(self):
        """Prints the properties of the BeamSection instance
        :returns: TODO

        """
        if self._type == 'rectangular':
            props = {'width':self._data[0], 'height':self._data[1]}
        else:
            props = 'undefined'

        return 'Properties: ' + str(props)

    def compute_properties(self):
        """Compute all the mechanical properties for the given section
        :returns: TODO

        """
        # Calculate the area
        self._area = self.calc_area()
        self._Iz, self._Iy = self.calc_inertia()

    def calc_area(self):
        """Calculate the area of the section
        :returns: TODO

        """
        type = self._type

        if type == 'rectangular':
            width = self._data[0]
            height = self._data[1]
            return width * height
        elif type == 'general':
            return self._data[0]

        elif type == 'circular':
            radius = self._adta[0]
            return np.pi * radius**2

    def calc_inertia(self):
        """Calculate the moment of inertia of the beam section
        :returns: Iz, Iy

        """
        type = self._type

        if type == 'rectangular':
            width = self._data[0]
            height = self._data[1]
            I_z = width * height**3 / 12.
            I_y = height * width**3 / 12.
            return I_z, I_y
        elif type == 'general':
            return self._data[1]

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam Section: {name}, type: {t}'.format(name=self._name, t=self._type)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam Section: {name}, type: {t}'.format(name=self._name, t=self._type)


