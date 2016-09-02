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

    def _generate_connectivity_matrix2D(self):
        """Generates the connectivity matrix for the model
        :returns: TODO


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

    def _generate_connectivity_matrix2D_a(self):
        """Generates the connectivity matrix 'a' (also known as incidence matrix or locator) such that,

            v = a*V

            where "v" correspond to the vector containing displacements at the same node
            from different elements
            
            "V" contains the global displacements of all DOF's
            "a" is the incidence matrix

        :returns: a numpy array

        """
        # Get the total DOF per node of each element and add it
        n_v = 0
        # FIXME: generalize for different types of elements
        for keys, elem in self.beams.items():
            n_v += elem._dof_per_node * elem._n_nodes

        # Total number of 'system displacements'
        n_v_sys = self.n_nodes * self.beams[0]._dof_per_node
        
        # create a zero matrix with the adequate size
        connectivity = np.zeros([n_v, n_v_sys])
        # Assemble the connectivity matrix
        for n_elem, elem in self.beams.items():
            n_dof = elem._dof_per_node
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



class Model2D(Model):
    """Subclass of the 'Model' class. It is intended to be used for the 2-dimensional
    models of frame structures."""

    def __init__(self, name):
        """TODO: to be defined1. """
        dimensionality = '2D'
        Model.__init__(self, name, dimensionality)

    def Node(self, x, y):
        """2D implementation of the Node.

        :*kwargs: TODO
        :returns: instance of Node

        """
        # A coordinate z=0 is passed to initiate the Node Instance
        node = Node(x=x, y=y, z=0., number=self.n_nodes)
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

        # dictionary containing the beams that use this node
        obj.beams = dict()

        return obj

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
        self.dof6 = 0. # rot z

        # Node 2
        self.dof7 = 0.  # trans x
        self.dof8 = 0.  # trans y
        self.dof12 = 0. # rot z

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


