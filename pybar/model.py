#!/usr/bin/env python
# encoding: utf-8
import numpy as np
from .utils import Local_Csys_two_points
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
        self.segments = dict()
        self.beam_sections = dict()
        self.materials = dict()
        self.n_nodes = 0
        self.n_segments = 0

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
        return str('Model: Name: {name}, Nodes: {n_nodes}, Segments: {n_segments}'.format(
            name=self._name, n_nodes=self.n_nodes, n_segments=self.n_segments))

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return str('Model: Name: {name}, Nodes: {n_nodes}, Segments: {n_segments}'.format(
            name=self._name, n_nodes=self.n_nodes, n_segments=self.n_segments))

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
        # Connectivity matrix for the segments
        conn_lines = np.zeros((len(self.segments), 3))
        # For different element types:
        count = 0
        for num, curr_line in self.segments.items():
            conn_lines[count, 0] = curr_line.number
            conn_lines[count, 1] = curr_line._node1.number
            conn_lines[count, 2] = curr_line._node2.number
            count += 1

        return conn_lines

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

    def Segment(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Segment2D(node1=node1, node2=node2, number=self.n_segments)
        self.segments[line.number] = line
        self.n_segments += 1

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

    def Segment(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Segment3D(node1=node1, node2=node2, number=self.n_segments)
        self.segments[line.number] = line
        self.n_segments += 1

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

        return obj

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

class Segment(object):

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
        # TODO: check that the number is not in use
        self.number = number
        # Calculate the length of the element
        delta_x = node2.x - node1.x
        delta_y = node2.y - node1.y
        delta_z = node2.z - node1.z
        self._length = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
        # calculate angle that the line forms with the horizontal
        self._alpha = np.arctan2(delta_y, delta_x)
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
        #T = self._localCSys.calc_tranformation_matrix(self._length, cx, cy, cz)
        #self.tranformation_matrix = T

    def assign_section(self, beam_section):
        """Assign a beam section instance to the segment

        :beam_section: a BeamSection instance
        :returns: self

        """
        self._beam_section = beam_section

        return self

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Segment {number}: ({n1},{n2})'.format(number=self.number, n1=self._node1, n2=self._node2)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Segment {number}'.format(number=self.number)

class Segment2D(Segment):

    """Line objects, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        Segment.__init__(self, node1, node2, number)
        # Define status of each degree of freedom, for each node
        # 1: free; 0: restrained
        # - Node 1
        self.dof1 = 1 # trans x
        self.dof2 = 1 # trans y
        self.dof6 = 1 # rot z

        # Node 2
        self.dof7 = 1  # trans x
        self.dof8 = 1  # trans y
        self.dof12 = 1 # rot z

        # Release rotation on the ends of the segment
        self.release_node_i = False # first node
        self.release_node_j = False # second node

        # Transformation matrix
        #T = self._localCSys.calc_tranformation_matrix(self._length, cx, cy, cz)
        #self.tranformation_matrix = T

class Segment3D(Segment):

    """Line objects, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        Segment.__init__(self, node1, node2, number)
        # Define status of each degree of freedom, for each node
        # 1: free; 0: restrained
        # - Node 1
        self.dof1 = 1 # trans x
        self.dof2 = 1 # trans y
        self.dof3 = 1 # trans z
        self.dof4 = 1 # rot x
        self.dof5 = 1 # rot y
        self.dof6 = 1 # rot z

        # Node 2
        self.dof7 = 1  # trans x
        self.dof8 = 1  # trans y
        self.dof9 = 1  # trans z
        self.dof10 = 1 # rot x
        self.dof11 = 1 # rot y
        self.dof12 = 1 # rot z

        # Transformation matrix
        #T = self._localCSys.calc_tranformation_matrix(self._length, cx, cy, cz)
        #self.tranformation_matrix = T

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
        self._area = self.calc_area()
        self._I33, self._I22 = self.calc_inertia()

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
        :returns: TODO

        """
        type = self._type

        if type == 'rectangular':
            width = self._data[0]
            height = self._data[1]
            I_33 = width * height**3 / 12.
            I_22 = height * width**3 / 12.
            return I_33, I_22

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


