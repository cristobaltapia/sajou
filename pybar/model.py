#!/usr/bin/env python
# encoding: utf-8
import numpy as np
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
        self.n_nodes = 0
        self.n_segments = 0

    def Node(self, position):
        """TODO: Docstring for Node.

        :*kwargs: TODO
        :returns: TODO

        """
        node = Node(position=position, number=self.n_nodes)
        self.nodes[node.number] = node
        self.n_nodes += 1

        return node

    def Line(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Line(node1=node1, node2=node2, number=self.n_segments)
        self.segments[line.number] = line
        self.n_segments += 1

        return line

    def __str__(self):
        """
        Printable string
        """
        str('Model: {name}, Nodes: {n_nodes}, Segments: {n_segments}'.format(
            name=self._name, n_nodes=self.n_nodes, n_segments=self.n_segments))

    def _ensemble_K(self):
        """Generates the global stiffness matrix

        :returns: TODO

        """
        pass

    def generate_connectivity_matrix(self):
        """Generates the connectivity matrix for the model
        :returns: TODO

        """
        C_lines = np.array()
        # For different element types:
        for curr_line in self.segments:
            pass


class Node(object):

    """2-Dimensional nodes"""

    def __init__(self, position, number):
        """
        position: tuple (x,y)
        """
        self.x = position[0]
        self.y = position[1]
        self.number = number

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Node {number}'.format(number=self.number)

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Node {number}: ({x},{y})'.format(number=self.number, x=self.x, y=self.y)

class Line(object):

    """2-Dimensional line, joining two nodes"""

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        # TODO: accept tuples with coordinates also
        self._node1 = node1
        self._node2 = node2
        # TODO: check that the number is not in use
        self.number = number
        # calculate angle that the line forms with the horizontal
        delta_x = node2.x - node1.x
        delta_y = node2.y - node1.y
        self._alpha = np.arctan2(delta_y, delta_x)

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Line {number}: ({n1},{n2})'.format(number=self.number, n1=self._node1, n2=self._node2)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Line {number}'.format(number=self.number)

class Material(object):

    """Material properties"""

    def __init__(self, table, type='isotropic'):
        """TODO: to be defined1.

        :table: TODO
        :type: TODO

        """
        self._table = table
        self._type = type


class Section(object):

    """Defines a section"""

    def __init__(self, material, table, type='rectangular'):
        """Initiates the section

        :material: material object
        :table: TODO
        :type: defines the type of cross-section
                - rectangular: table=(width, height,)
                - circular:    table=(r, )
                - triangular:  table=(base, height)

        """
        self._material = material
        self._table = table
        self._type = type

