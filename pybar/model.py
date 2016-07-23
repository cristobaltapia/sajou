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
        self.materials = dict()
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

    def Segment(self, node1, node2):
        """Define a line between two nodes.

        :node1: first node
        :node2: second node

        """
        line = Segment(node1=node1, node2=node2, number=self.n_segments)
        self.segments[line.number] = line
        self.n_segments += 1

        return line

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
        conn_lines = np.zeros((len(self.segments), 2))
        # For different element types:
        count = 0
        for num, curr_line in self.segments.items():
            conn_lines[count, 0] = curr_line._node1.number
            conn_lines[count, 1] = curr_line._node2.number
            count += 1

        return conn_lines


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

class Segment(object):

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
        return 'Segment {number}: ({n1},{n2})'.format(number=self.number, n1=self._node1, n2=self._node2)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Segment {number}'.format(number=self.number)

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

