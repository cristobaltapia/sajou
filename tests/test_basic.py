# -*- coding: utf-8 -*-

import numpy as np
from .context import pybar
model = pybar.model

import unittest


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_model_creation(self):
        # create model
        test_model = pybar.model.Model(name='Test Model', dimensionality="2D")
        self.assertEqual(test_model._name, 'Test Model')

    def test_node_generation(self):
        """Test the generation of nodes
        """
        test_model = pybar.model.Model(name='Test Model', dimensionality="2D")
        node_1 = test_model.Node((14.6, 17.8))

        self.assertEqual(node_1.x, 14.6)
        self.assertEqual(node_1.y, 17.8)

    def test_multiple_node_generation(self):
        """Test the generation of multiple nodes. Important is the correct definition of the id-number and coordinates

        """
        test_model = pybar.model.Model(name='Test Model', dimensionality="2D")
        nodes_id = np.arange(10)
        print("\n+--------------------------------------------------+")
        print("| Creating nodes:")
        for curr_node in nodes_id:
            x = 1*np.cos(2*np.pi*curr_node / 10.)
            y = 1*np.sin(2*np.pi*curr_node / 10.)
            new_node = test_model.Node((x, y))
            print("| --> " + str(new_node))
            self.assertEqual(new_node.number, curr_node)

        print("| Compare to the list of nodes of the model")
        list_model = np.array([ix for ix in test_model.nodes])
        for curr_node in nodes_id:
            self.assertEqual(list_model[curr_node], curr_node)
        print(test_model.nodes)
        print("+--------------------------------------------------+")

    def test_line_generation(self):
        """Test the generation of line elements between two nodes

        """
        test_model = pybar.model.Model(name='Test Model', dimensionality="2D")
        node_1 = test_model.Node((0,0))
        node_2 = test_model.Node((10,20.5))
        # Line element
        line = test_model.Line(node1=node_1, node2=node_2)

        self.assertEqual(line.number, 0)
        self.assertEqual(line._node1, node_1)
        self.assertEqual(line._node2, node_2)

    def test_multiple_line_generation(self):
        """Test the generation of multiple lines

        """
        test_model = pybar.model.Model(name='Test Model', dimensionality="2D")
        # create Node at (0,0)
        node_0 = test_model.Node((0,0))
        #
        nodes_id = np.arange(10)
        # Create lines between this node and the others
        print("\n+--------------------------------------------------+")
        print("| Creating multiple lines:")
        for curr_node in nodes_id:
            x = 1*np.cos(2*np.pi*curr_node / 10.)
            y = 1*np.sin(2*np.pi*curr_node / 10.)
            # create node
            new_node = test_model.Node((x, y))
            # create line
            new_line = test_model.Line(node1=node_0, node2=new_node)
            print('| --> ' + str(new_line))
            self.assertEqual(new_line.number, curr_node)

        # Print all lines
        print("| Segments in the model:")
        print(test_model.segments)
        print("+--------------------------------------------------+")



if __name__ == '__main__':
    unittest.main()
