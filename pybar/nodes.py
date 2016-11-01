#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This module contains the implementation of the nodes in 2D and 3D.
"""
import numpy as np

class Node(np.ndarray):
    """3-dimensional implementation of Nodes"""
    def __new__(cls, x, y, z, number):
        # A z-coordinate is incorporated if only two are given
        obj = np.asarray([x, y, z], dtype=np.float64).view(cls)
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
        # Similar with the reactions. These are, of course, added after the
        # solution is found.
        obj.reactions = dict()

        # dictionary containing the beams that use this node
        obj.beams = dict()
        #
        obj._hinge = False

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

    def add_hinge(self):
        """Method to add a hinge at the instance of the node
        :returns: TODO

        """
        self._hinge = True

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
        self.coords = [x, y]

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
