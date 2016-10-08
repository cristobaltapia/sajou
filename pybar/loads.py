#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Defines the different load types,thta can be applied to the elements, which are
transfered accordingly to the respective nodes.
"""
import numpy as np

class Load(object):

    """Defines the Load object."""

    def __init__(self):
        """Initialize the Load instance"""
        

class DistributedLoad(Load):

    """Docstring for DistributedLoad. """

    def __init__(self, elem, p1, p2=None, direction='z', coord_system='local'):
        """Apply a distributed load on the frame element.

        If 'p2' is given, then a linearly varying load is applied with value 'p1' at the
        fisrt node of the element and 'p2' at the second. Otherwise the load is uniformly
        distributed with value 'p1'.

        The direction can be established, being possible to use either 'z' (default) or 'x'. Also
        the coordinate system can be change between 'local' (default) or 'global'.

        :elem: beam element
        :p1: TODO
        :p2: TODO
        :direction: TODO
        :coord_system: TODO

        """
        Load.__init__(self)

        self._elem = elem
        self._p1 = p1
        self._p2 = p2
        self._direction = direction
        self._coord_system = coord_system

        # Detect if is a uniformly distributed load or not
        if p2 == None:
            is_uniform = True

        # Initialize transfer matrix
        # FIXME: make this dependant from the specific element.
        tr = np.zeros(6)
        # Calculate the transfer matrix for the axial load
        # (direction='x')
        if direction == 'x':
            if is_uniform:
                tr[0] = p1 * elem._length * 0.5
                tr[3] = p1 * elem._length * 0.5
            else:
                tr[0] = elem._length * (2.*p1 + p2) / 6.
                tr[3] = elem._length * (p1 + 2.*p2) / 6.

        elif direction == 'z':
            if is_uniform:
                tr[1] = -p1 * elem._length * 0.5
                tr[2] = p1 * elem._length / 12.
                tr[4] = -p1 * elem._length * 0.5
                tr[5] = -p1 * elem._length / 12.
            else:
                tr[1] = -elem._length * (7.*p1 + 3*p2) / 20.
                tr[2] = elem._length**2 * (p1/20. + p2/30.)
                tr[4] = -elem._length * (3.*p1 + 7*p2) / 20.
                tr[5] = -elem._length**2 * (p1/30. + p2/20.)

        self._transfer_matrix = tr
        # Calculate the transfer matrix in global coordinates, using the
        # transformation matrix
        T = elem.transformation_matrix
        self._transfer_matrix_global = np.dot(T.T, tr)

