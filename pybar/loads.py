#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Defines the different load types,thta can be applied to the elements, which are
transfered accordingly to the respective nodes.
"""
import numpy as np
import scipy.sparse as sparse

class Load(object):

    """Defines the Load object."""

    def __init__(self):
        """Initialize the Load instance"""
        self._type = ''
        

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
        self._direction = direction
        self._coord_system = coord_system
        self._type = 'Distributed Load'

        # Detect if distribution is a varying distributed load or not
        if p2 == None:
            self.is_uniform = True
            p2 = p1
            self._p2 = p1
        else:
            self.is_uniform = False
            self._p2 = p2

        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        load_v = np.zeros(6)
        # Calculate the transfer matrix for the axial load
        # (direction='x')
        if direction == 'x':
            load_v[0] = elem._length * (2.*p1 + p2) / 6.
            load_v[3] = elem._length * (p1 + 2.*p2) / 6.

        elif direction == 'z':
            load_v[1] = elem._length * (7.*p2 + 3.*p1) / 20.
            load_v[2] = elem._length**2 * (p2/20. + p1/30.)
            load_v[4] = elem._length * (3.*p2 + 7*p1) / 20.
            load_v[5] = -elem._length**2 * (p2/30. + p1/20.)

        self._loading_vector = load_v
        # Calculate the load vector in global coordinates, using the
        # transformation matrix
        T = elem.transformation_matrix
        self._load_vector_global = T.T.dot(load_v)

class DistributedMoment(Load):

    """Docstring for DistributedMoment. """

    def __init__(self, elem, m1, m2=None, direction='z', coord_system='local'):
        """Apply a distributed moment to a beam element

        :elem: TODO
        :m1: TODO
        :m2: TODO

        """
        Load.__init__(self)

        self._elem = elem
        self._m1 = m1
        self._m2 = m2
        self._direction = direction
        self._coord_system = coord_system
        self._type = 'Distributed Moment'
        self.is_uniform = True
        
        # Detect if distribution is a varying distributed load or not
        if m2 == None:
            self.is_uniform = True
            m2 = m1
            self._m2 = m1
        else:
            self.is_uniform = False
            self._m2 = m2

        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        load_v = np.zeros(6)

        # Calculate the transfer matrix for the axial load
        # (direction='x')
        if direction == 'z':
            load_v[1] = (m1 + m2) * 0.5
            load_v[2] = elem._length * (m1 - m2) / 12.
            load_v[4] = -(m1 + m2) * 0.5
            load_v[5] = -elem._length * (m1 - m2) / 12.

        else:
            # TODO: implement this in 3D
            pass

        self._transfer_matrix = load_v
        # Calculate the load vector in global coordinates, using the
        # transformation matrix
        T = elem.transformation_matrix
        self._load_vector_global = T.T.dot(load_v)

