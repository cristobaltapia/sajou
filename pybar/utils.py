#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

class CoordSys(object):

    """Class implementing coordinating systems"""

    def __init__(self, type='cartesian'):
        """TODO: to be defined1.

        :type: type of the coordinate system
            - 'cartesian'
            - 'cylindrical'
            - 'spherical'
        """
        self._type = type

    def from_cart_to_cylindric(self, cart_csys):
        """Transforms coordinates from cartesian to polar coordinates

        :cart_csys: TODO
        :returns: TODO

        """
        pass

    def from_cylindric_to_cart(self, cylindric_csys):
        """Trnsforms coordinates from polar to cartesian coordinates

        :polar_csys: TODO
        :returns: TODO

        """
        pass

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'CSys {t}'.format(t=self._type)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'CSys {t}'.format(t=self._type)


class Csys_two_points(CoordSys):

    """Class implementing local coordinate systems"""

    def __init__(self, point1, point2, type='cartesian'):
        """Initiates an instance

        :type: TODO

        """
        CoordSys.__init__(self)

        self._type = type
        self._point1 = point1
        self._point2 = point2
        ############################################################
        # Calculate unit vectors
        ############################################################
        # - Convert to numpy arrays
        p1 = np.array(point1)
        p2 = np.array(point2)
        # vector in z-direction
        vz = np.array([0,0,1])
        # Vector in the principal direction (between the two points)
        v1 = p2 - p1
        # cross product to obtain v3
        v2 = np.cross(v1, vz)
        # vector in the 2nd direction
        v3 = np.cross(v1, v2)
        # create matrix
        transformation_matrix = np.vstack((v1/np.linalg.norm(v1), v2/np.linalg.norm(v2), v3/np.linalg.norm(v3))).T

        self._v1 = transformation_matrix[0,:]
        self._v2 = transformation_matrix[1,:]
        self._v3 = transformation_matrix[2,:]

