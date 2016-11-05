#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse as sparse

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

    def calc_transformation_matrix(self, length, cx, cy, cz):
        """Calculates the transformation matrix for the current instance of local system
        :returns: TODO

        """
        # TODO: implement for 3D
        T = np.zeros([6,6], dtype=np.float64)

        T[0,0] = cx
        T[0,1] = cy
        T[1,0] = -cy
        T[1,1] = cx
        T[2,2] = 1
        T[3,3] = cx
        T[3,4] = cy
        T[4,3] = -cy
        T[4,4] = cx
        T[5,5] = 1

        # sparse form
        T_s = sparse.csr_matrix(T)

        return T_s

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

class Local_Csys_two_points(CoordSys):

    """Class implementing local coordinate systems"""

    def __init__(self, point1, point2, type='cartesian'):
        """Initiates an instance

        :type: TODO

        """
        CoordSys.__init__(self, type=type)

        self._point1 = point1
        self._point2 = point2
        ############################################################
        # Calculate unit vectors
        ############################################################
        # FIXME: no need to transform to np array
        # - Convert to numpy arrays
        p1 = np.array(point1, dtype=np.float64)
        p2 = np.array(point2, dtype=np.float64)
        # dist
        d = (p2 - p1)
        dist = np.linalg.norm(d)
        v1 = d / dist

        # vector in z-direction
        if d[2] == 0:
            vz = np.array([0, 0, -1], dtype=np.float64)
        elif d[1] == 0:
            vz = np.array([0, 1, 0], dtype=np.float64)
        elif d[0] == 0:
            vz = np.array([1, 0, 0], dtype=np.float64)
        else:
            vz = np.array([0, 0, -1], dtype=np.float64)

        # cross product to obtain v3
        v2 = np.cross(v1, vz)
        # vector in the 2nd direction
        v3 = np.cross(v1, v2)
        # store in the instance
        self.coord_system = np.array((v1,v2,v3), dtype=np.float64).T
        # create matrix
        transformation_matrix = np.vstack((v1/np.linalg.norm(v1), v2/np.linalg.norm(v2), v3/np.linalg.norm(v3))).T

        self._v1 = transformation_matrix[0,:]
        self._v2 = transformation_matrix[1,:]
        self._v3 = transformation_matrix[2,:]

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Coordinate System: \n {coord}'.format(coord=self.coord_system)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Coordinate System: \n {coord}'.format(coord=self.coord_system)

