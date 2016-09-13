#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains the different elements used and all its methods.
"""
import numpy as np
from .utils import Local_Csys_two_points

class Beam(object):
    """Beam object, joining two nodes. Parent class for the 2D and 3D implementation."""

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
        self.transformation_matrix = self._localCSys.calc_transformation_matrix(self._length, cx, cy, cz)

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
    """Beam object, joining two nodes"""

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
        self.dof3 = 0. # rot z

        # Node 2
        self.dof4 = 0. # trans x
        self.dof5 = 0. # trans y
        self.dof6 = 0. # rot z

        # Release rotation on the ends of the beam
        self.release_node_1 = False # first node
        self.release_node_2 = False # second node

        self._Ke = None

    def assemble_K(self, second_order=False):
        """This function assembles the stiffness matrix for one individual element. Optionally
        it can take the shear effect into account (second order effect).

        :second_order: boolean
        :returns: local stiffness matrix

        """
        # Modulus of elasticity
        E = self._beam_section._material._data[0]
        # Area of the section
        EA = self._beam_section._area * E
        # Inertias
        EI = self._beam_section._Iz * E
        # Length of the element
        L = self._length
        # TODO: Take account of the second order effects
        if second_order:
            G = self._material.G
            Ksy = 12.*E*Iz / (G*Asy*L*L)
            Ksz = 12.*E*Iy / (G*Asz*L*L)
        else:
            Ksy = 0.0
            Ksz = 0.0

        # Initialize stiffness matrix
        k = np.zeros((6,6), dtype=np.float64)

        k[0,0] = k[3,3] = EA / L
        k[1,1] = k[4,4] = 12. * EI / (L*L*L)
        k[2,2] = k[5,5] = 4. * EI / L
        k[2,1] = k[1,2] = 6. * EI / L**2
        k[3,0] = k[0,3] = - EA / L
        k[4,1] = k[1,4] = -12. * EI / (L*L*L)
        k[4,2] = k[2,4] = -6. * EI / L**2
        k[5,1] = k[1,5] = 6. * EI / L**2
        k[5,2] = k[2,5] = 2. * EI / L
        k[5,4] = k[4,5] = -6. * EI / L**2

        # transform to global coordinates
        #T = element.transformation_matrix

        #Ke = np.dot(T.T, np.dot(k,T))
        self._Ke_local = k

        # transform to global coordinates
        T = self.transformation_matrix

        Ke = np.dot(T.T, np.dot(k, T))

        self._Ke = Ke

        return k

    def assemble_sym_K(self):
        """This function assembles the stiffness matrix for one individual element. Optionally
        it can take the shear effect into account (second order effect).

        :returns: local stiffness matrix

        """
        from sympy.matrices import Matrix, zeros
        # Modulus of elasticity
        E = self._beam_section._material._data[0]
        # Area of the section
        EA = self._beam_section._area * E
        # Inertias
        EI = self._beam_section._Iz * E
        # Length of the element
        L = self._length

        # Initialize stiffness matrix
        k = zeros(6,6)

        k[0,0] = k[3,3] = EA / L
        k[1,1] = k[4,4] = 12. * EI / (L*L*L)
        k[2,2] = k[5,5] = 4. * EI / L
        k[2,1] = k[1,2] = 6 * EI / L**2
        k[3,0] = k[0,3] = - EA / L
        k[4,1] = k[1,4] = -12. * EI / (L*L*L)
        k[4,2] = k[2,4] = -6. * EI / L**2
        k[5,1] = k[1,5] = 6. * EI / L**2
        k[5,2] = k[2,5] = 2. * EI / L
        k[5,4] = k[4,5] = -6. * EI / L**2

        # transform to global coordinates
        #T = element.transformation_matrix

        #Ke = np.dot(T.T, np.dot(k,T))
        self._Ke_local = k

        # transform to global coordinates
        T = Matrix(self.transformation_matrix)

        aux = k.multiply(T)
        TT = T.T
        Ke = TT.multiply(aux)

        self._Ke = Ke

        return k

class Beam3D(Beam):
    """Beam object, joining two nodes"""

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

    def assemble_K(self, second_order=False):
        """This function assembles the stiffness matrix for one individual element. Optionally
        it can take the shear effect into account (second order effect).

        :element: segment instance
        :seecond_order: boolean
        :returns: local stiffness matrix

        """
        # Modulus of elasticity
        E = self._beam_section.material._data[0]
        # Area of the section
        Ax = self._beam_section._area
        # Shear areas
        Asy = self._beam_section._Sy
        Asz = self._beam_section._Sz
        # Inertias
        Iz = self._beam_section._I33
        Iy = self._beam_section._I22
        J = self._beam_section._J11
        # Length of the element
        Le = self._length
        # Take account of the second order effects
        if second_order:
            G = self._material.G
            Ksy = 12.*E*Iz / (G*Asy*Le*Le)
            Ksz = 12.*E*Iy / (G*Asz*Le*Le)
        else:
            Ksy = 0.0
            Ksz = 0.0

        # Initialize stiffness matrix
        k = np.zeros((12,12))

        k[0,0] = k[6,6]   = E*Ax / Le
        k[1,1] = k[7,7]   = 12.*E*Iz / ( Le*Le*Le*(1.+Ksy) )
        k[2,2] = k[8,8]   = 12.*E*Iy / ( Le*Le*Le*(1.+Ksz) )
        k[3,3] = k[9,9]   = G*J / Le
        k[4,4] = k[10,10] = (4.+Ksz)*E*Iy / ( Le*(1.+Ksz) )
        k[5,5] = k[11,11] = (4.+Ksy)*E*Iz / ( Le*(1.+Ksy) )

        k[4,2] = k[2,4] = -6.*E*Iy / ( Le*Le*(1.+Ksz) );
        k[5,1] = k[1,5] = 6.*E*Iz / ( Le*Le*(1.+Ksy) );
        k[6,0] = k[0,6] = -k[0,0];

        k[11,7] = k[7,11] = k[7,5] = k[5,7] = -k[5,1];
        k[10,8] = k[8,10] = k[8,4] = k[4,8] = -k[4,2];
        k[9,3]  = k[3,9]  = -k[3,3];
        k[10,2] = k[2,10] = k[4,2];
        k[11,1] = k[1,11] = k[5,1];

        k[7,1]  = k[1,7]  = -k[1,1];
        k[8,2]  = k[2,8]  = -k[2,2];
        k[10,4] = k[4,10] = (2.-Ksz)*E*Iy / ( Le*(1.+Ksz) );
        k[11,5] = k[5,11] = (2.-Ksy)*E*Iz / ( Le*(1.+Ksy) );

        # transform to global coordinates
        T = self._localCSys.transformation_matrix


        return k

