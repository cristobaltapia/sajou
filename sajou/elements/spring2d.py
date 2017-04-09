#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines a 2-dimensional spring element
"""
import numpy as np
import scipy.sparse as sparse
from numpy import cumsum
from sajou import loads
from sajou.elements.element import Element
from sajou.utils import Local_Csys_two_points


class Spring2D(Element):
    """Defines a 2D spring element and all its methods

    Parameters
    ----------

    node1: Node instance
        first node
    node2: Node instance
        second node
    number: int
        number of the element

    Attributes
    ----------
    node1: Node instance
        first node
    node2: Node instance
        second node
    efs: dict
        Element freedom signature. Defines the degree of freedom active (dof)
        in each node.
    nefmt: dict
        Node freedom map table of the element.
    n_active_dof: int
        number of active *dof*
    transformation_matrix: ndarray
        transformation matrix for the element, from local to global
    """

    def __init__(self, node1, node2, number):
        Element.__init__(self, number)
        self._kind = 'Spring2D'

        self._node1 = node1
        self._node2 = node2

        # Element nodal connectivity:
        self._nodal_connectivity = {0: node1, 1: node2}

        # Element Freedom Signature:
        self.efs = {
            0: np.array([1, 2, 0], dtype=np.int),
            1: np.array([1, 2, 0], dtype=np.int)
        }

        # Node freedom map table of the element (will be automatically
        # calculated when the element stiffness matrix is assebled)
        self.nefmt = None

        # Total number of active DOFs in the element (will be updated
        # when the element stiffness matrix is assembled)
        self.n_active_dof = 4
        # length of the stiffness matrix
        self._k_size = 4

        # Make the corresponding nodes aware that they belong to this
        # element instance
        node1.append_element(self, 0)
        node2.append_element(self, 1)

        # Calculate the length of the element
        self._length = np.linalg.norm(node2 - node1)
        # Local coordinate system
        self._localCSys = Local_Csys_two_points(point1=node1, point2=node2,
                                                type='cartesian')
        # Directive cosines
        delta = node2 - node1
        cx = delta[0] / self._length
        cy = delta[1] / self._length
        cz = delta[2] / self._length

        # Transformation matrix
        self.transformation_matrix = self._localCSys.calc_transformation_matrix(
            self._length, cx, cy, cz)
        # FIXME: the transformation matrix has to be a method withon the Element
        # class to avoid the following correction:
        self.transformation_matrix = self.transformation_matrix[np.ix_(
            [0, 1, 3, 4], [0, 1, 3, 4])]

        # Loads applied to the frame element
        self._loads = []
        # Initialize the loading vector
        self._load_vector_e = np.zeros(self.n_active_dof)

        # Stiffness matrix (global coordinates)
        self._Ke = None
        # Vector to calculate section forces
        # (The size of this vector depends on the order of the polynom
        # describing the section forces [...and deflections TODO])
        self._poly_sec_force = np.zeros((4, 3))

    def assign_elastic_constant(self, elastic_K):
        """Assigns the elastic constant to the spring element.

        Parameters
        ----------

        elastic_K: float
            The elastic constant of the spring

        """
        self._elastic_constant = elastic_K

        return 1

    def assemble_Ke(self):
        """
        Assemble the element stiffness matrix 'Ke' in local and global coordinates.

        :returns: local stiffness matrix

        """

        # Generate the element stiffness matrix
        k = self._assemble_Ke()
        # Generate sparse matrix
        ke_local = sparse.csr_matrix(k)
        # Store as property of the object
        self._Ke_local = ke_local
        # Transform to global coordinates:
        Te = self.transformation_matrix
        Ke = Te.T @ ke_local @ Te

        self._Ke = Ke

        # Generate the Element Node Freedom Map Table
        self._generate_element_node_freedom_map_dict()
        # Calculate total number of active DOFs in the element
        sum_dof = np.array([np.sum(nfs) for node_i, nfs in self.efs.items()])
        self.n_active_dof = np.sum(sum_dof)

        return Ke

    def _assemble_Ke(self):
        """Assemble the element stiffness matrix in local coordinates for the Bernoulli beam."""
        # Spring constant
        elastic_K = self._elastic_constant

        # Initialize stiffness matrix
        ke = np.zeros((4, 4), dtype=np.float64)

        ke[0, 0] = elastic_K
        ke[2, 2] = elastic_K
        ke[0, 2] = -elastic_K
        ke[2, 0] = -elastic_K

        return ke

    def _generate_element_node_freedom_map_dict(self):
        """
        Generate the Node Freedom Map Table of the element.

        The Node Freedom Map Table of the element is a dictionary that contains the index
        to which each node's first active DOF contributes within the element.

        It is analogous to the function __generate_node_freedom_map_dict__() from the
        Model class.

        :returns: TODO

        """
        # Obtain the number of active DOFs in each node:
        # -
        # Not sure about this. It seems that this is not really
        # necessary, since the total number of DOF and the size of the
        # element stiffness matrix is not changed. Thus, the node
        # freedom map table remains constant, even when some DOFs are
        # not used (e.g. release ends)
        # -
        # TODO: make this a dictionary
        enfmt = np.array([0, 2], dtype=np.int)

        self.enfmt = enfmt

        return enfmt
