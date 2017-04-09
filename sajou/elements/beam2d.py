#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines a 2-dimensional Bernoulli beam element
"""
import numpy as np
import scipy.sparse as sparse
from numpy import cumsum
from sajou import loads
from sajou.elements.element import Element
from sajou.utils import Local_Csys_two_points


class Beam2D(Element):
    """
    Two-dimensional Bernoulli Beam element.

    This beam element connects two nodes and is based on the Bernoulli beam theory.
    Optionally the rotations on the end nodes can be released to create hinges.

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
    release_end_1: bool
        Release rotation on the node 1 of the beam
    release_end_2: bool
        Release rotation on the node 2 of the beam
    transformation_matrix: ndarray
        transformation matrix for the element, from local to global
    _k_size: int
        size of the stiffness matrix
    _nodal_connectivity: dict
        defines the order of the nodes in the element
    _length: float
        length of the element
    _localCSys: Csys
        local coordinate system of the element
    _Ke: ndarray
        element stiffness matrix (local coordinates)
    _load_vector_e: ndarray
        Load vector (global coordinates)
    _poly_sec_force: ndarray
        Vector to calculate section forces
    _beam_section: BeamSection
        Beam section
    _loads: list[Load]
        Loads applied to the frame element

    Todo
    ----

    fixed problem with release ends and distribuetd loads

    """

    def __init__(self, node1, node2, number):
        # Instatiate the Element parent class
        Element.__init__(self, number)
        self._kind = 'Beam2D'
        # TODO: accept tuples with coordinates also
        self._node1 = node1
        self._node2 = node2

        # Element nodal connectivity:
        self._nodal_connectivity = {0: node1, 1: node2}

        # Element Freedom Signature:
        self.efs = {
            0: np.array([1, 2, 3], dtype=np.int),
            1: np.array([1, 2, 3], dtype=np.int)
        }

        # Node freedom map table of the element (will be automatically
        # calculated when the element stiffness matrix is assebled)
        self.nefmt = None

        # Total number of active DOFs in the element (will be updated
        # when the element stiffness matrix is assembled)
        self.n_active_dof = 6
        # length of the stiffness matrix
        self._k_size = 6

        # Make the corresponding nodes aware that they belong to this
        # element instance
        node1.append_element(self, 0)
        node2.append_element(self, 1)

        # Calculate the length of the element
        self._length = np.linalg.norm(node2 - node1)
        # Local coordinate system
        self._localCSys = Local_Csys_two_points(
            point1=node1, point2=node2, type='cartesian')
        # Directive cosines
        delta = node2 - node1
        cx = delta[0] / self._length
        cy = delta[1] / self._length
        cz = delta[2] / self._length

        # Transformation matrix
        self.transformation_matrix = self._localCSys.calc_transformation_matrix(
            self._length, cx, cy, cz)

        # Stiffness matrix (global coordinates)
        self._Ke = None
        # Release rotation on the ends of the beam
        self.release_end_1 = False  # first node
        self.release_end_2 = False  # second node

        # Beam section
        self._beam_section = None
        # Loads applied to the frame element
        self._loads = []
        # Initialize the loading vector
        self._load_vector_e = np.zeros(self.n_active_dof)
        # Vector to calculate section forces
        # (The size of this vector depends on the order of the polynom
        # describing the section forces [...and deflections TODO])
        self._poly_sec_force = np.zeros((4, 3))

    def assemble_Ke(self, second_order=False):
        """
        Assemble the element stiffness matrix 'Ke' in local and global coordinates.

        :second_order: boolean
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
        sum_dof = np.array([np.sum((nfs >= 1)) for node_i, nfs in self.efs.items()])
        self.n_active_dof = np.sum(sum_dof)

        return Ke

    def distributed_load(self, **kwargs):
        """Assign a DistributedLoad object to the frame current element.

        The parameters are the same as defined for the class DistributedLoad()

        Parameters
        ----------

        p1: float
            value of the distributed load at the start node
        p2: float
            value of the distributed load at the end node
        direction: str
            direction of load application. Can be 'x' or 'y'
        coord_system:
            coordinate system to which the *direction* parameter applies. It can
            be 'local' or 'global'

        Returns
        -------

        a DistributedLoad instance:

        """
        dist_load = loads.DistributedLoad(elem=self, **kwargs)

        # Add this DistributedLoad instance to the list of loads of the
        # element
        self._loads.append(dist_load)
        # Append the load vector (in global coordinates)
        self._load_vector_e += dist_load._load_vector_global
        self._poly_sec_force += dist_load._poly_sec_force

        return 1

    def distributed_moment(self, **kwargs):
        """Assign a DistributedLoad object to the frame current element.

        Parameters
        ----------

        m1: float
            moment applied at the start end of the beam
        m2: float
            moment applied at the end of the beam
        direction: str
            direction of application of the moment. Only 'z' is allowed in Beam2D
        coord_system: str
            coordinate system to which the 'direction' parameter applies

        Returns
        -------

        returns: TODO

        Note
        ----

        The parameters are the same as defined for the class DistributedMoment()

        """
        dist_moment = loads.DistributedMoment(elem=self, **kwargs)

        # Add this DistributedLoad instance to the list of loads of the
        # element
        self._loads.append(dist_moment)
        # Append the load vector (in global coordinates)
        self._load_vector_e += dist_moment._load_vector_global
        self._poly_sec_force += dist_moment._poly_sec_force

        return 1

    def release_end(self, which):
        """
        Release the rotation DOF of one or both ends of the beam element.

        Parameters
        ----------

        which: int, str
            specifies the end that should be released. It can be '1', '2' or 'both'

        Returns
        -------
        TODO: bool

        """
        # Set the respective property to 'True'
        if which == 1:
            self.release_end_1 = True
            n_new_dof = self._nodal_connectivity[0]._add_dof(self, 1)
            self.efs[0][2] = n_new_dof
        elif which == 2:
            self.release_end_2 = True
            n_new_dof = self._nodal_connectivity[1]._add_dof(self, 1)
            self.efs[1][2] = n_new_dof
        elif which == 'both':
            self.release_end_1 = True
            self.release_end_2 = True
            n_new_dof_1 = self._nodal_connectivity[0]._add_dof(self, 1)
            n_new_dof_2 = self._nodal_connectivity[1]._add_dof(self, 1)
            self.efs[0][2] = n_new_dof_1
            self.efs[1][2] = n_new_dof_2

        return 1

    def assign_section(self, beam_section):
        """Assign a beam section instance to the beam

        Parameters
        ----------

        beam_section: BeamSection instance
            section to be assigned

        Returns
        -------
        self: the same Beam2D instance

        """
        self._beam_section = beam_section

        return self

    def _assemble_Ke(self):
        """
        Assemble the element stiffness matrix in local coordinates for the Bernoulli beam.
        :returns: TODO

        """
        # Modulus of elasticity
        E = self._beam_section._material._data[0]
        # Area of the section
        EA = self._beam_section._area * E
        # Inertias
        EI = self._beam_section._Iz * E
        # Length of the element
        L = self._length

        # Initialize stiffness matrix
        ke = np.zeros((6, 6), dtype=np.float64)

        ke[0, 0] = ke[3, 3] = EA / L
        ke[1, 1] = ke[4, 4] = 12. * EI / (L * L * L)
        ke[2, 2] = ke[5, 5] = 4. * EI / L
        ke[2, 1] = ke[1, 2] = 6. * EI / L**2
        ke[3, 0] = ke[0, 3] = -EA / L
        ke[4, 1] = ke[1, 4] = -12. * EI / (L * L * L)
        ke[4, 2] = ke[2, 4] = -6. * EI / L**2
        ke[5, 1] = ke[1, 5] = 6. * EI / L**2
        ke[5, 2] = ke[2, 5] = 2. * EI / L
        ke[5, 4] = ke[4, 5] = -6. * EI / L**2

        return ke

    def _generate_element_node_freedom_map_dict(self):
        """
        Generate the Node Freedom Map Table of the element.

        The Node Freedom Map Table of the element is a dictionary that contains
        the index to which each node's first active DOF contributes to within
        the element.

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
        #n_active_dof = [sum(nfs) for n_node, nfs in self.efs.items()]
        # Obtain the cumulative sum
        #enfmt = cumsum(n_active_dof, dtype=np.int) - n_active_dof[0]
        # TODO: make this a dictionary
        enfmt = np.array([0, 3], dtype=np.int)

        self.enfmt = enfmt

        return enfmt

