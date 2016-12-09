#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse as sparse
from pybar import loads
from pybar.utils import Local_Csys_two_points
from pybar.elements.element import Element

class Beam2D(Element):
    """
    Two-dimensional Bernoulli Beam element.

    This beam element connects two nodes and is based on the Bernoulli beam theory.
    Optionally the rotations on the end nodes can be released to create hinges.

    TODO: fixed problem with release ends and distribuetd loads
    """

    def __init__(self, node1, node2, number):
        """TODO: to be defined1.

        :node1: first node
        :node2: second node
        :number: number of the line

        """
        # Instatiate the Element parent class
        Element.__init__(self, number)
        # TODO: accept tuples with coordinates also
        self._node1 = node1
        self._node2 = node2

        # Element nodal connectivity:
        self._nodal_connectivity = {
                0:node1,
                1:node2
                }

        # Element Freedom Signature:
        self.efs = {
                0: np.array([1, 1, 1], dtype=np.int),
                1: np.array([1, 1, 1], dtype=np.int)
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
        self._length = np.linalg.norm(node2-node1)
        # Local coordinate system
        self._localCSys = Local_Csys_two_points(point1=node1, point2=node2, type='cartesian')
        # Directive cosines
        delta = node2 - node1
        cx = delta[0] / self._length
        cy = delta[1] / self._length
        cz = delta[2] / self._length

        # Transformation matrix
        self.transformation_matrix = self._localCSys.calc_transformation_matrix(self._length, cx, cy, cz)

        # Stiffness matrix (global coordinates)
        self._Ke = None
        # Load vector (global coordinates)
        self._load_vector_e = None
        # Vector to calculate section forces
        self._poly_sec_force = None
        # Release rotation on the ends of the beam
        self.release_end_1 = False # first node
        self.release_end_2 = False # second node

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

        # Generate the element stiffness matrix depending on the use of
        # hinges or not (release_end_1 and release_end_2)
        if self.release_end_1 == False and self.release_end_2 == False:
            k = self._assemble_Ke()
        # If both ends are released
        elif self.release_end_1 == True and self.release_end_2 == True:
            k = self._assemble_Ke_end_release_both()
        # If the rotation is release in the first node
        elif self.release_end_1 == True:
            k = self._assemble_Ke_end_release1()
        # If the rotation is release in the second node
        elif self.release_end_2 == True:
            k = self._assemble_Ke_end_release2()

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
        sum_dof = np.array( [np.sum(nfs) for node_i, nfs in self.efs.items()] )
        self.n_active_dof = np.sum(sum_dof)

        return Ke

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
        ke = np.zeros((6,6), dtype=np.float64)

        ke[0,0] = ke[3,3] = EA / L
        ke[1,1] = ke[4,4] = 12. * EI / (L*L*L)
        ke[2,2] = ke[5,5] = 4. * EI / L
        ke[2,1] = ke[1,2] = 6. * EI / L**2
        ke[3,0] = ke[0,3] = - EA / L
        ke[4,1] = ke[1,4] = -12. * EI / (L*L*L)
        ke[4,2] = ke[2,4] = -6. * EI / L**2
        ke[5,1] = ke[1,5] = 6. * EI / L**2
        ke[5,2] = ke[2,5] = 2. * EI / L
        ke[5,4] = ke[4,5] = -6. * EI / L**2

        return ke

    def _assemble_Ke_end_release1(self):
        """
        Assembles the element stiffness matrix in local and global coordinates for the
        case when the release_end_1 option is set to 'True'.
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
        ke = np.zeros((6,6), dtype=np.float64)

        ke[0,0] = ke[3,3] = EA / L
        ke[1,1] = ke[4,4] = 3. * EI / (L*L*L)
        ke[5,5] = 3. * EI / L
        ke[3,0] = ke[0,3] = - EA / L
        ke[4,1] = ke[1,4] = -3. * EI / (L*L*L)
        ke[5,1] = ke[1,5] = 3. * EI / L**2
        ke[5,4] = ke[4,5] = -3. * EI / L**2

        # The Element Freedom signature has to be changed for the
        # respective node
        self.efs[0] = np.array([1, 1, 0], dtype=np.int)

        return ke

    def _assemble_Ke_end_release2(self):
        """
        Assembles the element stiffness matrix in local and global coordinates for the
        case when the release_end_2 option is set to 'True'.
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
        ke = np.zeros((6,6), dtype=np.float64)

        ke[0,0] = ke[3,3] = EA / L
        ke[1,1] = ke[4,4] = 3. * EI / (L*L*L)
        ke[2,2] = 3. * EI / L
        ke[2,1] = ke[1,2] = 3. * EI / L**2
        ke[2,4] = ke[4,2] = -3. * EI / L**2
        ke[3,0] = ke[0,3] = - EA / L
        ke[4,1] = ke[1,4] = -3. * EI / (L*L*L)

        # The Element Freedom signature has to be changed for the
        # respective node
        self.efs[1] = np.array([1, 1, 0], dtype=np.int)

        return ke

    def _assemble_Ke_end_release_both(self):
        """
        Assembles the element stiffness matrix in local and global coordinates for the
        case when the release_end_2 option is set to 'True'.
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
        ke = np.zeros((6,6), dtype=np.float64)

        ke[0,0] = ke[3,3] = EA / L
        ke[3,0] = ke[0,3] = - EA / L

        # The Element Freedom signature has to be changed for the
        # respective node
        self.efs[0] = np.array([1, 1, 0], dtype=np.int)
        self.efs[1] = np.array([1, 1, 0], dtype=np.int)

        return ke

    def _calc_condensed_displacements(self, displ):
        """
        Calculate the nodal displacements taking into account the condensend DOF due to
        the end_release options.

        If the condensed displacemtns are denoted by the subindex 'j' and the remaining
        DOFs are denoted with the letter 'i', then the condensed displacements are
        calculated according to


                v_j = -K_jj^-1 . K_ji . v_i

        :displ: TODO
        :returns: TODO

        """
        if self.release_end_1 == True and self.release_end_2 == True:
            j_n = [2, 5]
            i_n = [0, 1, 3, 4]
        elif self.release_end_1 == True:
            j_n = [2]
            i_n = [0, 1, 3, 4, 5]
        elif self.release_end_2 == True:
            j_n = [5]
            i_n = [0, 1, 2, 3, 4]

        ke = self._assemble_Ke()
        # To sparse
        ke_local = sparse.csr_matrix(ke)
        # Transform to global coordinates:
        Te = self.transformation_matrix
        Ke = Te.T @ ke_local.todense() @ Te
        # Transform to global coordinates
        k_jj = Ke[j_n,j_n]
        k_ji = Ke[j_n,i_n]
        from numpy.linalg import inv
        #print(k_jj.todense())
        #print(k_ji.todense())
        v_j = k_jj**(-1) @ k_ji @ displ[i_n]

        return v_j

    def get_index_array_of_node(self, node):
        """
        Get an array containing the indices of the used DOFs of the given node of the
        element.

        The array has the following form:
                [0, 1, 2] --> all DOFs are used
                [0, 2] --> DOF 0 and 2 are used only (ux and rz)

        This array is used to know exactly which DOFs should be used to assemble the global
        stiffness matrix or to retrieve the corresponding displacements.

        For example:

            i_global_node_1 = e.get_index_list_of_node(n_node_ele) + nfat[global_node_number]

        :node: the number of the node in the element (element number of the node: 0, 1,
        2... )
        :returns: array of indices

        """
        efs = self.efs[node]

        return np.arange(len(efs))[efs>0]

    def _generate_element_node_freedom_map_dict(self):
        """
        Generate the Node Freedom Map Table of the element.

        The Node Freedom Map Table of the element is a dictionary that contains the index 
        to which each node's first active DOF contributes within the element.

        It is analogous to the function __generate_node_freedom_map_dict__() from the
        Model class.

        :returns: TODO

        """
        from numpy import cumsum
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

        self._Ke_local = k

        # transform to global coordinates
        T = Matrix(self.transformation_matrix)

        aux = k.multiply(T)
        TT = T.T
        Ke = TT.multiply(aux)

        self._Ke = Ke

        return k

    def distributed_load(self, **kwargs):
        """Assign a DistributedLoad object to the frame current element.

        The parameters are the same as defined for the class DistributedLoad()

        :p1: TODO
        :p2: TODO
        :direction: TODO
        :coord_system: TODO
        :returns: TODO

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

        The parameters are the same as defined for the class DistributedMoment()

        :p1: TODO
        :p2: TODO
        :direction: TODO
        :coord_system: TODO
        :returns: TODO

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

        :which: specifies the end that should be released. It can be '1', '2' or 'both'
        :returns: TODO

        """
        # Set the respective property to 'True'
        if which == 1:
            self.release_end_1 = True
        elif which == 2:
            self.release_end_2 = True
        elif which == 'both':
            self.release_end_1 = True
            self.release_end_2 = True

        return 1

    def assign_section(self, beam_section):
        """Assign a beam section instance to the beam

        :beam_section: a BeamSection instance
        :returns: self

        """
        self._beam_section = beam_section

        return self

