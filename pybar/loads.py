#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Defines the different load types,thta can be applied to the elements, which are
transfered accordingly to the respective nodes.
"""
import numpy as np
import scipy.sparse as sparse
from .utils import Local_Csys_two_points

class Load(object):

    """Defines the Load object."""

    def __init__(self):
        """Initialize the Load instance"""
        self._type = ''
        self._load_vector_global = None
        
class DistributedLoad(Load):

    """Docstring for DistributedLoad. """

    def __init__(self, elem, p1, p2=None, direction='y', coord_system='local'):
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
        # Assign a coordinate system to the load
        if coord_system == 'local':
            self._localCSys = elem._localCSys
        elif coord_system == 'global':
            self._localCSys = Local_Csys_two_points(point1=(0.,0.,0.), point2=(1.,0.,0.))

        # Detect if distribution is a varying distributed load or not
        if p2 == None:
            self.is_uniform = True
            p2 = p1
            self._p2 = p1
        else:
            self.is_uniform = False
            self._p2 = p2

        # If the distributed load is given in local coordinates
        if coord_system == 'local':
            # Generate loading vector
            load_v, poly_sec_force = self._calc_loading_vector_local(
                        p1,
                        p2,
                        elem._length,
                        direction)
        # Else, if the distributed load is given in global coordinates
        elif coord_system == 'global':
            # Generate loading vector
            load_v, poly_sec_force = self._calc_loading_vector_global(
                    p1,
                    p2,
                    elem._length,
                    direction)

        self._load_vector_global = load_v
        self._poly_sec_force = poly_sec_force

    def _calc_loading_vector_local(self, p1, p2, length, direction):
        """ Generate the loading vector, when the distributed load is in local coords.
        Also returns the matrix used for the calculation of the sectional forces.

        :returns: TODO

        """
        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        n_dof = self._elem._ndof
        load_v = np.zeros(n_dof)

        # Load vector for the axial load
        # (direction='x')
        if direction == 'x':
            load_v[0] = length * (2.*p1 + p2) / 6.
            load_v[3] = length * (p1 + 2.*p2) / 6.
            # Generate matrix used for the calculation of section forces
            poly_sec_force = self._generate_section_force_poly(
                    p1,
                    p2,
                    length,
                    direction)
        # Load vector for the transversal load
        # (direction='z')
        elif direction == 'y':
            load_v[1] = length * (7.*p1 + 3.*p2) / 20.
            load_v[2] = length**2 * (p1/20. + p2/30.)
            load_v[4] = length * (3.*p1 + 7*p2) / 20.
            load_v[5] = -length**2 * (p1/30. + p2/20.)
            # Generate matrix used for the calculation of section forces
            poly_sec_force = self._generate_section_force_poly(
                    p1,
                    p2,
                    length,
                    direction)

        self._loading_vector = load_v
        # Calculate the load vector in global coordinates, using the
        # transformation matrix
        T = self._elem.transformation_matrix
        # Rotate
        load_vector_global = T.T.dot(load_v)

        return load_vector_global, poly_sec_force

    def _calc_loading_vector_global(self, p1, p2, length, direction):
        """ Generate the loading vector, when the distributed load is in global coords.

        :returns: TODO

        """
        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        n_dof = self._elem._ndof
        load_v = np.zeros(n_dof)
        poly_sec_force = np.zeros((4, 3))

        # transformation matrix
        T = self._elem.transformation_matrix

        # the load has to be decomposed in their respective local
        # components:
        if direction == 'x':
            # x-component in local coordinates
            p1_x = p1 * T[0,0]
            p2_x = p2 * T[0,0]
            load_v_aux, poly_sec_force_aux = self._calc_loading_vector_local(p1_x, p2_x, length, 'x')
            load_v += load_v_aux
            poly_sec_force += poly_sec_force_aux
            # y-component in local coordinates
            p1_y = p2 * T[0,1]
            p2_y = p1 * T[0,1]
            load_v_aux, poly_sec_force_aux = self._calc_loading_vector_local(p1_y, p2_y, length, 'y')
            load_v += load_v_aux
            poly_sec_force += poly_sec_force_aux

        elif direction == 'y':
            # x-component in local coordinates
            p1_x = p1 * T[0,1]
            p2_x = p2 * T[0,1]
            load_v_aux, poly_sec_force_aux = self._calc_loading_vector_local(p1_x, p2_x, length, 'x')
            load_v += load_v_aux
            poly_sec_force += poly_sec_force_aux
            # y-component in local coordinates
            p1_y = p2 * T[0,0]
            p2_y = p1 * T[0,0]
            load_v_aux, poly_sec_force_aux = self._calc_loading_vector_local(p1_y, p2_y, length, 'y')
            load_v += load_v_aux
            poly_sec_force += poly_sec_force_aux

        return load_v, poly_sec_force

    def _generate_section_force_poly(self, p1, p2, length, direction):
        """Generate the matrix used to calculate the section forces.

        This matrix has a shape (4xn_dof) and will be used to calculate the sectional
        forces produced by this LoadDistribution instance.
        It will then be added in the Element object to contain every contribution made to
        the element.

                [ N ]  
                [ V ] = [ 1, x, x**2, x**3 ] * S
                [ M ]

            where 'S' is the matrix created here, 'x' is the position along the element in
            local direction 'x'.

        :p1: TODO
        :p2: TODO
        :length: TODO
        :direction: TODO
        :returns: TODO

        """
        # TODO: implement for the 3D case
        # Initialize matrix
        m_sec_force = np.zeros((4, 3))
        # Determine in which case we are
        if direction == 'x':
            m_sec_force[:, 0] = np.array([0., -p1, (p1-p2)/(2*length), 0. ])
        # For the case in which the loading direction is 'y'
        elif direction == 'y':
            m_sec_force[:, 1] = np.array([0., p1, (p2-p1)/(2*length), 0.                 ])
            m_sec_force[:, 2] = np.array([0., 0., p1*0.5,             (p2-p1)/(6*length) ])

        return m_sec_force

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

        if coord_system == 'local':
            # Generate loading vector
            load_v, poly_sec_force = self._calc_loading_vector_local(
                        m1,
                        m2,
                        elem._length,
                        direction)
        # Else, if the distributed load is given in global coordinates
        elif coord_system == 'global':
            # Generate loading vector
            load_v, poly_sec_force = self._calc_loading_vector_global(
                    m1,
                    m2,
                    elem._length,
                    direction)

        self._load_vector_global = load_v
        self._poly_sec_force = poly_sec_force

    def _calc_loading_vector_local(self, m1, m2, length, direction):
        """ Generate the loading vector, when the distributed load is in local coords.
        Also returns the matrix used for the calculation of the sectional forces.

        :returns: TODO

        """
        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        n_dof = self._elem._ndof
        load_v = np.zeros(n_dof)

        # Load vector for the moment load
        # (direction='z')
        # FIXME:
        if direction == 'z':
            load_v[1] = -(m1 + m2) * 0.5
            load_v[2] = -length * (m1 - m2) / 12.
            load_v[4] = (m1 + m2) * 0.5
            load_v[5] = length * (m1 - m2) / 12.
            # Generate matrix used for the calculation of section forces
            poly_sec_force = self._generate_section_force_poly(
                    m1,
                    m2,
                    length,
                    direction)

        self._loading_vector = load_v
        # Calculate the load vector in global coordinates, using the
        # transformation matrix
        T = self._elem.transformation_matrix
        # Rotate
        load_vector_global = T.T.dot(load_v)

        return load_vector_global, poly_sec_force

    def _calc_loading_vector_global(self, m1, m2, length, direction):
        """ Generate the loading vector, when the distributed load is in global coords.

        :returns: TODO

        """
        # Initialize loading vector
        # FIXME: make this dependant from the specific element.
        # (thinking in 3D case)
        n_dof = self._elem._ndof
        load_v = np.zeros(n_dof)
        poly_sec_force = np.zeros((4, 3))

        # transformation matrix
        T = self._elem.transformation_matrix

        # the load has to be decomposed in their respective local
        # components:
        if direction == 'z':
            # TODO: 3D case
            # z-component in local coordinates
            load_v_aux, poly_sec_force_aux = self._calc_loading_vector_local(m1, m2, length, 'z')
            load_v += load_v_aux
            poly_sec_force += poly_sec_force_aux
            # y-component in local coordinates
            # ...

        elif direction == 'y':
            # x-component in local coordinates
            # TODO: 3D case
            pass

        return load_v, poly_sec_force

    def _generate_section_force_poly(self, m1, m2, length, direction):
        """Generate the matrix used to calculate the section forces.

        This matrix has a shape (4xn_dof) and will be used to calculate the sectional
        forces produced by this LoadDistribution instance.
        It will then be added in the Element object to contain every contribution made to
        the element.

                [ N ]  
                [ V ] = [ 1, x, x**2, x**3 ] * S
                [ M ]

            where 'S' is the matrix created here, 'x' is the position along the element in
            local direction 'x'.

        :p1: TODO
        :p2: TODO
        :length: TODO
        :direction: TODO
        :returns: TODO

        """
        # TODO: implement for the 3D case
        # Initialize matrix
        m_sec_force = np.zeros((4, 3))

        # Determine in which case we are
        if direction == 'z':
            # Effect over the moment
            m_sec_force[:, 2] = np.array([0., m1 + (m2-m1)/length, 0., 0. ])
        # For the case in which the loading direction is 'y'
        elif direction == 'y':
            # TODO:
            pass

        return m_sec_force

