#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In this module functions for handling the stiffness matrices are created.
"""
import numpy as np

def assemble_Ke_2D(element, second_order=False):
    """This function assembles the stiffness matrix for one individual element. Optionally
    it can take the shear effect into account (second order effect).

    :element: segment instance
    :seecond_order: boolean
    :returns: local stiffness matrix

    """
    # Modulus of elasticity
    E = element._beam_section._material._data[0]
    # Area of the section
    EA = element._beam_section._area * E
    # Inertias
    EI = element._beam_section._Iz * E
    # Length of the element
    L = element._length
    # TODO: Take account of the second order effects
    if second_order:
        G = element._material.G
        Ksy = 12.*E*Iz / (G*Asy*L*L)
        Ksz = 12.*E*Iy / (G*Asz*L*L)
    else:
        Ksy = 0.0
        Ksz = 0.0

    # Initialize stiffness matrix
    k = np.zeros((6,6))

    k[0,0] = k[3,3] = EA / L
    k[1,1] = k[4,4] = 12. * EI / L**3
    k[2,2] = k[5,5] = 4. * EI / L
    k[2,1] = k[1,2] = 6 * EI / L**2
    k[3,0] = k[0,3] = - EA / L
    k[4,1] = k[1,4] = -12. * EI / L**3
    k[4,2] = k[2,4] = -6. * EI / L**2
    k[5,1] = k[1,5] = 6. * EI / L**2
    k[5,2] = k[2,5] = 2. * EI / L
    k[5,4] = k[4,5] = -6. * EI / L**2

    # transform to global coordinates
    #T = element.transformation_matrix

    #Ke = np.dot(T.T, np.dot(k,T))

    return k

def assemble_Ke_3D(element, second_order=False):
    """This function assembles the stiffness matrix for one individual element. Optionally
    it can take the shear effect into account (second order effect).

    :element: segment instance
    :seecond_order: boolean
    :returns: local stiffness matrix

    """
    # Modulus of elasticity
    E = element._beam_section.material._data[0]
    # Area of the section
    Ax = element._beam_section._area
    # Shear areas
    Asy = element._beam_section._Sy
    Asz = element._beam_section._Sz
    # Inertias
    Iz = element._beam_section._I33
    Iy = element._beam_section._I22
    J = element._beam_section._J11
    # Length of the element
    Le = element._length
    # Take account of the second order effects
    if second_order:
        G = element._material.G
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
    T = element._localCSys.transformation_matrix


    return k
