#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Define the classes and methods to work with sections.
"""
import numpy as np


class BeamSection(object):
    """Defines a beam section"""

    def __init__(self, name, material, data, type='rectangular'):
        """Initiates the section

        :name: name of the section
        :material: material of the section defined as an instance of Material object
        :data: properties of the section:
        :type: defines the type of cross-section
                - rectangular: data=(width, height,)
                - circular:    data=(r, )
                - I-section:  data=(H, h_f, w_web, w_f)
                - general:    data=(A, I_3,)

        """
        self._name = name
        self._material = material
        self._data = data
        self._type = type
        self._area = 0
        self._Iz = 0
        self._Iy = 0
        self._Jx = 0
        self.compute_properties()

    def print_properties(self):
        """Prints the properties of the BeamSection instance
        :returns: TODO

        """
        if self._type == 'rectangular':
            props = {'width': self._data[0], 'height': self._data[1]}
        else:
            props = 'undefined'

        return 'Properties: ' + str(props)

    def compute_properties(self):
        """Compute all the mechanical properties for the given section
        :returns: TODO

        """
        # Calculate the area
        self._area = self.calc_area()
        self._Iz, self._Iy = self.calc_inertia()

    def calc_area(self):
        """Calculate the area of the section
        :returns: TODO

        """
        type = self._type

        if type == 'rectangular':
            width = self._data[0]
            height = self._data[1]
            return width * height
        elif type == 'general':
            return self._data[0]

        elif type == 'circular':
            radius = self._data[0]
            return np.pi * radius**2

    def calc_inertia(self):
        """Calculate the moment of inertia of the beam section
        :returns: Iz, Iy

        """
        type = self._type

        if type == 'rectangular':
            width = self._data[0]
            height = self._data[1]
            I_z = width * height**3 / 12.
            I_y = height * width**3 / 12.
            return I_z, I_y
        elif type == 'general':
            return self._data[1], 0

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam Section: {name}, type: {t}'.format(name=self._name,
                                                        t=self._type)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam Section: {name}, type: {t}'.format(name=self._name,
                                                        t=self._type)
