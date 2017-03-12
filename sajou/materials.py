#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This module contains the classes and methods used for the creation of material
properties.
"""


class Material(object):
    """Material properties

    Parameters
    ----------

    name: str
        name of the material
    table: tuple
        properties of the material
    type: str
        type of the material:

        - 'isotropic': data = (E, )
        - 'orthotropic': data = (E_1, E_2, E_3)

    Todo
    ----

    create function to print information of the material

    """

    def __init__(self, name, data, type='isotropic'):
        self._name = name
        self._data = data
        self._type = type

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Material: {name}'.format(name=self._name)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Material: {name}'.format(name=self._name)
