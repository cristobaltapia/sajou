#!/usr/bin/env python
# -*- coding: utf-8 -*-


class Element(object):
    """Generic Element object. Parent class for the 2D and 3D implementation."""

    def __init__(self, number):
        """TODO: to be defined1.

        :number: number of the line

        """
        # TODO: check that the number is not in use
        self.number = number

    def __str__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam {number}: (N{n1}, N{n2})'.format(
            number=self.number, n1=self._node1.number, n2=self._node2.number)

    def __repr__(self):
        """
        Returns the printable string for this object
        """
        return 'Beam {number}'.format(number=self.number)
