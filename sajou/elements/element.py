#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np


class Element(object):
    """Generic Element object. Parent class for the 2D and 3D implementation.

    Parameters
    ----------

    number: int
        number identifying the element in the model

    """

    def __init__(self, number):
        # TODO: check that the number is not in use
        self.number = number

    def get_index_array_of_node(self, node):
        """
        Get an array containing the indices of the used DOFs of the given node of the
        element.

        Parameters
        ----------

        node: Node instance
            the number of the node in the element (element number of the node: 0, 1,
            2, ... )

        Returns
        -------

        array of indices:

        Note
        ----

        The array has the following form::

                [0, 1, 2] --> all DOFs are used
                [0, 2] --> DOF 0 and 2 are used only (ux and rz)

        This array is used to know exactly which DOFs should be used to assemble the global
        stiffness matrix or to retrieve the corresponding displacements.

        Example
        ---------

        It would be used like this::

            i_global_node_1 = e.get_index_list_of_node(n_node_ele) + nfat[global_node_number]


        """
        efs = self.efs[node]

        return np.arange(len(efs))[efs > 0]

    def __str__(self):
        """Returns the printable string for this object"""
        return 'Element {number}: (N{n1}, N{n2}), type: {kind}'.format(
            number=self.number, n1=self._node1.number, n2=self._node2.number,
            kind=self._kind)

    def __repr__(self):
        """Returns the printable string for this object"""
        return 'Element {number}, type: {kind}'.format(number=self.number,
                                                       kind=self._kind)
