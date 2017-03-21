#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains an extension to the MarkerStyle class of Matplotlib.

It implements an additional set of markers which are used for the
representation of different attributes of a structural analysis model.

In addition to the standard markers, the following are implemented:

    +---------+----------------------------------+
    | marker  | description                      |
    +=========+==================================+
    | `'ap'`  | arc with arrow (anti clock-wise) |
    +---------+----------------------------------+
    | `'an'`  | arc with arrow (clock-wise)      |
    +---------+----------------------------------+
    | `'psx'` | horizontal pinned support        |
    +---------+----------------------------------+
    | `'psy'` | vertical pinned support          |
    +---------+----------------------------------+
    | `'rsx'` | horizontal rolling support       |
    +---------+----------------------------------+
    | `'rsy'` | vertical rolling support         |
    +---------+----------------------------------+
    | `'es'`  | encastrated support              |
    +---------+----------------------------------+
    | `'rex'` | 'rolling_encastrated_x',         |
    +---------+----------------------------------+
    | `'rey'` | 'rolling_encastrated_y',         |
    +---------+----------------------------------+

Examples
--------

.. plot:: sajou_examples/markers.py

.. literalinclude:: ../../docs/sajou_examples/markers.py

"""
import math
import numpy as np

from matplotlib import markers as mpl_mk
from matplotlib.path import Path
from matplotlib.transforms import IdentityTransform, Affine2D


class MarkerStyle(mpl_mk.MarkerStyle):
    """Defines custom marker styles for Sajou

    Attributes
    ----------
    markers : list of known markes
    fillstyles : list of known fillstyles
    filled_markers : list of known filled markers.

    Parameters
    ----------
    marker : string or array_like, optional, default: None
        See the descriptions of possible markers in the module docstring.
    fillstyle : string, optional, default: 'full'
        'full', 'left", 'right', 'bottom', 'top', 'none'

    """

    markers_custom = {
        'ap': 'arc_arrow_positive',
        'an': 'arc_arrow_negative',
        'psx': 'pinned_support_x',
        'psy': 'pinned_support_y',
        'rsx': 'rolling_support_x',
        'rsy': 'rolling_support_y',
        'es': 'encastrated_support',
        'rex': 'rolling_encastrated_x',
        'rey': 'rolling_encastrated_y',
    }

    markers = {**markers_custom, **mpl_mk.MarkerStyle.markers}

    def __init__(self, marker=None, fillstyle=None):
        mpl_mk.MarkerStyle.__init__(self, marker=None, fillstyle=None)
        self.set_fillstyle(fillstyle)
        self.set_marker(marker)

    def _set_arc_arrow(self, sign):
        self._transform = Affine2D().translate(0, 0)
        # Define the arc between -45° and 135°
        path_arc = Path.arc(theta1=-45., theta2=135.)
        # get the vertices and codes
        verts_m = path_arc.vertices
        codes_m = path_arc.codes
        # size of the arrow
        arrow_head_s = 0.4

        # get position of the tip
        if sign == 'positive':
            pos_tip = np.arctan2((verts_m[-1, 1] - verts_m[-3, 1]),
                                 (verts_m[-1, 0] - verts_m[-3, 0]))
            tip_arrow = verts_m[-1, :]

        elif sign == 'negative':
            pos_tip = np.arctan2((verts_m[0, 1] - verts_m[2, 1]),
                                 (verts_m[0, 0] - verts_m[2, 0]))
            tip_arrow = verts_m[0, :]

        # define the angle of the tip arrow
        ang_arrow = np.deg2rad(25.)
        ang_p1 = pos_tip - np.pi + ang_arrow
        ang_p2 = pos_tip - np.pi - ang_arrow

        # Position of the first line of the arrow
        p1 = tip_arrow + arrow_head_s * np.array(
            [np.cos(ang_p1), np.sin(ang_p1)])
        # Position of the second line of the arrow
        p2 = tip_arrow + arrow_head_s * np.array(
            [np.cos(ang_p2), np.sin(ang_p2)])
        # Add to the vertices and code
        aux = np.array([p1, tip_arrow, p2])
        verts_mp = np.vstack((verts_m, aux))
        aux_code = np.array([1, 2, 2])
        codes_mp = np.hstack((codes_m, aux_code))

        # Create path
        self._path = Path(verts_mp, codes_mp)
        self._joinstyle = 'miter'
        self._transform.scale(3.0, 3.0)

    def _set_arc_arrow_positive(self):
        return self._set_arc_arrow('positive')

    def _set_arc_arrow_negative(self):
        return self._set_arc_arrow('negative')

    def _set_pinned_support(self, rot):
        self._transform = Affine2D().translate(0, 0).rotate_deg(rot)
        x = [
            0, 1, -1, 0, -1.2, 1.2, -1.2, -0.8, -0.8, -0.4, -0.4, 0, 0, 0.4,
            0.4, 0.8, 0.8, 1.2
        ]
        y = [
            0, -1, -1, 0, -1, -1, -1.5, -1, -1.5, -1, -1.5, -1, -1.5, -1, -1.5,
            -1, -1.5, -1
        ]
        xy = list(zip(x, y))

        codes = [
            Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.MOVETO,
            Path.LINETO, Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
            Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO, Path.MOVETO,
            Path.LINETO, Path.MOVETO, Path.LINETO
        ]
        self._filled = False
        self._path = Path(xy, codes)
        self._transform.scale(2.0, 2.0)
        self._joinstyle = 'miter'

    def _set_pinned_support_x(self):
        return self._set_pinned_support(0.)

    def _set_pinned_support_y(self):
        return self._set_pinned_support(90.)

    def _set_rolling_support(self, rot):
        self._transform = Affine2D().translate(0, 0).rotate_deg(rot)

        x = [0, 1, -1, 0, -1.2, 1.2]
        y = [0, -1, -1, 0, -1.5, -1.5]
        xy = list(zip(x, y))
        codes = [
            Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.MOVETO,
            Path.LINETO
        ]

        self._filled = False
        self._path = Path(xy, codes)
        self._transform.scale(2.0, 2.0)
        self._joinstyle = 'miter'

    def _set_rolling_support_x(self):
        return self._set_rolling_support(0)

    def _set_rolling_support_y(self):
        return self._set_rolling_support(90)

    def _set_encastrated_support(self):
        self._transform = Affine2D().translate(0, 0)
        x = [
            -1, 1, 1, -1, -1, -1.2, -0.8, -0.8, -0.4, -0.4, 0, 0, 0.4, 0.4,
            0.8, 0.8, 1.2
        ]
        y = [
            0, 0, -1, -1, 0, -1.5, -1, -1.5, -1, -1.5, -1, -1.5, -1, -1.5, -1,
            -1.5, -1
        ]
        xy = list(zip(x, y))
        codes = [
            Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
            Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO, Path.MOVETO,
            Path.LINETO, Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO,
            Path.MOVETO, Path.LINETO
        ]

        self._filled = False
        self._path = Path(xy, codes)
        self._transform.scale(2.0, 2.0)
        self._joinstyle = 'miter'

    def _set_rolling_encastrated(self, rot):
        self._transform = Affine2D().translate(0, 0).rotate_deg(rot)
        x = [-1, 1, 1, -1, -1, -1.2, 1.2]
        y = [0, 0, -1, -1, 0, -1.5, -1.5]
        xy = list(zip(x, y))
        codes = [
            Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
            Path.MOVETO, Path.LINETO
        ]

        self._filled = False
        self._path = Path(xy, codes)
        self._transform.scale(2.0, 2.0)
        self._joinstyle = 'miter'

    def _set_rolling_encastrated_x(self):
        return self._set_rolling_encastrated(0)

    def _set_rolling_encastrated_y(self):
        return self._set_rolling_encastrated(90)
