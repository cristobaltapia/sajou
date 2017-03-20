#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extends the Line2D class of matplotlib to use the custom markers
"""

from matplotlib import lines
from .markers_mpl import MarkerStyle


class Line2D(lines.Line2D):
    """Extends the Line2D class from matplotlib."""
    markers = MarkerStyle.markers

    def __init__(self,
                 xdata,
                 ydata,
                 linewidth=None,  # all Nones default to rc
                 linestyle=None,
                 color=None,
                 marker=None,
                 markersize=None,
                 markeredgewidth=None,
                 markeredgecolor=None,
                 markerfacecolor=None,
                 markerfacecoloralt='none',
                 fillstyle=None,
                 antialiased=None,
                 dash_capstyle=None,
                 solid_capstyle=None,
                 dash_joinstyle=None,
                 solid_joinstyle=None,
                 pickradius=5,
                 drawstyle=None,
                 markevery=None,
                 **kwargs):

        lines.Line2D.__init__(self,
                              xdata=xdata,
                              ydata=ydata,
                              linewidth=linewidth,  # all Nones default to rc
                              linestyle=linestyle,
                              color=color,
                              markersize=markersize,
                              markeredgewidth=markeredgewidth,
                              markeredgecolor=markeredgecolor,
                              markerfacecolor=markerfacecolor,
                              markerfacecoloralt=markerfacecoloralt,
                              fillstyle=fillstyle,
                              antialiased=antialiased,
                              dash_capstyle=dash_capstyle,
                              solid_capstyle=solid_capstyle,
                              dash_joinstyle=dash_joinstyle,
                              solid_joinstyle=solid_joinstyle,
                              pickradius=pickradius,
                              drawstyle=drawstyle,
                              markevery=markevery,
                              **kwargs)

        self._marker = MarkerStyle(marker, fillstyle)


