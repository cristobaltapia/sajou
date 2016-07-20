#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib as pyplot
"""
This fil defines all the functions needed to plot the created structures and the results
obtained with the program.
"""

def plot_geometry(model, ax):
    """Plots the geometry of the model passed

    :model: TODO
    :ax: a matplotlib axis object
    :returns: TODO

    """
    nodes = model.nodes
