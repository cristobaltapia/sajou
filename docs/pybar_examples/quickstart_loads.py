#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File to test the module pybar.
Example of the Book "Mechanics of Structures - Variational and computational methods"
page 284
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybar as pb
from pybar.plot import Display

m = pb.Model(name='test model', dimensionality='2D')

# add nodes
n1 = m.node(0., 0.)
n2 = m.node(0., 2000.)
n3 = m.node(2000., 2000.)
n4 = m.node(2000., 0.)

# add segment
b1 = m.beam(node1=n1, node2=n2)
b2 = m.beam(node1=n2, node2=n3)
b3 = m.beam(node1=n3, node2=n4)

# create material
mat = m.material(name='Wood', data=(200e3, ), type='isotropic')

# create beam section
section1 = m.beam_section(name='rectangular 1', material=mat, data=(
    32e2, 2356e4), type='general')
section2 = m.beam_section(name='rectangular 2', material=mat, data=(
    66e2, 5245e4), type='general')

# add beam section to the beams
b1.assign_section(section1)
b2.assign_section(section1)
b3.assign_section(section2)

# Add border conditions
m.bc(node=n1, v1=0., v2=0., r3=0.)
m.bc(node=n4, v1=0., v2=0.)

# Add load
m.load(node=n3, f2=-5e3)
m.load(node=n2, f2=-10e3)

# Distributed load
b1.distributed_load(p1=-1, p2=-2, direction='y', coord_system='local')

fig = plt.figure(figsize=(6., 5.5))
ax = fig.add_subplot(111)

disp = Display(theme='dark')

ax = disp.plot_geometry(ax=ax, model=m)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_smart_bounds(True)
ax.spines['left'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


plt.tight_layout()
