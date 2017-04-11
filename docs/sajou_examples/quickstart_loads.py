#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File to test the module sajou.
Example of the Book "Mechanics of Structures - Variational and computational methods"
page 284
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sajou as sj

m = sj.Model(name='test model', dimensionality='2D')

# add nodes
n1 = m.node(0., 0.)
n2 = m.node(0., 2000.)
n3 = m.node(1500., 2500.)
n4 = m.node(3000., 2000.)
n5 = m.node(3000., 0.)

# add segment
b1 = m.beam(node1=n1, node2=n2)
b2 = m.beam(node1=n2, node2=n3)
b3 = m.beam(node1=n3, node2=n4)
b4 = m.beam(node1=n4, node2=n5)

# create material
mat = m.material(name='Wood', data=(12e3, ), type='isotropic')

# create beam section
section1 = m.beam_section(name='rectangular 1', material=mat, data=(
    300, 100), type='rectangular')

# add beam section to the beams
b1.assign_section(section1)
b2.assign_section(section1)
b3.assign_section(section1)

# Add border conditions
m.bc(node=n1, v1=0., v2=0.)
m.bc(node=n5, v1=0., v2=0.)

# Add load
m.load(node=n3, f2=-10e3)

# Distributed load
b1.distributed_load(p1=-1, p2=-2, direction='y', coord_system='local')
b2.distributed_load(p1=-1, direction='y', coord_system='global')
b3.distributed_load(p1=-1, direction='y', coord_system='global')

# release end
b2.release_end(which=2)

fig = plt.figure(figsize=(6., 5.5))
ax = fig.add_subplot(111)
disp = sj.Display(theme='dark')
ax = disp.plot_geometry(ax=ax, model=m)

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_smart_bounds(True)
ax.spines['left'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlim(xmin=-800, xmax=3500)

plt.tight_layout()
plt.show()
