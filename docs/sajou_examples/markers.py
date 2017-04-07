#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Example to show the new marker styles"""
import matplotlib.pyplot as plt
from sajou.plot.lines_mpl import Line2D

fig = plt.figure(figsize=(12, 3))
ax = fig.add_subplot(111)
markers = ['ap', 'an', 'psx', 'rsx', 'es', 'rex', 'rc']

for ix, mark in enumerate(markers):
    marker = Line2D([ix], [0], marker=mark, fillstyle='none', color='k')
    ax.add_line(marker)

ax.set_xlim(-1, len(markers))
ax.set_ylim(-1, 1)
plt.show()
