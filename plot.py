#!/usr/bin/env python3

import math

import matplotlib.pyplot as plt
import numpy as np
from shapely import wkt

SHOW_FIGURE_BEING_ANALYZED = True
SHOW_MEDIAL_GRAPH_SAMPLE_POINTS = True
SHOW_MEDIAL_GRAPH_VERTICES_AND_LINE_SEGMENTS = True

fig, ax = plt.subplots()

# Output file from C++ medial axis code
pts = np.loadtxt("voronoi_points.csv", skiprows=1, delimiter=',')


if SHOW_MEDIAL_GRAPH_VERTICES_AND_LINE_SEGMENTS:
  ma = np.loadtxt("voronoi_edges.csv", dtype=int, skiprows=1, delimiter=',')
  for mal in ma:
    print(mal)
    ax.plot((pts[mal[0]][0], pts[mal[1]][0]), (pts[mal[0]][1], pts[mal[1]][1]), '-o')


if SHOW_MEDIAL_GRAPH_SAMPLE_POINTS:
  mgp = np.loadtxt("medial_graph_points.csv", dtype=float, skiprows=1, delimiter=',')
  ax.scatter(mgp[:,0], mgp[:,1], c=mgp[:,2], cmap='Greens')
  # The following is a sanity check that shows the numeric distances of each
  # sample point from the edge of the figure and draws a circle around each
  # point. If the circles are all in the figure and cover it entirely, then
  # things are probably good. TODO(r-barnes): Make this pretty and uncommented
  # for i in range(len(mgp[:,0])):
  #   ax.text(mgp[i,0], mgp[i,1], s="{0:.2f}".format(math.sqrt(mgp[i,2])))
  #   circle = plt.Circle((mgp[i,0], mgp[i,1]), math.sqrt(mgp[i,2]))
  #   ax.add_patch(circle)


if SHOW_FIGURE_BEING_ANALYZED:
  fig4 = wkt.loads(open("data/fig4_data.wkt").read())
  for geom in fig4.geoms:
    xs, ys = geom.exterior.xy
    ax.plot(xs, ys, '-ok', lw=4)


plt.show()