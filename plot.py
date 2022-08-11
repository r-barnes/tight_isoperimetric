#!/usr/bin/env python3

import glob
import math
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from shapely import wkt

SHOW_FIGURE_BEING_ANALYZED = True
SHOW_MEDIAL_GRAPH_SAMPLE_POINTS = False
SHOW_MEDIAL_GRAPH_VERTICES_AND_LINE_SEGMENTS = False
SHOW_POLY_HOLES = False

def read_split_table_file(filename: str) -> List[np.ndarray]:
  """Reads a file with the format
      x y
      x y
      x y
      #
      x y
      x y
      x y
  Into a list of numpy arrays
  """
  with open(filename, "r") as fin:
    arrays = []
    array = []
    for line in fin.readlines():
      if line.startswith("#"):
        arrays.append(np.array(array))
        array = []
        continue
      array.append([float(x) for x in line.split()])
  return arrays


fig, ax = plt.subplots()
ax.axis('equal')

# Output file from C++ medial axis code
if SHOW_MEDIAL_GRAPH_VERTICES_AND_LINE_SEGMENTS:
  pts = np.loadtxt("voronoi_points.csv", skiprows=1, delimiter=',')
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
  for i in range(len(mgp[:,0])):
    ax.text(mgp[i,0], mgp[i,1], s="{0:.2f}".format(math.sqrt(mgp[i,2])))
    circle = plt.Circle((mgp[i,0], mgp[i,1]), math.sqrt(mgp[i,2]), alpha=0.2)
    ax.add_patch(circle)


if SHOW_FIGURE_BEING_ANALYZED:
  fig4 = wkt.loads(open("data/fig4_data.wkt").read())
  for geom in fig4.geoms:
    xs, ys = geom.exterior.xy
    ax.plot(xs, ys, '-ok', lw=4)


plt.show()


if SHOW_POLY_HOLES:
  colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
  for step in range(255):
    poly_filename = f"polygon_w_holes_output_step{step}_*.out"
    print(poly_filename)
    fig, ax = plt.subplots()
    ax.axis('equal')
    for i, filename in enumerate(glob.glob(poly_filename)):
      for x in read_split_table_file(filename):
        ax.plot(x[:,0], x[:,1], colors[i])
    plt.show()


# if SHOW_POLY_HOLES:
#   with open("poly_holes.csv", "r") as fin:
#     holes_data = fin.readlines()
#   holes_data = [x.strip() for x in holes_data]
#   # Identify shape by string prefix on line and remove said prefix
#   circles = [x[7:] for x in holes_data if x.startswith("circle")]
#   lines = [x[6:] for x in holes_data if x.startswith("line")]
#   # Split up coordinates
#   circles = [[float(c) for c in x.split(',')] for x in circles]
#   lines = [[float(c) for c in x.split(',')] for x in lines]
#   print(circles)
#   print(lines)
#   for x, y, sr in circles:
#     print(x,y,sr)
#     circle = plt.Circle((x, y), math.sqrt(sr), alpha=0.01)
#     ax.add_patch(circle)
