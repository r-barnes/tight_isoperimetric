#!/usr/bin/env python3

import math

import matplotlib.pyplot as plt
import numpy as np
from shapely import wkt

fig, ax = plt.subplots()

# Output file from C++ medial axis code
pts = np.loadtxt("build/voronoi_points.csv", skiprows=1, delimiter=',')


# # Output file from C++ medial axis code
# ma = np.loadtxt("build/voronoi_edges.csv", dtype=int, skiprows=1, delimiter=',')
# for mal in ma:
#   print(mal)
#   ax.plot((pts[mal[0]][0], pts[mal[1]][0]), (pts[mal[0]][1], pts[mal[1]][1]), '-o')

# Output file from C++ medial graph code
# mgp = np.loadtxt("build/medial_graph_points.csv", dtype=float, skiprows=1, delimiter=',')
# ax.scatter(mgp[:,0], mgp[:,1]) #, c=mgp[:,2], cmap='Greens')
# for i in range(len(mgp[:,0])):
#   ax.text(mgp[i,0], mgp[i,1], s="{0:.2f}".format(math.sqrt(mgp[i,2])))
#   circle = plt.Circle((mgp[i,0], mgp[i,1]), math.sqrt(mgp[i,2]))
#   ax.add_patch(circle)


# mgp = np.loadtxt("/z/circs.csv", dtype=float, skiprows=1, delimiter=',')
# ax.scatter(mgp[:,0], mgp[:,1]) #, c=mgp[:,2], cmap='Greens')
# for i in range(len(mgp[:,0])):
#   ax.text(mgp[i,0], mgp[i,1], s="{0:.2f}".format(mgp[i,2]))
#   circle = plt.Circle((mgp[i,0], mgp[i,1]), mgp[i,2])
#   ax.add_patch(circle)

# union_lines = np.loadtxt("/z/arcs.csv", dtype=float, delimiter=',')
# ax.scatter(union_lines[:,0], union_lines[:,1])


union_lines = np.loadtxt("build/poly_holes.csv", dtype=float, delimiter=',')
ax.scatter(union_lines[:,0], union_lines[:,1])

# ea = np.loadtxt("build/v_edgeattributes.txt")
# for eal in ea:
#   ax.plot((eal[0], eal[2]), (eal[1], eal[3]))



fig4 = wkt.loads(open("fig4_data.wkt").read())
for geom in fig4.geoms:
  xs, ys = geom.exterior.xy
  ax.plot(xs, ys, '-ok', lw=4)


plt.show()