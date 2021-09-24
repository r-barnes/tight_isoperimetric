#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from shapely import wkt

fig, ax = plt.subplots()

# Output file from C++ medial axis code
pts = np.loadtxt("build/voronoi_points.csv", skiprows=1, delimiter=',')

fig4 = wkt.loads(open("fig4_data.wkt").read())
for geom in fig4.geoms:
  xs, ys = geom.exterior.xy
  ax.plot(xs, ys, '-ok', lw=4)

# # Output file from C++ medial axis code
# ma = np.loadtxt("build/voronoi_edges.csv", dtype=int, skiprows=1, delimiter=',')
# for mal in ma:
#   print(mal)
#   ax.plot((pts[mal[0]][0], pts[mal[1]][0]), (pts[mal[0]][1], pts[mal[1]][1]), '-o')

# Output file from C++ medial graph code
mgp = np.loadtxt("build/medial_graph_points.csv", dtype=float, skiprows=1, delimiter=',')
sampled = mgp[mgp[:,2]>80,:]
ax.scatter(sampled[:,0], sampled[:,1], c=sampled[:,2], cmap='Greens')


# ea = np.loadtxt("build/v_edgeattributes.txt")
# for eal in ea:
#   ax.plot((eal[0], eal[2]), (eal[1], eal[3]))



plt.show()