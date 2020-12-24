"""
Record all the molecule-cycle interactions.
"""

from ice7analysis import *
import pickle
import networkx as nx
import random
import sys

coordfile = sys.argv[1]
cyclefile = sys.argv[2]

with open(cyclefile, "rb") as f:
    bucket = pickle.load(f)
cycles = bucket["cycles"]
weights = bucket["weight"]


comeus, cell = load_nx3a(open(coordfile))
cellmat = np.diag(cell)
rcom = comeus[:,:3] @ np.linalg.inv(cellmat)

Nmol = len(comeus)
R = np.zeros([Nmol,3,3])
for i in range(Nmol):
    e = comeus[i,3:6]
    R[i] = quat2rotmat(euler2quat(e))

water = tip4picesites()
charges = tip4picecharges()

dg = hbn(rcom, cellmat, R, water)

# 分子内座標
waters = water @ R

intr = np.zeros([Nmol, len(cycles)])
dist = np.zeros_like(intr)
mord = np.zeros([Nmol, len(cycles), 3])
for center in range(Nmol):
    print(center)
    for c, cycle in enumerate(cycles):
        # main network
        # どれが一番近い?
        e, d, morder = cycle_energy(np.array(cycle), center, rcom, cellmat, waters, dg)
        intr[center, c] = e
        dist[center, c] = d
        mord[center, c] = morder[:3]

if len(sys.argv) > 3:
    with open(sys.argv[3], "wb") as f:
        pickle.dump((intr, dist, mord), f)
