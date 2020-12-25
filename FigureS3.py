"""
Illustration of the cycles that pass through a water molecule in ice Ih.
"""

import yaplotlib as yap

s = ""
s += yap.SetPalette(1,255,255,255) #black bg
for i in range(30): # -1 .. 0  blue to white
    s += yap.SetPalette(i+10, i*255//30, i*255//30, 255)
for i in range(30): # 0 .. 1  white to black
    s += yap.SetPalette(i+40, (30-i)*255//30, (30-i)*255//30, (30-i)*255//30)
for i in range(30): # 1 .. 2  black to red
    s += yap.SetPalette(i+70, i*255//30, 0, 0)
s += yap.SetPalette(5,150,0,0) #black bg
s += yap.SetPalette(6,0,150,150) #black bg


import pickle

with open("1h-1000.cycles5.pickle", "rb") as f:
    bucket = pickle.load(f)
cycles = bucket["cycles"]
weights = bucket["weight"]

from ice7analysis import *

coord = "q/1h-1000.q.nx3a"
comeus, cell = load_nx3a(open(coord))
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



cpos = np.array([0.5,0.5,0.5])
center = np.argsort(np.linalg.norm(cpos-rcom, axis=1))[0]

w = water @ R[center] + rcom[center] @ cellmat
s += yap.Size(0.5)
s += yap.Color(5) #red
s += yap.Circle(w[3])

s += yap.Size(0.25)
s += yap.Color(6) #cyan?
s += yap.Circle(w[0])
s += yap.Circle(w[1])

sumw_center = defaultdict(float)

for cycle, weight in zip(cycles, weights):
    if center in cycle:
        for i in range(len(cycle)):
            a, b = cycle[i-1], cycle[i]
            sumw_center[a,b] += weight

for a,b in sumw_center:
    weight = sumw_center[a,b]
    palette = int((weight+1)*30)
    if palette < 0:
        palette = 0
    elif palette >= 90:
        palette = 89
    s += yap.Color(palette + 10)
    pa = rcom[a]
    pb = rcom[b]
    d = pb - pa
    d -= np.floor(d+0.5)
    s += yap.Line(pa@cellmat, (pa+d)@cellmat)

with open("FigureS2.yap", "w") as f:
    f.write(s)
