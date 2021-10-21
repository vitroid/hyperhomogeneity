#!/usr/bin/env python3

# double the cell

import numpy as np
import sys

def doubl(nx3a, box3, nx=2,ny=2,nz=2):
    nx3a[:, :3] /= box3
    vs = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                a = nx3a.copy()
                a[:,:3] += np.array([ix,iy,iz])
                vs.append(a)
    dbl = np.vstack(vs)
    dbl[:,:3] *= box3
    return dbl, box3*np.array([nx,ny,nz])


def load(fh):
    while True:
        line = fh.readline()
        if line[:5] == "@BOX3":
            line = fh.readline()
            box3 = np.array([float(x) for x in line.split()])
        elif line[:5] == "@NX3A":
            line = fh.readline()
            N = int(line)
            nx3a = []
            for i in range(N):
                line = fh.readline()
                nx3a.append([float(x) for x in line.split()])
            nx3a = np.array(nx3a)
            return nx3a, box3


nx3a, box3 = load(sys.stdin)
# print(nx3a)
n, b = doubl(nx3a, box3)
# print(n)
print("@BOX3")
[print(x, end=" ") for x in b]
print()
print("@NX3A")
print(n.shape[0])
for d in n:
    [print(x, end=" ") for x in d]
    print()
