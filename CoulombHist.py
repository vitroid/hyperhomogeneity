"""
Record the cumulative Coulomb interaction on a molecule-basis.
"""

from ice7analysis import *
import pickle
import sys


linear = np.linspace(3, 13, 501)
ra = (-200,-40)
binw = 0.5
nbin = int((ra[1] -ra[0])/binw)
values = np.linspace(ra[0], ra[1], nbin+1)
# print(linear)

coord = sys.argv[1]
outfile = sys.argv[2]

HH = np.zeros([linear.shape[0], values.shape[0]])
comeus, cell = load_nx3a(open(coord))
print(coord)
cellmat = np.diag(cell)
d_e = accum0(comeus, cellmat, range(len(comeus)), maxdist=13.2, LJ=(0,0))
for sf in steppify(d_e, linear):
    isf = ((sf-ra[0])/binw+0.5).astype(int)
    assert np.all(0 <= isf), np.min(sf)
    assert np.all(isf < nbin), np.max(sf)
    for i, v in enumerate(isf):
        if v < values.shape[0]:
            HH[i,v] += 1
H2 = HH, linear, values
with open(outfile, 'wb') as f:
    pickle.dump(H2, f)
