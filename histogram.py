from ice7analysis import *
import glob
import sys

linear = np.linspace(3, 13, 501)
values = np.linspace(-160,-20,281)
# print(linear)
H2 = dict()
infile = sys.argv[1]
outfile = sys.argv[2]
HH = np.zeros([linear.shape[0], values.shape[0]])
comeus, cell = load_nx3a(open(infile))
cellmat = np.diag(cell)
d_e = accum0(comeus, cellmat, range(len(comeus)), maxdist=13.2)
for sf in steppify(d_e, linear):
    isf = ((sf+160)*2+0.5).astype(int)
    for i, v in enumerate(isf):
        if v < values.shape[0]:
            HH[i,v] += 1
H2 = HH, linear, values


import pickle
with open(outfile, 'wb') as f:
    pickle.dump(H2, f)
