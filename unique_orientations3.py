"""
Find the crystallographically identical molecular arrangements in the ice lattice.
"""

import sys
from ice7analysis import *

import pairlist as pl
import networkx as nx
import itertools as it

basename=sys.argv[1]
reprfile=sys.argv[2]

comeus, cell = load_nx3a(open(basename))
Rcellmat = np.diag(cell)
Rrpos = comeus[:,:3] @ np.linalg.inv(Rcellmat)
Nmol = len(comeus)
RR = np.zeros([Nmol,3,3])
for i in range(Nmol):
    e = comeus[i,3:6]
    RR[i] = quat2rotmat(euler2quat(e))


water = tip4picesites()
g = nx.Graph(hbn(Rrpos, Rcellmat, RR, water))


# fortran match function
from match import match

def match_pointcloud(p1, p2, thres=0.8):
    def _match(qS,qL,thres2):
        err = 0
        for i,x in enumerate(qS):
            found = False
            for j,y in enumerate(qL):
                d = x-y
                if d@d < thres2:
                    found = True
                    err += d@d
                    break
            if not found:
                return -1
        return err

    assert p1.shape[0] <= p2.shape[0]
    # Mirror the template.0=None 1=x 2=xy 3=y
    bestm = -1
    beste = 9999
    err0 =  match(p1, p2, thres**2)
    if err0 > 0 and err0 < beste:
        bestm = 0
        beste = err0
    p1x = p1.copy()
    p1x[:,0] = -p1[:,0]
    err1 = match(p1x,p2, thres**2)
    if err1 > 0 and err1 < beste:
        bestm = 1
        beste = err1
    p1x[:,1] = -p1[:,1]
    err2 = match(p1x,p2, thres**2)
    if err2 > 0 and err2 < beste:
        bestm = 2
        beste = err2
    p1x[:,0] = p1[:,0]
    err3 = match(p1x,p2, thres**2)
    if err3 > 0 and err3 < beste:
        bestm = 3
        beste = err3
    return bestm, beste



oris = dict()
repr = dict()
for i, c in enumerate(g):
    for n1,n2 in it.combinations(g[c], 2):
        reloc = Rrpos - Rrpos[c]
        reloc -= np.floor(reloc+0.5)
        Rc = nodes2rotmat(c,n1,n2,Rrpos,Rcellmat)
        pos = reloc @ Rcellmat @ Rc.T
        found = False
        best = 9999
        bestl = None
        bestm = -1
        for label in oris:
            mirror, score = match_pointcloud(oris[label], pos)
            if mirror >= 0:
                if best > score:
                    best = score
                    bestl = label
                    bestm = mirror
        print(best, "{0}/{1}".format(i, len(Rrpos)))
        if best != 9999:
            repr[c,n1,n2] = bestl, bestm
        else:
            newlabel = c,n1,n2
            oris[newlabel] = pos[np.linalg.norm(pos, axis=1)<6.0]
            repr[newlabel] = newlabel, 0
            print(newlabel) #, oris[label])

sss = dict()
for i in repr:
    ri = repr[i][0][0]
    if ri not in sss:
        sss[ri] = 0
    sss[ri] +=1
sss = dict()
for i in repr:
    ri = repr[i][0]
    if ri not in sss:
        sss[ri] = 0
    sss[ri] +=1
print(len(sss))

import pickle
with open(reprfile, "wb") as f:
    pickle.dump(repr, f)
