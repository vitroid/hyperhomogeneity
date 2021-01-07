"""
cycles5

Colculate the weights for cycles using SVD.
"""
from ice7analysis import *
import pickle
import sys
import random
from cycless.dicycles import dicycles_iter
import scipy.sparse.linalg

# Pseudo-inversion of a matrix.
def AtoX(A):
    Q1, S, Q2T = np.linalg.svd(A)
    S[np.abs(S)<1e-12] = 0
    rank = np.sum(np.abs(S) > 0)
    Sd = np.zeros_like(A).T
    Sd[:rank, :rank] = np.diag(1/S[:rank])
    # A†
    Ad = Q2T.T @ Sd @ Q1.T
    b = np.ones(A.shape[0])
    x = Ad @ b
    return x, rank


def cover_by_cycles(g, Nmol):
    HBcode = {(d,a): x for x, (d,a) in enumerate(g.edges())}

    A = []
    cycles = []
    for s in range(4, 20):
        lastA = len(A)
        for cycle in dicycles_iter(g,s, vec=True):
            cycles.append(cycle)
            row = np.zeros(Nmol*2)
            for j in range(len(cycle)):
                edge = HBcode[cycle[j-1], cycle[j]]
                row[edge] = 1.0
            A.append(row)
        if lastA==len(A):
            continue
        lastA = len(A)
        AT = np.array(A).T
        x, rank = AtoX(AT)
        print(f"cycle size {s}")
        print(f"number of cycles {len(cycles)}")
        b = AT @ x
        print(b)
        if np.allclose(b, np.ones_like(b)):
            # sum of weights is one at every edge
            break
    assert np.allclose(b, np.ones_like(b))
    return cycles, x

if __name__ == "__main__":
    coord = sys.argv[1]
    comeus, cell = load_nx3a(open(coord))
    cellmat = np.diag(cell)
    rcom = comeus[:,:3] @ np.linalg.inv(cellmat)
    Nmol = len(comeus)
    R = np.zeros([Nmol,3,3])
    for i in range(Nmol):
        e = comeus[i,3:6]
        R[i] = quat2rotmat(euler2quat(e))
    water = tip4picesites()
    g = hbn(rcom, cellmat, R, water)

    cycles, weights =cover_by_cycles(g, Nmol)

    if len(sys.argv) > 2:
        with open(sys.argv[2], "wb") as f:
            pickle.dump({"cycles":cycles, "weight":weights}, f)
