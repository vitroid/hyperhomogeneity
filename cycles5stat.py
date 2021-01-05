import pickle
import sys
from ice7analysis import *
from collections import defaultdict



def cycle_energy(cycle, center, rcom, cellmat, waters, dg):
    # 小さめに限定
    # assert len(cycle) < 10

    # 酸素の座標は一括抽出できる。
    # 酸素の電荷はO位置ではなくM位置(3番目)
    oxygens = waters[cycle, 3]

    # cycleの中にcenterを含む場合はoxygensとhydrogensからそれらを除去しておく必要がある。
    if center in cycle:
        oxygens = waters[cycle[cycle!=center], 3]

    # 水素はどっちか判別必要。
    which = []
    N = len(cycle)
    #print(cycle)
    for i, w in enumerate(cycle):
        nxt = cycle[(i+1) % N]
        which.append(dg[w][nxt]["H"])
    hydrogens = waters[cycle, which]

    if center in cycle:
        which = np.array(which)
        hydrogens = waters[cycle[cycle!=center], which[cycle!=center]]

    # 分子の相対位置
    d = rcom[cycle] - rcom[center]
    d -= np.floor(d + 0.5)
    d = d @ cellmat

    if center in cycle:
        d = d[cycle!=center]
    dmin = np.min(np.linalg.norm(d, axis=1))

    if center in cycle:
        dmin = 0.0
    # print(np.linalg.norm(d, axis=1))
    # coulomb相互作用の和
    ec = 0
    charges = [1.0, 1.0, 0.0, -2.0] # 中心分子の電荷。中心分子は半分ではない。
    for i in (0,1,3):
        # 中心分子の原子位置
        ac = waters[center, i]

        # cycle上の半酸素位置との差
        rel = oxygens - ac + d
        L   = np.linalg.norm(rel, axis=1)
        ec -= np.sum(1/L) * charges[i]

        # cycle上の水素位置との差
        rel = hydrogens - ac + d
        L   = np.linalg.norm(rel, axis=1)
        ec += np.sum(1/L) * charges[i]

    unitcharge = tip4picecharges()[0]
    return ec * unitcharge**2 * qe0, dmin



# read repr.pickle
with open(sys.argv[2], "rb") as f:
    repr = pickle.load(f)

# distance-energy pairs for different arrangements
d_e = defaultdict(list)

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
charges = tip4picecharges()

dg = hbn(rcom, cellmat, R, water)

# Intramolecular coordinates
waters = water @ R

# .cycles5.pickle; Cycles and weights
bucket = pickle.load(open(sys.argv[3], "rb"))
cycles_raw = bucket["cycles"]
weight     = bucket["weight"]


for center in range(0,Nmol,1):
    print(center)
    n1, n2 = [x for x in dg.successors(center)]
    if (center, n1, n2) in repr:
        ori = repr[center, n1, n2][0]
    else:
        ori = repr[center, n2, n1][0]
    # cycle by cycle
    eps = []
    ds  = []
    for cycle, wei in zip(cycles_raw, weight):
        ep, d = cycle_energy(np.array(cycle), center, rcom, cellmat, waters, dg)
        eps.append(ep*wei)
        ds.append(d)

    ds = np.array(ds)
    order = np.argsort(ds)
    eps = np.array(eps)
    d_e[ori].append([ds[order],np.cumsum(eps[order])])

if len(sys.argv) > 4:
    with open(sys.argv[4], "wb") as f:
        pickle.dump(d_e, f)
