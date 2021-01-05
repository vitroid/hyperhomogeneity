"""
Interaction with a cycle does not depend on its size.
"""


from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import glob

#plt.rc(usetex = True)
from matplotlib import rc

# rc('text', usetex=True) # No use on MacPro
import pickle

import networkx as nx
import random

random.seed(1)

fig = plt.figure(figsize=(8,10))
gs = plt.GridSpec(5,2)
plt.subplots_adjust(wspace=0, hspace=0)

ices = [
        ["Ih",  "1h" ],
        ["III", "3" ],
        ["V",   "5" ],
        ["VI",  "6" ],
        ["VII", "7" ],
]


from ice7analysis import *


for panel, (ice, basename) in enumerate(ices):
    print(ice)
    ra = (-110, -40)
    nbin = ra[1] - ra[0]
    bins = [np.zeros(nbin) for i in range(20)]
    binsw = [np.zeros(nbin) for i in range(20)]
    binw = 1.0
    for cycles5 in glob.glob(f"{basename}-*.cycles5.pickle"):
        with open(cycles5, "rb") as f: # non-quenched is ok
            bucket = pickle.load(f)
            cycles = bucket["cycles"]
            weights = bucket["weight"]
        cyclesintr = cycles5.replace("cycles5", "cyclesintr")
        print(cyclesintr)
        with open(cyclesintr, "rb") as f:
            intr, dist, mord = pickle.load(f)

        coord = cycles5.replace("cycles5.pickle", "q.nx3a")
        coord = "q/" + coord

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

        # 分子内座標
        waters = water @ R

        for center in range(Nmol):
            sum_cen = defaultdict(float)
            for c in range(intr.shape[1]):
                cycle, weight = cycles[c], weights[c]
                e, d, morder = intr[center, c], dist[center, c], mord[center, c]

                if d == 0:
                    bin = int(np.floor(e+0.5)-ra[0])
                    if 0 <= bin < nbin:
                        binsw[len(cycle)][bin] += weight
                        bins[len(cycle)][bin] += 1.0 #重みをかけるのをやめ、規格化して、純粋にサイズと相互作用の関係だけに絞る。

    axL = plt.subplot(gs[panel:panel+1, :1])
    axL.annotate(f"ice {ice}",
                 xy=(0.7,0.8),
                 fontsize=14,
                 xycoords='axes fraction',
                 horizontalalignment='left')
    for nmem in range(4,17):
        total = np.sum(bins[nmem])
        if total > 0:
            print(nmem)
            axL.plot(range(*ra), bins[nmem] / (total*binw), label=f"{nmem}-member", color=cm.hsv((nmem-4)/(17-4)))
    axL.set_xlim(*ra)
    axL.yaxis.set_visible(False)
    if ice == "VII":
        axL.set_xlabel(r"Coulomb Interaction / kJ mol$^{-1}}$", fontsize=14)
    else:
        axL.xaxis.set_visible(False)

    axR = plt.subplot(gs[panel:panel+1, -1:])
    # fig = plt.figure() bbox_to_anchor=(1.05, 1), loc='upper left',
    axR.annotate("ice {0}".format(ice), xy=(0.7,0.8), fontsize=14,xycoords='axes fraction', horizontalalignment='left')
    for nmem in range(4,20):
        total = np.sum(bins[nmem])
        if total > 0:
            print(nmem)
            axR.plot(range(*ra),
                     binsw[nmem],
                     label=f"{nmem}-member",
                     color=cm.hsv((nmem-4)/(17-4)))
    axR.set_xlim(*ra)
    axR.yaxis.set_visible(False)
    if ice == "V":
        axR.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
    if ice == "VII":
        axR.set_xlabel(r"Coulomb Interaction / kJ mol$^{-1}}$", fontsize=14)
    else:
        axR.xaxis.set_visible(False)
plt.show()
fig.savefig("FigureS4.pdf", bbox_inches="tight")
