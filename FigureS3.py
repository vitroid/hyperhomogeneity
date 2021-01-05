"""
Interactions with cycles.
"""


from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib import pyplot as plt
from ice7analysis import *

import pickle

import networkx as nx
import random
import glob

random.seed(1)

fig = plt.figure(figsize=(5,10))
gs = plt.GridSpec(5,1)
plt.subplots_adjust(wspace=0, hspace=0)

ices = [
        ["Ih",  "1h", plt.subplot(gs[0:1, :]), 4.56],
        ["III", "3",  plt.subplot(gs[1:2, :]), 4.24],
        ["V",   "5",  plt.subplot(gs[2:3, :]), 3.96],
        ["VI",  "6",  plt.subplot(gs[3:4, :]), 4.00],
        ["VII", "7",  plt.subplot(gs[4:5, :]), 2.98],
]





for ice, basename, ax, r1 in ices:
    print(ice)
    Ymain = []
    Ysub  = []
    Yfar  = []
    Ycen  = []
    Yall  = []
    Ynear = []
    for cycles5 in glob.glob(f"{basename}-*.cycles5.pickle"):
        with open(cycles5, "rb") as f: # non-quenched is ok
            bucket = pickle.load(f)
            cycles = bucket["cycles"]
            weights = bucket["weight"]
        cyclesintr = cycles5.replace("cycles5", "cyclesintr")
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

        if ice in ("VII",):
            g = nx.Graph(dg)
            compos = [compo for compo in nx.connected_components(g)]
            assert len(compos) == 2
            assert len(compos[0]) == len(compos[1])

        for center in range(Nmol):
            sum_cen = 0.0
            sum_main = 0.0
            sum_sub  = 0.0
            sum_far  = 0.0
            sum_near = 0.0
            sum_all  = 0.0
            if ice in ("VII",):
                # if center is not in the 0th component,
                if center in compos[1]:
                    # swap them to let center be in the 0th compo.
                    compos = compos[::-1]
            for c in range(intr.shape[1]):
                cycle, weight = cycles[c], weights[c]
                e, d, morder = intr[center, c], dist[center, c], mord[center, c]

                closest = morder[0]
                if morder[0] == center:
                    closest = morder[1]

                sum_all += e*weight
                if d < 3:
                    sum_near += e*weight
                if d > 3: # VIの場合、どれが副ネットワークの第一近接かいまいちはっきりしない。
                    sum_far += e*weight
                elif morder[0] == center:
                    sum_cen += e*weight
                elif ice in ("VII",):
                    if closest in compos[0]:
                        # primary network
                        sum_main += e*weight
                    else:
                        # secondary network
                        sum_sub += e*weight
                else:
                    sum_main += e*weight
            Ymain.append(sum_main)
            Ysub.append(sum_sub)
            Ycen.append(sum_cen)
            Yfar.append(sum_far)
            Yall.append(sum_all)
            Ynear.append(sum_near)

    ax.annotate("ice {0}".format(ice), xy=(0.7,0.8), fontsize=14,xycoords='axes fraction', horizontalalignment='left')
    H1 = np.histogram(Yall, bins=100, range=(-180,20), density=True)
    ax.plot((H1[1][:-1]+H1[1][1:])/2, H1[0], label=r"$\Sigma$", color="black")
    H1 = np.histogram(Ycen, bins=100, range=(-180,20), density=True)
    ax.fill_between((H1[1][:-1]+H1[1][1:])/2, 0, H1[0], label="0", facecolor="blue", alpha=0.3)
    H1 = np.histogram(Ymain, bins=100, range=(-180,20), density=True)
    ax.fill_between((H1[1][:-1]+H1[1][1:])/2, 0, H1[0], label="1", facecolor="red", alpha=0.3)
    H1 = np.histogram(Ynear, bins=100, range=(-180,20), density=True)
    ax.fill_between((H1[1][:-1]+H1[1][1:])/2, 0, H1[0], label="0+1", facecolor="purple", alpha=0.3)
    H1 = np.histogram(Yfar, bins=100, range=(-180,20), density=True)
    ax.plot((H1[1][:-1]+H1[1][1:])/2, H1[0], label=r"2$\cdots$", color="green")
    ax.axvline(linestyle='--', color='gray')
    # if ice in ("VII",):
    #     H1 = np.histogram(Ysub, bins=100, range=(-180,20), density=True)
    #     ax.plot((H1[1][:-1]+H1[1][1:])/2, H1[0], label="1'", color="cyan")
    ax.set_xlim(-180,20)
    ax.yaxis.set_visible(False)
    if ice == "Ih":
        ax.legend(fontsize=12)
    if ice == "VII":
        ax.set_xlabel(r"Interaction Energy / kJ mol$^{-1}}$", fontsize=18)

    # dists = None
    # for seed in range(1000,1030):
    #     try:
    #         with open(f"{basename}-{seed}.couhi.pickle", "rb") as f:
    #             H2 = pickle.load(f)
    #     except:
    #         continue
    #     print(f"{basename}-{seed}.couhi.pickle")
    #     dist, x, y = H2
    #     if dists is None:
    #         dists = dist
    #     else:
    #         dists += dist
    # sliced = dists[x>r1, :][0,:]
    # ax.plot(y, sliced/ np.sum(sliced)*2, color="orange", linewidth=1) #高さてきとう
fig.savefig("FigureS3.pdf")
