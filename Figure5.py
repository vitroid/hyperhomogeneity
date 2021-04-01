"""
Divergence of the interaction at the surface of ice.

分布の表現方法を変える。ただし、そのためにはもっと大量の構造が必要。
"""

import os
import glob

from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm

from ice7analysis import *
from collections import defaultdict

distances = [3, 6, 13]
hists = defaultdict(list)

layerlabels = [(0,"top layer"), (1,"2nd"), (3,"4th"), (7,"8th")]

basename = "1cs"

coords = glob.glob(f"q/{basename}*.q.nx3a")
for coord in coords:
    if os.path.getsize(coord) == 0:
        continue
    print(coord)
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

    dg = hbn(rcom, cellmat, R, water, icetest=False)

    # 分子内座標
    waters = water @ R

    thick = cellmat[0,0]*2/32
    layer = np.floor((rcom@cellmat)[:,2] / thick + 0.5)
    # print(layer)

    for lay, layerlabel in layerlabels:
        print(f"Layer {lay}")
        bot = (lay <= layer) & (layer < lay+1)
        top = ((31-lay) <= layer) & (layer < (32-lay))
        surface = top | bot
        # print(lay, np.nonzero(surface))

        d_e = accum0(comeus, cellmat, np.nonzero(surface)[0], maxdist=13.2)
        N = len(d_e)

        for d0 in distances:
            for i, (d,e) in enumerate(d_e):
                hists[lay, d0].append(e[d<=d0][-1])

fig  = plt.figure(figsize=(5,7*3/4))
grid = plt.GridSpec(3, 1, height_ratios=[7,7,7], wspace=0, hspace=0)
linear = np.linspace(2.5,13.5,1000)

for panel, d in enumerate(distances):
    main_ax = fig.add_subplot(grid[panel,0]) #, xticks=[3,8,13] )#,
    for i, (lay, layerlabel) in enumerate(layerlabels):
        H = np.histogram(hists[lay,d], bins=100, range=(-150,-60))
        main_ax.plot((H[1][:-1]+H[1][1:])/2, H[0]/np.max(H[0]), color=cm.jet(i/3), label=layerlabel)
    if panel==0:
        main_ax.legend()
    main_ax.annotate(f"r={d*0.1:.1f} nm", xy=(0.95, 0.8), fontsize=14,xycoords='axes fraction', horizontalalignment='right')
    main_ax.set_yticks([])
    if panel==2:
        main_ax.set_xlabel(r"$I_i(r)$ / kJ mol$^{-1}$", fontsize=14)
    else:
        main_ax.set_xticks([])

    #
    # # 距離4でのエネルギーの値でソートする。
    #
    # e4 = np.zeros(N)
    # for i, (d,e) in enumerate(d_e):
    #     e4[i] = e[d<4.0][-1]
    # order = np.argsort(e4)
    # redro = np.zeros(N)
    # for i in range(N):
    #     redro[order[i]] = i
    # for i in range(N):
    #     j = redro[i]
    #     avg,sd,cnt = stepgraph(d_e[i:i+1], linear) # single data
    #     main_ax.plot(linear, avg, color=cm.coolwarm(j/N), lw=0.5)
    #
    # # center panel
    # main_ax.set_xlabel(r"$r$ / 0.1 nm",fontsize=18)
    main_ax.set_xlim(-150,-60)
    # main_ax.label_outer()
    # main_ax.set_ylim(ra)
    # main_ax.annotate(f"{layerlabel}", xy=(0.95,0.8), fontsize=18,xycoords='axes fraction', horizontalalignment='right')
    # main_ax.tick_params(labelsize=14)
    #
    # acc,sd,cnt = stepgraph(d_e, linear)
    #
    # if panel==2:
    #     main_ax.set_ylabel(r"$I_i(r)$"+r" / kJ mol$^{-1}$", fontsize=18)

fig.savefig("FigureE.pdf", bbox_inches="tight")
