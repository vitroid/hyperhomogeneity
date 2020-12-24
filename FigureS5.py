from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm
import glob
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
from ice7analysis import *

# fig, ax = plt.subplots(1,3)
fig  = plt.figure(figsize=(5,14))
linear = np.linspace(2.5,13.5,1000)

import pickle


ices = [
        ["Ih", "1h"],
        ["Ic", "1c"],
        ["III", "3"],
        ["V", "5"],
        ["VI", "6"],
        ["VII", "7"],
        ["XVI", "16"],
        ["empty sI", "CS1"],
        ["2D", "2D2"],
]
grid = plt.GridSpec(len(ices), 3, height_ratios=[7]*len(ices), width_ratios=[1, 3, 1], wspace=0, hspace=0)


linear = np.linspace(0.0,20,1000)

for panel, (ice, num) in enumerate(ices):
    ra = [-180, -120]
    main_ax = fig.add_subplot(grid[panel,1], xticks=[0,6,12])
    hist0   = fig.add_subplot(grid[panel,0], xticklabels=[], sharey=main_ax, yticks=(-180,-160,-140))
    hist13  = fig.add_subplot(grid[panel,2], xticklabels=[], sharey=main_ax)

    # 配向ごとに別個のd_eを準備する。
    d_e = defaultdict(list)

    for cycle5stat in glob.glob(f"{num}-*.cycles5stat.pickle"):
        with open(cycle5stat, "rb") as f:
            d_e_elem = pickle.load(f)
        for ori in d_e_elem:
            d_e[ori] += d_e_elem[ori]

    print(ice, len(d_e))
    for ori in d_e:
        acc,sd,cnt = stepgraph(d_e[ori], linear)
        bandplot1(linear,acc,sd, plt=main_ax, label=f"{ori}")
        # print(ori, sd[1], sd[-1])
        e13 = []
        for d,e in d_e[ori]:
            e13.append(e[d<13][-1])
        H = np.histogram(e13, bins=40, range=(-180, -100))
        hist13.plot(H[0], (H[1][1:]+H[1][:-1])/2, label=f"{ori}") #normed=True, histtype="step",

        e0 = []
        for d,e in d_e[ori]:
            e0.append(e[d>0][0])
        H = np.histogram(e0, bins=40, range=(-180, -100))
        hist0.plot(H[0], (H[1][1:]+H[1][:-1])/2, label=f"{ori}") #normed=True, histtype="step",

    main_ax.set_xlabel("Distance / 0.1 nm",fontsize=18)
    main_ax.set_xlim(0,13)
    main_ax.label_outer()
    main_ax.set_ylim(ra)
    if ice == "2D":
        name = f"{ice} ice"
    elif num == "CS1":
        name = f"{ice}"
    else:
        name = f"ice {ice}"
    main_ax.annotate(name, xy=(0.95,0.8), fontsize=18,xycoords='axes fraction', horizontalalignment='right')
    main_ax.tick_params(labelsize=14)

    if num == "3":
        d_e = defaultdict(list)
        with open("9.cycles5stat.pickle", "rb") as f:
            d_e_elem = pickle.load(f)
        for ori in d_e_elem:
            d_e[ori] += d_e_elem[ori]
        for ori in d_e:
            acc,sd,cnt = stepgraph(d_e[ori], linear)
            main_ax.plot(linear, acc, "k-")
#             bandplot1(linear,acc,sd, plt=main_ax, label=f"{ori}")
#             print(ori, sd[1], sd[-1])
#             e13 = []
#             for d,e in d_e[ori]:
#                 e13.append(e[d<13][-1])
#             H = np.histogram(e13, bins=40, range=(-180, -100))
#             hist13.plot(H[0], (H[1][1:]+H[1][:-1])/2, label=f"{ori}") #normed=True, histtype="step",

#             e0 = []
#             for d,e in d_e[ori]:
#                 e0.append(e[d>0][0])
#             H = np.histogram(e0, bins=40, range=(-180, -100))
#             hist0.plot(H[0], (H[1][1:]+H[1][:-1])/2, label=f"{ori}") #normed=True, histtype="step",


    # left
    if panel==len(ices)//2:
        hist0.set_ylabel(r"$I_i^c(r)$" + r" / kJ mol$-1$", fontsize=18)
    #    hist0.fill(HH[0,:], ygauge, color=cm.viridis(panel/4))
    #     ixmax = np.argmin(np.abs(xgauge-xmax))
    #     hist0.plot(HH[ixmax,:], ygauge, 'k-', lw=0.5)

    hist0.invert_xaxis()
    hist0.set_xlim(None,0)
    hist0.set_ylim(ra)
    hist0.tick_params(labelsize=14)

    # right
    #     hist13.fill(HH[-1,:], ygauge, color=cm.viridis(panel/4))
    hist13.set_xlim(0,None)
    hist13.label_outer()
    hist13.set_ylim(ra)

plt.show()
fig.savefig("FigureS5.pdf", bbox_inches="tight")
