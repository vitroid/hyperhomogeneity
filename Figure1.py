"""
For resubmission

Figure 1: Interaction distribution with nearby molecules.
"""

from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm
from ice7analysis import *
import glob

# fig, ax = plt.subplots(1,3)
fig  = plt.figure(figsize=(5,7))
# fig  = plt.figure()
grid = plt.GridSpec(4, 1, height_ratios=[7,7,7,7], wspace=0, hspace=0)
linear = np.linspace(2.5,13.5,1000)

import pickle



conv = dict()
for line in """
# ice c xmax ymax ylast
1h 3.6 5.5 7.4 2.1
3 3.5 4.3 12.3 3.5
5 3.1 4.0 13.3 4.3
6 3.9 4.4 13.7 3.5
7 7.4 3.1 19.5 2.6
""".splitlines():
    if len(line) == 0:
        continue
    if line[0] == "#":
        continue
    cols = line.split()
    conv[cols[0]] = [float(x) for x in cols[1:]]

panels = [0]*4
for i in range(4):
    panels[i] = fig.add_subplot(grid[i,0]) #, xticks=[3,8,13])

ices = ["1h", "3", "5", "6", "7"]
icef = ["1h", "3", "5m", "6m", "7"]
ices = [x for x in ices[::-1]]
icef = [x for x in icef[::-1]]
for num, ice in enumerate(ices):
    rel = num / (len(ices)-1)
    comeus, cell = load_nx3a(open(f"q/{ice}-1000.q.nx3a"))
    cellmat = np.diag(cell)
    d_e = accum0(comeus, cellmat, range(50), maxdist=13.2)
    N = len(d_e)
    ra = [-150,-30]

    HH = None
    for file in glob.glob(f"{ice}-*.hist.pickle"):
        with open(file, "rb") as f:
            print(file)
            H2 = pickle.load(f)
            if HH is None:
                HH = np.zeros_like(H2[0])
            HH += H2[0]
            xgauge, ygauge = H2[1:]

    # left
    # if panel==2:
    for i, distance in enumerate((3, 6, 13)):
        a = np.argmin((xgauge-distance)**2)
        # panels[i].fill(ygauge, HH[a,:]/np.max(HH[a,:]), color=cm.jet(rel), alpha=0.5)
        panels[i].plot(ygauge, HH[a,:]/np.max(HH[a,:]), color=cm.jet(rel))
        # ixmax = np.argmin(np.abs(xgauge-xmax))
        # plt.plot(ygauge, HH[ixmax,:], 'k-', lw=0.5)

        # plt.invert_xaxis()
        # plt.ylim(None,0)
        panels[i].set_xlim(ra)
        panels[i].tick_params(labelsize=14)
        if num==0:
            panels[i].annotate(f"{distance*0.1:.1f} nm", xy=(0.95, 0.8), fontsize=14,xycoords='axes fraction', horizontalalignment='right')
        panels[i].set_yticks([])
        panels[i].set_xticks([])

    # farthest data from sd_ice
    i = 3
    distance=50
    with open(f"sd_ice/{icef[num]}.d") as f:
        D = [[float(v) for v in line.split()] for line in f.readlines() if line[0] != "#"]
        D = np.array(D)
        DD = D[D[:,0]==50.0, :]
        r, avg, sd, _ = DD[0]
        Y = np.exp(-(ygauge-avg)**2/(2*sd**2))
        # panels[i].fill(ygauge, Y, color=cm.jet(rel), alpha=0.5)
        panels[i].plot(ygauge, Y, color=cm.jet(rel))
        # ixmax = np.argmin(np.abs(xgauge-xmax))
        # plt.plot(ygauge, HH[ixmax,:], 'k-', lw=0.5)

        # plt.invert_xaxis()
        # plt.ylim(None,0)
        panels[i].set_xlim(ra)
        panels[i].tick_params(labelsize=14)
        if num==0:
            panels[i].annotate(f"{distance*0.1:.1f} nm", xy=(0.95, 0.8), fontsize=14,xycoords='axes fraction', horizontalalignment='right')
        panels[i].set_yticks([])
        panels[i].set_xlabel(r"$I_i(r)$ / kJ mol$^{-1}$", fontsize=14)

fig.savefig("Figure1.pdf", bbox_inches="tight")
