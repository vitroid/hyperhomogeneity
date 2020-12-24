from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm
from ice7analysis import *
import glob

# fig, ax = plt.subplots(1,3)
fig  = plt.figure(figsize=(5,10))
grid = plt.GridSpec(5, 3, width_ratios=[1, 3, 1], height_ratios=[7,7,7,7,11], wspace=0, hspace=0)
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


for panel, (ice, ra) in enumerate([["1h", [-150,-80]],
                                         ["3",  [-150,-80]],
                                         ["5",  [-150,-80]],
                                         ["6",  [-150,-80]],
                                         ["7",  [-130,-20]]]):
    comeus, cell = load_nx3a(open(f"q/{ice}-1000.q.nx3a"))
    cellmat = np.diag(cell)
    d_e = accum0(comeus, cellmat, range(50), maxdist=13.2)
    N = len(d_e)

    main_ax = fig.add_subplot(grid[panel,1], xticks=[3,8,13])
    hist3   = fig.add_subplot(grid[panel,0], xticklabels=[], sharey=main_ax, yticks=(-140,-100,-60))
    hist13  = fig.add_subplot(grid[panel,2], xticklabels=[], sharey=main_ax)

    # 距離4でのエネルギーの値でソートする。

    e4 = np.zeros(N)
    for i, (d,e) in enumerate(d_e):
        e4[i] = e[d<4.0][-1]
    order = np.argsort(e4)
    redro = np.zeros(N)
    for i in range(N):
        redro[order[i]] = i
    for i in range(N):
        j = redro[i]
        avg,sd,cnt = stepgraph(d_e[i:i+1], linear) # single data
        main_ax.plot(linear, avg, color=cm.coolwarm(j/N), lw=0.5)

    convergence, xmax, ymax, ylast = conv[ice]
    # center panel
    main_ax.set_xlabel("Distance / 0.1 nm")
    main_ax.set_xlim(3,13)
    main_ax.label_outer()
    main_ax.set_ylim(ra)
    main_ax.annotate("ice {0}".format(ice), xy=(0.95,0.8), fontsize=18,xycoords='axes fraction', horizontalalignment='right')
    main_ax.annotate("c = {0:.1f}".format(convergence), xy=(0.95,0.6), fontsize=14,xycoords='axes fraction', horizontalalignment='right')
    main_ax.plot([xmax, xmax], [0, -1000], "g-", lw=1.0)
    main_ax.tick_params(labelsize=14)

    if ice == "1h":
        coord11, cell11 = load_nx3a(open("11.nx3a"))
        cellmat11 = np.diag(cell11)
        rpos11 = coord11[:, :3]
        d_e11 = accum0(coord11, cellmat11, range(2), maxdist=13.2)
        N11 = len(d_e11)
        for i in range(N11):
            avg,sd,cnt = stepgraph(d_e11[i:i+1], linear) # single data
            main_ax.plot(linear, avg, "k", lw=1.0)



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
    if panel==2:
        hist3.set_ylabel(r"$I_i(r)$ / kJ mol$-1$")
    hist3.fill(HH[0,:], ygauge, color=cm.viridis(panel/4))
    ixmax = np.argmin(np.abs(xgauge-xmax))
    hist3.plot(HH[ixmax,:], ygauge, 'k-', lw=0.5)

    hist3.invert_xaxis()
    hist3.set_xlim(None,0)
    hist3.set_ylim(ra)
    hist3.tick_params(labelsize=14)

    # right
    hist13.fill(HH[-1,:], ygauge, color=cm.viridis(panel/4))
    hist13.set_xlim(0,None)
    hist13.label_outer()
    hist13.set_ylim(ra)



#for a in ax.flat:
#    # Hide x labels and tick labels for top plots and y ticks for right plots.
#    a.label_outer()
#    a.set_ylim(-150,-80)
fig.savefig("Figure1.pdf", bbox_inches="tight")
