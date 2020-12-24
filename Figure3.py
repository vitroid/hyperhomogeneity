from ice7analysis import *
from matplotlib import pyplot as plt
import pickle
import glob
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)


fig  = plt.figure(figsize=(8,8))
gs = plt.GridSpec(5,3)

ices = [["Ih", "1h", plt.subplot(gs[:1, -1:])],
        ["III", "3", plt.subplot(gs[1:2, -1:])],
        ["V", "5", plt.subplot(gs[2:3, -1:])],
        ["VI", "6", plt.subplot(gs[3:4, -1:])],
        ["VII", "7", plt.subplot(gs[:4, :2])]]


plt.subplots_adjust(wspace=0, hspace=0)
linear = np.linspace(0.0,20,1000)

water = tip4picesites()
charges = tip4picecharges()
AA,BB = tip4piceLJ()

every=5

for ice, basename, ax in ices:
    coord = f"q/{basename}-1000.q.nx3a"
    comeus, cell = load_nx3a(open(coord))
    cellmat = np.diag(cell)
    rcom = comeus[:,:3] @ np.linalg.inv(cellmat)

    Nmol = len(comeus)
    R = np.zeros([Nmol,3,3])
    for i in range(Nmol):
        e = comeus[i,3:6]
        R[i] = quat2rotmat(euler2quat(e))


    dg = hbn(rcom, cellmat, R, water)

    # 分子内座標
    waters = water @ R

    # fig = plt.figure()
    if ice == "VII":
        # 先に、分子単位での積算。

        d_e = accum0(comeus, cellmat, range(0,Nmol,every), LJ=(0.0,0.0))
        acc,sd,cnt = stepgraph(d_e, linear)
        bandplot1(linear,acc,smooth(sd), plt=ax, label="Molecule")

        #plt.legend()
        ax.set_xlim(0,13)
        bottom = int(acc[-1]/10+0.5)*10
        ax.set_ylim(bottom-20,bottom+30)


        # SDをinset
        ax2 = plt.subplot(gs[4:5, -1:])
        # Manually set the position and relative size of the inset axes within ax1
        ip = InsetPosition(ax, [0.5,0.75,0.5,0.25])
        ax2.set_axes_locator(ip)
        ax2.set_xlim(0,13)
        ax2.set_ylim(0,16)
        ax2.set_yticks([0,10])
        ax2.set_xlabel("Distance / 0.1 nm",fontsize=14)
        ax2.set_ylabel("SD",fontsize=14)
        ax2.plot(linear, smooth(sd))
        i1 = np.argmax(smooth(sd))
        r1 = linear[i1]
        v1 = smooth(sd)[i1]
        print(f"{ice} r_1 {r1}")
        ax2.plot([r1,r1],[0,v1], color='gray')


        d_e = []
        for cycles5 in glob.glob(f"{basename}-*.cycles5.pickle"):
            with open(cycles5, "rb") as f: # non-quenched is ok
                bucket = pickle.load(f)
                cycles_raw = bucket["cycles"]
                weights = bucket["weight"]
            cyclesintr = cycles5.replace("cycles5", "cyclesintr")
            print(cycles5)
            with open(cyclesintr, "rb") as f:
                intr, dist, mord = pickle.load(f)
            for center in range(Nmol):
                # cycle by cycle
                eps = []
                ds  = []
                for c, (cycle, wei) in enumerate(zip(cycles_raw, weights)):
                    ep = intr[center, c]
                    d  = dist[center, c]
                    eps.append(ep*wei)
                    ds.append(d)

                ds = np.array(ds)
                order = np.argsort(ds)
                eps = np.array(eps)
                # plt.plot(ds[order],np.cumsum(eps[order]), label=label)
                # print(np.cumsum(eps[order])[-1])
                d_e.append([ds[order],np.cumsum(eps[order])])
        if len(d_e) > 0:
            acc,sd,cnt = stepgraph(d_e, linear)
            bandplot1(linear,acc,smooth(sd), plt=ax, label=f"Cycle ({ice})")

        # Create a set of inset Axes: these should fill the bounding box allocated to
        # them.
        ax2.annotate(f"ice {ice}", (0.6, 0.6), xycoords="axes fraction",fontsize=18)
        ax2.plot(linear, smooth(sd))


        #plt.ylim(-150,-130)
        #plt.legend()
        ax.set_xlabel("Distance / 0.1 nm",fontsize=18)
        ax.set_ylabel(r"Coulomb interaction / kJ mol$^{-1}$",fontsize=18)
    else:
        # 先に、分子単位での積算。

        d_e = accum0(comeus, cellmat, range(0,Nmol,every), LJ=(0.0,0.0))
        acc,sd,cnt = stepgraph(d_e, linear)

        ax.set_xlim(0.0,13)
        ax.set_xlabel("Distance / 0.1 nm",fontsize=14)
        if ice == "Ih":
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
        else:
            # no tick no label
            ax.xaxis.set_visible(False)

        ax.set_ylim(0,16)
        ax.yaxis.tick_right()
        ax.set_yticks([0,10])
        if ice == "III":
            ax.yaxis.set_label_position('right')
            ax.set_ylabel(r"SD / kJ mol$^{-1}$",fontsize=18)

        ax.plot(linear, smooth(sd))
        i1 = np.argmax(smooth(sd[linear>3.5]))
        r1 = linear[linear>3.5][i1]
        v1 = smooth(sd[linear>3.5])[i1]
        print(f"{ice} r_1 {r1}")
        ax.plot([r1,r1],[0,v1], color='gray')


#         assert False, "遅すぎるので,cycleintr.pickleを使え。"

        d_e = []
        for cycles5 in glob.glob(f"{basename}-*.cycles5.pickle"):
            with open(cycles5, "rb") as f: # non-quenched is ok
                bucket = pickle.load(f)
                cycles_raw = bucket["cycles"]
                weights = bucket["weight"]
            cyclesintr = cycles5.replace("cycles5", "cyclesintr")
            with open(cyclesintr, "rb") as f:
                intr, dist, mord = pickle.load(f)
            print(cycles5)
            for center in range(Nmol):
                # cycle by cycle
                eps = []
                ds  = []
                for c, (cycle, wei) in enumerate(zip(cycles_raw, weights)):
                    ep = intr[center, c]
                    d  = dist[center, c]
                    eps.append(ep*wei)
                    ds.append(d)

                ds = np.array(ds)
                order = np.argsort(ds)
                eps = np.array(eps)
                # plt.plot(ds[order],np.cumsum(eps[order]), label=label)
                # print(np.cumsum(eps[order])[-1])
                d_e.append([ds[order],np.cumsum(eps[order])])
        if len(d_e) > 0:
            acc,sd,cnt = stepgraph(d_e, linear)
            ax.annotate(f"ice {ice}", (0.6, 0.6), xycoords="axes fraction",fontsize=18)
            ax.plot(linear, smooth(sd))
#             bandplot1(linear,acc,smooth(sd), plt=ax, label=f"Cycle ({ice})")

        # Create a set of inset Axes: these should fill the bounding box allocated to
        # them.

# plt.show()
fig.savefig("Figure3.pdf", bbox_inches="tight")
