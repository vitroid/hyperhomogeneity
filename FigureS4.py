from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.cm as cm



from ice7analysis import *

basename = "1cs"

coord = f"{basename}.nx3a"
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
print(layer)

# fig, ax = plt.subplots(1,3)
fig  = plt.figure(figsize=(5,10))
grid = plt.GridSpec(5, 1, height_ratios=[7,7,7,7,7], wspace=0, hspace=0)
linear = np.linspace(2.5,13.5,1000)

import pickle

layerlabels = ["top layer", "2nd", "3rd", "4th", "5th"]

for lay, layerlabel in enumerate(layerlabels):
    panel = lay
    if panel == 0:
        ra = (-100,-30)
    else:
        ra = (-150,-80)
    bot = (lay <= layer) & (layer < lay+1)
    top = ((31-lay) <= layer) & (layer < (32-lay))
    surface = top | bot
    print(lay, np.nonzero(surface))


    d_e = accum0(comeus, cellmat, np.nonzero(surface)[0], maxdist=13.2)
    N = len(d_e)


    main_ax = fig.add_subplot(grid[panel,0], xticks=[3,8,13], yticks=(-140,-100,-60,-20))

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

    # center panel
    main_ax.set_xlabel("Distance / 0.1 nm")
    main_ax.set_xlim(3,13)
    main_ax.label_outer()
    main_ax.set_ylim(ra)
    main_ax.annotate(f"{layerlabel}", xy=(0.95,0.8), fontsize=18,xycoords='axes fraction', horizontalalignment='right')
#     main_ax.annotate("c = {0:.1f}".format(convergence), xy=(0.95,0.6), fontsize=14,xycoords='axes fraction', horizontalalignment='right')
#     main_ax.plot([xmax, xmax], [0, -1000], "g-", lw=1.0)
    main_ax.tick_params(labelsize=14)

    acc,sd,cnt = stepgraph(d_e, linear)

#     # SDをinset
#     ax2 = plt.subplot(grid[panel,0:1])
#     # Manually set the position and relative size of the inset axes within ax1
#     ip = InsetPosition(main_ax, [0.5,0.6,0.5,0.4])
#     ax2.set_axes_locator(ip)
#     ax2.set_xlim(0,13)
#     ax2.set_ylim(0,10)
#     ax2.set_xlabel("Distance / 0.1 nm")
#     ax2.set_ylabel("SD")
#     print(1)
#     ax2.plot(linear, smooth(sd))
#     print(2)

#     HH, xgauge, ygauge = H2[ice]

#     # left
    if panel==2:
        main_ax.set_ylabel(r"$I_i(r)$"+r" / kJ mol$-1$")
#     hist3.fill(HH[0,:], ygauge, color=cm.viridis(panel/4))
#     ixmax = np.argmin(np.abs(xgauge-xmax))
#     hist3.plot(HH[ixmax,:], ygauge, 'k-', lw=0.5)

#     hist3.invert_xaxis()
#     hist3.set_xlim(None,0)
#     hist3.set_ylim(ra)
#     hist3.tick_params(labelsize=14)

#     # right
#     hist13.fill(HH[-1,:], ygauge, color=cm.viridis(panel/4))
#     hist13.set_xlim(0,None)
#     hist13.label_outer()
#     hist13.set_ylim(ra)



#for a in ax.flat:
#    # Hide x labels and tick labels for top plots and y ticks for right plots.
#    a.label_outer()
#    a.set_ylim(-150,-80)
fig.savefig("FigureS4.pdf", bbox_inches="tight")
