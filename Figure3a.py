"""
molecular view of a model 2D ice
"""

import random
import numpy as np
from icegraph import ice_graph
import svgwrite


def draw_line(dwg, a, b, **kwargs):
    gr = dwg.g()
    d = b-a
    L = np.linalg.norm(d)
    e = d / L
    w = kwargs["stroke_width"]
    dd = e*L
    f = np.array([e[1], -e[0]])
    A = a+dd+f*w/2
    B = a+dd-f*w/2
    C = a-f*w/2
    D = a+f*w/2
    gr.add(dwg.polygon([A,B,C,D], stroke=svgwrite.rgb(0,0,0), fill=kwargs["stroke"]))
    return gr


def address(a, N=6):
    return np.array([a//N, a%N], dtype=np.float)


def draw_mol_svg(dwg, dg, N=6, zoom=50):
    gr = dwg.g()
    # oxygen
    for v in range(N*N):
        xy = address(v, N)
        gr.add(dwg.circle(center=(xy+0.5)*zoom, r=0.2*zoom,
                           fill=svgwrite.rgb(100,50,50,'%'),
                           stroke=svgwrite.rgb(0,0,0,'%')))
    for donor, acceptor in dg.edges():
        don = address(donor, N)
        acc = address(acceptor, N)
        d = acc - don
        d -= np.floor(d/N+0.5)*N
        gr.add(dwg.circle(center=(don+d/3+0.5)*zoom, r=0.1*zoom,
                        fill=svgwrite.rgb(50,100,100,'%'),
                        stroke=svgwrite.rgb(0,0,0,'%')))
        L = np.linalg.norm(d)
        s1 = (0.2**2-0.05**2)**0.5
        s2 = (0.1**2-0.05**2)**0.5
        e = d/L
        # color=random_color(c+1)
        gr.add(draw_line(dwg, (don+e*s1+0.5)*zoom, (don+e*(L/3-s2)+0.5)*zoom,
                         stroke=svgwrite.rgb(100,100,100,'%'),
                         stroke_width=zoom*0.1))
    dwg.add(gr)



random.seed(111)

# 格子点の数はN^2、 水素の数はその2倍。
N=6

# 別ファイルで128でも計算したが、ほぼ同じ傾向が得られたので、そこまで大きくする必要はなさそうだ。

# O-Oの間の水素の位置。
Hoffset = 1/3

#酸素の座標。酸素は格子点上にある。
Opos = np.zeros([N*N, 2])
X,Y = np.mgrid[0:N,0:N]
# 酸素原子の位置はずっと固定。
Opos[:, 0] = X.reshape(N*N)
Opos[:, 1] = Y.reshape(N*N)


cell = np.array([N,N])
rOpos = Opos / cell
cellmat = np.diag(cell)

g = ice_graph(N, rOpos)

zoom=50

dwg = svgwrite.Drawing(size=(zoom*N,zoom*N))
draw_mol_svg(dwg, g, zoom=zoom, N=N)
svg = dwg.tostring()
with open("Figure3a.svg", "w") as f:
    f.write(svg)
