"""
cycle view of a model 2D ice
"""


import svgwrite
from collections import defaultdict
import numpy as np
import random
import networkx as nx
from cycles5 import cover_by_cycles
from icegraph import ice_graph


def address(a, N=6):
    return np.array([a//N, a%N], dtype=np.float)




def arrow(dwg, a, b, **kwargs):
    gr = dwg.g()
    d = b-a
    L = np.linalg.norm(d)
    e = d / L
    w = kwargs["stroke_width"]
    h = w*3
    if h > L:
        h = L
    dd = e*(L-h)
    f = np.array([e[1], -e[0]])
    A = b
    B = a+dd+f*w
    C = a+dd-f*w
    gr.add(dwg.polygon([A,B,C], stroke_width=0, fill=kwargs["stroke"]))
    A = a+dd+f*w/2
    B = a+dd-f*w/2
    C = a-f*w/2
    D = a+f*w/2
    gr.add(dwg.polygon([A,B,C,D], stroke_width=0, fill=kwargs["stroke"]))
    return gr


def draw_arrows_svg(dwg, graph, zoom=50, N=6):
    # cycles
    # 事前に、外周をまたぐ部分で分割しておいたほうがいい。
    gr = dwg.g()
    for edge in graph.edges():
        a,b = edge
        A = address(a)
        B = address(b)
        d = B - A
        wr = np.floor(d/N+0.5)*N
        gr.add(arrow(dwg, (A+0.5)*zoom, (B - wr +0.5)*zoom,
                        stroke=svgwrite.rgb(100,100,100,'%'),
                        stroke_width=zoom*0.1))
        if not np.allclose(wr, np.zeros_like(wr)):
            gr.add(arrow(dwg, (A+wr+0.5)*zoom, (B +0.5)*zoom,
                            stroke=svgwrite.rgb(100,100,100,'%'),
                            stroke_width=zoom*0.1))
    dwg.add(gr)

    gr = dwg.g()
    for v in graph:
        A = address(v)
        gr.add(dwg.circle(center=(A+0.5)*zoom, r=0.1*zoom,
                           fill=svgwrite.rgb(0,0,0,'%'),
                           ))
    dwg.add(gr)


def draw_path(dwg, path, cyclic=False, **kwargs):
    magic = 0.552284749831 # for 90 degrees only
    p = [list(x) for x in path]
    if cyclic:
        return dwg.polygon(points=p, **kwargs)
    return dwg.polyline(points=p, **kwargs)


def rearrange_path(p):
    """
    周期境界をまたぐ部分で分割する。
    """
    paths = []
    chain = []
    for i in range(0,p.shape[0]):
        a = p[i-1]
        b = p[i]
        d = b-a
        wr = np.floor(d/N+0.5)*N
        if not np.allclose(wr, np.zeros_like(wr)):
            # spans the boundary
            chain.append(b-wr)
            paths.append(chain)
            chain = [a+wr]
        chain.append(b)
    paths.append(chain)
    for path in paths:
        print(path)
    print("==")
    if np.allclose(paths[-1][-1],paths[0][0]) and len(paths) > 1:
        paths[0] = paths[-1][:-1] + paths[0]
        paths.pop(-1)
    for path in paths:
        print(path)
    print("==")
    return paths


def random_color(c):
    import colorsys
    tau = (1+5**0.5)/2.0
    return colorsys.hsv_to_rgb((c/tau)%1.0, 1.0, 1.0)


def draw_cycles_svg(dwg, cycles, zoom=50):
    count = defaultdict(int)
    for cycle in cycles:
        count[len(cycle)] += 1
    print(count)

    # cycles
    # 事前に、外周をまたぐ部分で分割しておいたほうがいい。
    for c, path in enumerate(cycles):
        gr = dwg.g()
        p = []
        for node in path:
            xy = address(node)
            p.append((xy))
        p.append(p[0]) #close path
        p = np.array(p, dtype=np.float)
        color=np.array(random_color(c+1))
        paths = rearrange_path(p)
        for chain in paths:
            chain = np.array(chain)
            gr.add(draw_path(dwg, (chain+0.5)*zoom, cyclic=np.allclose(chain[0], chain[-1]),
                           stroke=svgwrite.rgb(*(color*100),'%'),
                           stroke_width=zoom*0.4,
                           stroke_linejoin="round",
                          fill="none"
                           ))
            for i in range(len(chain)-1):
                v1 = chain[i]
                v2 = chain[i+1]
                gr.add(arrow(dwg, (v1+0.5)*zoom, (v2+0.5)*zoom,
                                stroke=svgwrite.rgb(100,100,100,'%'),
                                stroke_width=zoom*0.1))
            for v in chain:
                gr.add(dwg.circle(center=(v+0.5)*zoom, r=0.1*zoom,
                           fill=svgwrite.rgb(0,0,0,'%'),
                           ))
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
draw_arrows_svg(dwg, g)
cycles, weights = cover_by_cycles(g, N*N)
draw_cycles_svg(dwg, cycles)
svg = dwg.tostring()
with open("Figure3bc.svg", "w") as f:
    f.write(svg)
