"""
cycles5

あらかじめ小さいすべての環を列挙して、それで特異値分解して有向グラフを環の線形和にする。
"""
from ice7analysis import *
import pickle
import sys
import random
from cycless.dicycles import dicycles_iter



def AtoX(A):
    # 環の個数のほうが十分多い場合の近似
    # 特異値分解
    Q1, S, Q2T = np.linalg.svd(A)
    # print(Q1.shape, S.shape, Q2T.shape)
    # ほぼ0の要素を0にする
    S[np.abs(S)<1e-12] = 0
    rank = np.sum(np.abs(S) > 0)
    # print(S,r)
    # SS = Q1.T @ A @ Q2T.T
    # 対角行列S†の準備
    Sd = np.zeros_like(A).T
    Sd[:rank, :rank] = np.diag(1/S[:rank])
    # print(SS)
    # print(Sd.shape)
    # A†
    Ad = Q2T.T @ Sd @ Q1.T
    #print(Ad.shape)
    b = np.ones(A.shape[0])
    x = Ad@b
    # print(A@x)
    # print(x)
    return x, rank


def cover_by_cycles(g, Nmol):
    HBcode = {(d,a): x for x, (d,a) in enumerate(g.edges())}

    A = []
    cycles = []
    for s in range(4, 30):
        print(f"cycle size {s}")
        lastA = len(A)
        for cycle in dicycles_iter(g,s, vec=True):
            cycles.append(cycle)
            row = np.zeros(Nmol*2)
            for j in range(len(cycle)):
                edge = HBcode[cycle[j-1], cycle[j]]
                row[edge] = 1.0
            A.append(row)
        if lastA==len(A):
            continue
        lastA = len(A)
        AT = np.array(A).T
        x, rank = AtoX(AT)
        print(f"number of cycles {len(cycles)}")
        b = AT @ x
        print(b)
        if np.allclose(b, np.ones_like(b)):
            break
    return cycles, x

if __name__ == "__main__":
    coord = sys.argv[1]
    comeus, cell = load_nx3a(open(coord))
    cellmat = np.diag(cell)
    rcom = comeus[:,:3] @ np.linalg.inv(cellmat)
    Nmol = len(comeus)
    R = np.zeros([Nmol,3,3])
    for i in range(Nmol):
        e = comeus[i,3:6]
        R[i] = quat2rotmat(euler2quat(e))
    water = tip4picesites()
    g = hbn(rcom, cellmat, R, water)

    # 氷7でうまく分割できない問題を解決するために、成分に分ける。
    # 氷がうまくいかないのは、脱分極がちゃんとできていないから。

    # for compo in nx.connected_components(nx.Graph(g)):
    #     subgraph = g.subgraph(compo) # compo is a set of node labels.
    #
    #     for node in subgraph:
    #         assert subgraph.in_degree(node) == 2
    #         assert subgraph.out_degree(node) == 2
    #     A = GenerateA(subgraph)
    #     A = np.array(A).T
    #     x, rank = AtoX(A)
    #     print(x, rank, A.shape)
    #     b = A @ x
    #     print(b)
    cycles, weights =cover_by_cycles(g, Nmol)

    if len(sys.argv) > 2:
        with open(sys.argv[2], "wb") as f:
            pickle.dump({"cycles":cycles, "weight":weights}, f)
