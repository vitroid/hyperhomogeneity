"""
Common functions.

It contains many functions that are not used for historical reasons.
"""

import numpy as np
from math import sin,cos,pi,sqrt, acos
from matplotlib import pyplot as plt
from collections import defaultdict
import yaplotlib as yp
import vapory as pov
import networkx as nx
import pairlist as pl
import itertools as it

#
# SI unit 2019
Na=6.02214076e23    #Avogadro
ee=1.60217662e-19   #Elementary charge
E0=8.8541878128e-12 #Vacuum permittivity
kB=1.38064852e-23
qe0 = ee**2/(4*pi*E0)*Na*1e7 # in kJ/mol
UJ = 4.184          #1 cal in J


def load_cell(file):
    """
    Special function to read the cell size from Tanaka's format.
    """
    return np.array([float(x) for x in file.readline().split()[4:7]])


def load_nx3a(file):
    """
    Centers of mass and Euler angles.
    """
    coms = []
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if "@BOX3" in line:
            line = file.readline()
            cell = np.array([float(x) for x in line.split()])
        elif "@NX3A" in line:
            line = file.readline()
            nmol = int(line)
            for i in range(nmol):
                line = file.readline()
                comeu = [float(x) for x in line.split()[:6]]
                coms.append(comeu)
    return np.array(coms), cell


def load_coms(file):
    """
    Centers of mass and Euler angles; Tanaka's format.
    """
    coms = []
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if line[0] == "#":
            continue
        cols = line.split()
        if len(cols) == 7:
            comeu = [float(x) for x in line.split()[1:7]]
            coms.append(comeu)
        else:
            print(line, end="")
    return np.array(coms)

def load_coms2(file):
    """
    Another format of Tanaka's
    """
    coms = []
    cell = None
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if line[0] == "#":
            cols = line.split()
            cell = np.array([float(x) for x in cols[2:5]])
            continue
        cols = line.split()
        if len(coms) == 0 and len(cols) == 2:
            # first special line without # is for cubic cell
            cell = float(cols[1])
            cell = np.array([cell,cell,cell])
            continue
        if len(cols) == 7:
            comeu = [float(x) for x in line.split()[1:7]]
            coms.append(comeu)
        else:
            print(line, end="")
    return np.array(coms), cell





def quat2rotmat(q):
    """
    Quaternionを回転行列(後置形式)にする。
    """
    a, b, c, d = q
    sp11 = (a * a + b * b - (c * c + d * d))
    sp12 = -2.0 * (a * d + b * c)
    sp13 = 2.0 * (b * d - a * c)
    sp21 = 2.0 * (a * d - b * c)
    sp22 = a * a + c * c - (b * b + d * d)
    sp23 = -2.0 * (a * b + c * d)
    sp31 = 2.0 * (a * c + b * d)
    sp32 = 2.0 * (a * b - c * d)
    sp33 = a * a + d * d - (b * b + c * c)
    return np.array([[sp11, sp12, sp13], [sp21, sp22, sp23], [sp31, sp32, sp33]]).T




def euler2quat(e):
    """
    Euler角をQuaternionにする。
    """
    ea, eb, ec = e
    a = cos(ea / 2) * cos((ec + eb) / 2)
    b = sin(ea / 2) * cos((ec - eb) / 2)
    c = sin(ea / 2) * sin((ec - eb) / 2)
    d = cos(ea / 2) * sin((ec + eb) / 2)
    return np.array((a, b, c, d))


def steppify(xy, linear):
    """
    xy is a list of x-to-y data sequences
    """
    for x,y in xy:
        N = x.shape[0]
        stepx=np.zeros(N*2-1)
        stepy=np.zeros(N*2-1)
        stepx[0::2] = x
        stepx[1:2*N-1:2] = x[1:] - 1e-8
        stepy[0::2] = y
        stepy[1::2] = y[:-1]
        yield np.interp(linear, stepx, stepy)


def stepgraph(xy, linear):
    """
    xy is a list of x-to-y data sequences

    return average line and SDs
    """
    acc = np.zeros_like(linear)
    sd = np.zeros_like(linear)

    for sf in steppify(xy, linear):
        acc += sf
        sd += sf**2

    cnt = len(xy)
    acc /= cnt
    sd /= cnt
    sd -= acc**2
    sd = sd**0.5
    return acc, sd, cnt

def bandplot1(linear, avg, sd, plt=plt, **kwarg):
    """
    xy is a list of x-to-y data sequences
    """
    plt.fill_between(linear, avg-sd, avg+sd, alpha=0.5)
    plt.plot(linear, avg, **kwarg)


def bandplot0(xy, range=(2.5, 13), **kwarg):
    """
    xy is a list of x-to-y data sequences
    """
    linear = np.linspace(range[0], range[1], 1000)
    acc, sd, cnt = stepgraph(xy, linear)
    bandplot1(linear, acc,sd, **kwarg)
    plt.axhline()


def tip5psites():
    """
    TIP5P Geometry
    """
    L1 = 0.9572
    L2 = 0.70
    theta=104.52 * pi/180
    phi  = 109.47*pi/180


    hy = L1*sin(theta/2)
    hz = L1*cos(theta/2)
    mx = L2*sin(phi/2)
    mz = L2*cos(phi/2)

    # HHOMM
    sites = np.array([[0.0, hy,  hz],
                      [0.0,-hy,  hz],
                      [0.0, 0.0, 0.0],
                      [mx, 0.0, -mz],
                      [-mx,0.0, -mz]])
    sites -= (sites[0]+sites[1]+sites[3]*0)/18
    return sites

def tip5pcharges():
    """
    TIP5P Charge
    """
    q = 0.2410
    return np.array([q,q,0,-q,-q])

def tip5pLJ():
    """
    TIP5P integrated LJ parameters
    from Chaplin
    http://www1.lsbu.ac.uk/water/water_models.html
    """
    eps = 0.6694e3
    sig = 3.12
    AA = 4*eps*sig**12 / 1000
    BB = 4*eps*sig**6 / 1000
    return AA,BB



def tip4psites():
    """
    TIP4P Geometry
    """
    L1 = 0.9572
    L2 = 0.15
    theta=104.52 * pi/180


    hy = L1*sin(theta/2)
    hz = L1*cos(theta/2)
    mz = L2

    # HHOM
    sites = np.array([[0.0, hy,  hz],
                      [0.0,-hy,  hz],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, mz]])
    sites -= (sites[0]+sites[1]+sites[3]*0)/18
    return sites

def tip4pcharges():
    """
    TIP4P Charge
    """
    q = 0.52
    return np.array([q,q,0,-2*q])

# TIP4P/2005 water
def tip4picesites():
    """
    TIP4P/Ice Geometry
    """
    L1 = 0.9572
    L2 = 0.1577
    theta=104.52 * pi/180


    hy = L1*sin(theta/2)
    hz = L1*cos(theta/2)
    mz = L2

    # HHOM
    sites = np.array([[0.0, hy,  hz],
                      [0.0,-hy,  hz],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, mz]])
    sites -= (sites[0]+sites[1]+sites[3]*0)/18
    return sites


# TIP4P/2005 water
def tip4p2005sites():
    """
    TIP4P/2005 Geometry
    """
    L1 = 0.9572
    L2 = 0.1546
    theta=104.52 * pi/180


    hy = L1*sin(theta/2)
    hz = L1*cos(theta/2)
    mz = L2

    # HHOM
    sites = np.array([[0.0, hy,  hz],
                      [0.0,-hy,  hz],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, mz]])
    sites -= (sites[0]+sites[1]+sites[3]*0)/18
    return sites

def tip4picecharges():
    """
    TIP4P/Ice Charge
    """
    q = 0.5897
    return np.array([q,q,0,-2*q])

def tip4p2005charges():
    """
    TIP4P/2005 Charge
    """
    q = 0.5564
    return np.array([q,q,0,-2*q])


def tip4pLJ():
    """
    TIP4P integrated LJ parameters
    from SKlog
    http://www.sklogwiki.org/SklogWiki/index.php/TIP4P_model_of_water
    """
    eps = 78 * kB * Na
    sig = 3.154
    AA = 4*eps*sig**12 / 1000
    BB = 4*eps*sig**6 / 1000
    print(AA,BB)
    AA = 6e5 * UJ  # kJ/mol/AA**12
    BB = 610 * UJ  # kJ/mol/AA**6
    return AA,BB


def tip4piceLJ():
    """
    TIP4P/Ice integrated LJ parameters
    from SKlog
    http://www.sklogwiki.org/SklogWiki/index.php/TIP4P/Ice_model_of_water
    """
    eps = 106.1 * kB * Na
    sig = 3.1668
    AA = 4*eps*sig**12 / 1000
    BB = 4*eps*sig**6 / 1000
    return AA,BB


def tip4p2005LJ():
    """
    TIP4P/2005 integrated LJ parameters
    from SKlog
    http://www.sklogwiki.org/SklogWiki/index.php/TIP4P/Ice_model_of_water
    """
    eps = 93.2 * kB * Na
    sig = 3.1589
    AA = 4*eps*sig**12 / 1000
    BB = 4*eps*sig**6 / 1000
    return AA,BB



def interactions4px(target, atoms, wrap, charges, AA, BB):
    """
    target: atoms of the central molecule
    atoms: atomic positions of all molecules
    wrap: PBC wrap vector between target and others
    """
    Cs = np.zeros(atoms.shape[0])
    LJs = np.zeros(atoms.shape[0])

    for j in range(atoms.shape[0]):
        # coulobs
        coulomb = 0.0
        for a in (0,1,3):
            for b in (0,1,3):
                dd = atoms[j,a] - target[b] - wrap[j]
                L1 = np.linalg.norm(dd)
                coulomb += charges[a]*charges[b] / L1
        # to kJ/mol
        coulomb *= qe0
        dd = atoms[j,2] - target[2] - wrap[j]
        Q = dd @ dd
        LJ = AA/Q**6 - BB/Q**3
        Cs[j] = coulomb
        LJs[j] = LJ
    return Cs, LJs

def interactions4p(target, atoms, wrap, charges, AA, BB):
    """
    target: central molecule
    atoms: atomic positions of all molecules
    wrap: PBC wrap vector between target and others
    """
    return interactions4px(atoms[target], atoms, wrap, charges, AA, BB)




# 基本にたちかえって、その場ですべて計算する。
def interaction4ps(i, j, comeus, cellmat, water, charges, AA, BB):
    comi = comeus[i,:3]
    euli = comeus[i,3:6]
    Ri   = quat2rotmat(euler2quat(euli))
    atomi = water @ Ri
    comj = comeus[j,:3]
    eulj = comeus[j,3:6]
    Rj   = quat2rotmat(euler2quat(eulj))
    atomj = water @ Rj

    d = comj - comi
    d -= np.floor( d @ np.linalg.inv(cellmat) + 0.5 ) @ cellmat
    ev = d / np.linalg.norm(d)
    cf = 0.0
    for a in (0,1,3):
        for b in (0,1,3):
            # 電荷間相対位置ベクトル
            dd = atomj[b] - atomi[a] + d
            L1 = np.linalg.norm(dd)
            f = charges[a]*charges[b] / L1
            cf += f
    cf *= qe0
    dd = atomj[2] - atomi[2] + d
    Q = dd @ dd
    f = AA/Q**6 - BB/Q**3
    ljf = f
    # force strength along com-com direction (two scalars)
    return cf, ljf



def interaction4pt(Ri, Rj, dij, water, charges, AA, BB):
    atomi = water @ Ri
    atomj = water @ Rj

    ev = dij / np.linalg.norm(dij)
    cf = 0.0
    for a in (0,1,3):
        for b in (0,1,3):
            # 電荷間相対位置ベクトル
            dd = atomj[b] - atomi[a] + dij
            L1 = np.linalg.norm(dd)
            f = charges[a]*charges[b] / L1
            cf += f
    cf *= qe0
    dd = atomj[2] - atomi[2] + dij
    Q = dd @ dd
    f = AA/Q**6 - BB/Q**3
    ljf = f
    # force strength along com-com direction (two scalars)
    return cf, ljf





# 基本にたちかえって、その場ですべて計算する。
def force4p(i, j, comeus, cellmat, water, charges, AA, BB):
    comi = comeus[i,:3]
    euli = comeus[i,3:6]
    Ri   = quat2rotmat(euler2quat(euli))
    atomi = water @ Ri
    comj = comeus[j,:3]
    eulj = comeus[j,3:6]
    Rj   = quat2rotmat(euler2quat(eulj))
    atomj = water @ Rj

    d = comj - comi
    d -= np.floor( d @ np.linalg.inv(cellmat) + 0.5 ) @ cellmat
    ev = d / np.linalg.norm(d)
    cf = 0.0
    for a in (0,1,3):
        for b in (0,1,3):
            # 電荷間相対位置ベクトル
            dd = atomj[b] - atomi[a] + d
            L1 = np.linalg.norm(dd)
            f = charges[a]*charges[b] / L1**3 * dd # vector
            cf += f @ ev
    cf *= qe0
    dd = atomj[2] - atomi[2] + d
    Q = dd @ dd
    f = (12*AA/Q**7 - 6*BB/Q**4) * dd # vector
    ljf = f @ ev
    # force strength along com-com direction (two scalars)
    return cf, ljf



def double_system(cell, pos):
    Nmol = pos.shape[0]
    newpos = np.zeros([Nmol*8, 6])
    j = 0
    for x in (0.0, 1.0):
        for y in (0.0, 1.0):
            for z in (0.0, 1.0):
                D = np.array([x,y,z]) @ cell
                newpos[j:j+Nmol, :3] = D + pos[:,:3]
                newpos[j:j+Nmol, 3:6] = pos[:,3:6]
                j += Nmol
    return cell*2, newpos



# RI[i]をかけることで正しい格子点に変換されているかどうか。
# これにより、けっこう大きな変位があることがわかる。
# old name: test_rpos
def rpos2loc(pos, ep=None, d0=0.41, d1=2.0, dz=1.2):
    assert ep is None
    msg  = "{0} {1}".format(pos,ep)
    if (pos[2] > dz and abs(pos[1]) > d0) or (pos[2] < -dz and abs(pos[1]) < d0):
#        assert ep <  -10
        loc = 0
        if pos[2] < -dz:
            loc += 0
            assert not -d1 < pos[0] < d1, msg
            assert -d0 < pos[1] < d0, msg
            if pos[0] > d1:
                loc += 1
        elif pos[2] > +dz:
            loc += 2
            assert not -d1 < pos[1] < d1, msg
            assert -d0 < pos[0] < d0, msg
            if pos[1] > d1:
                loc += 1
        else:
            assert False, msg
    else:
        loc = 4
        if pos[2] > dz:
            loc += 2
            assert not -d1 < pos[0] < d1, msg
            assert -d0 < pos[1] < d0, msg
            if pos[0] > d1:
                loc += 1
        elif pos[2] < -dz:
            loc += 0
            assert not -d1 < pos[1] < d1, msg
            assert -d0 < pos[0] < d0, msg
            if pos[1] > d1:
                loc += 1
        else:
            assert False, msg
    #if loc == 6:
    #    print(pos @ R45)
    #    assert False
    return loc
#loc is a 3-bit number
# top bit==0 : on the main network.
#          4 : on the other network.
# 2nd bit==0 : bottom side (oxygen side)
#          2 : top side (hydrogen side)
# bottom bit=0: minor side of the axis
#            1: major side of the axis


def orientation(Ri, Rj):
    # Ri is the rotmat of i, not inverted.
    # iの重心から見たjの配向
    z = np.array([0.0, 0.0, 1.0])
    zj = z @ Rj
    zz = zj @ np.linalg.inv(Ri)
    if zz[2] < -0.9:
        ori = 4
    elif zz[2] > 0.9:
        ori = 5
    elif zz[0] > 0.5:
        if zz[1] > 0.5:
            ori = 0
        elif zz[1] < -0.5:
            ori = 3
        else:
            assert False
    elif zz[0] < -0.5:
        if zz[1] > 0.5:
            ori = 1
        elif zz[1] < -0.5:
            ori = 2
        else:
            assert False
    else:
        assert False

    # print(zz) # dipole direction
    return ori


def location_orientation(rcomi, Ri, rcomj, Rj, cellmat):
    # Ri is the rotmat of i, not inverted.
    #まず、相対位置をおまかに判定する。
    # iの重心から見たjの重心位置。
    d = rcomj - rcomi
    d -= np.floor(d+0.5)
    D = d @ cellmat
    # i の座標系で見た、ljの位置
    jpos = D @ np.linalg.inv(Ri)
    # 位置コード
    loc = rpos2loc(jpos)
    ori = orientation(Ri,Rj)
    return loc,ori

_tau = 2**0.5/2
R45 = np.array([[ _tau,_tau,0],
                [-_tau,_tau,0],
                [0,0,1]])



def draw_water_relative_to(comeu, origin, cell):
    d = comeu[:3] - origin[:3]
    celli = np.linalg.inv(cell)
    d -= np.floor( d @ celli + 0.5) @ cell

    R =  quat2rotmat(euler2quat(comeu[3:6]))
    Ri = np.linalg.inv(quat2rotmat(euler2quat(origin[3:6])))

    water = (tip4picesites() @ R + d) @ Ri
    s = ""
    s += yp.Color(2) # white
    s += draw_water(water)
    d = d @ Ri
    s += yp.Color(0) # black
    o = np.zeros(3)
    z = np.array([0,0,d[2]])
    yz = np.array([0,d[1],d[2]])
    xyz = d
    s += yp.Line(o,z)
    s += yp.Line(z,yz)
    s += yp.Line(yz,xyz)

    return s


def draw_water(water, povray=False):
    s = ""
    if povray:
        return [pov.Sphere(water[0], "RH", pov.Material("MATH")),
            pov.Sphere(water[1], "RH", pov.Material("MATH")),
            pov.Sphere(water[2], "RO", pov.Material("MATO")),
            pov.Cylinder(water[0], water[2], "ROH", pov.Material("MATOH")),
            pov.Cylinder(water[1], water[2], "ROH", pov.Material("MATOH"))]
    else:
        s += yp.Color(5)
        s += yp.Size(0.2)
        s += yp.Circle(water[0])
        s += yp.Circle(water[1])
        s += yp.Color(4)
        s += yp.Size(0.4)
        s += yp.Circle(water[2])
        s += yp.Color(2)
        d0 = water[0] - water[2]
        L0 = np.linalg.norm(d0)
        s += yp.Line(water[2] + d0/L0*0.4, water[0] - d0/L0*0.2)
        d1 = water[1] - water[2]
        L1 = np.linalg.norm(d1)
        s += yp.Line(water[2] + d1/L1*0.4, water[1] - d1/L1*0.2)
        # draw a tetrahedron
        y = water[1] - water[0]
        z = (water[1] + water[0]) / 2 - water[2]
        x = np.cross(y,z)
        x /= np.linalg.norm(x)
        y /= np.linalg.norm(y)
        z /= np.linalg.norm(z)
        com = (water[1] + water[0] + water[2]*16)/18
        R = 2.76/2
        a = y*(2/3)**0.5 + z*(1/3)**0.5
        b = -y*(2/3)**0.5 + z*(1/3)**0.5
        c = x*(2/3)**0.5 - z*(1/3)**0.5
        d = -x*(2/3)**0.5 - z*(1/3)**0.5
        for e,f in it.combinations((a,b,c,d), 2):
            s += yp.Line(com+e*R,com+f*R)
    return s



# outer product of two vectors; return None if vector is too small
def op(i, j, check=True):
    if check and (np.linalg.norm(i) < 0.001 or np.linalg.norm(j) < 0.001):
        return None
    a = np.cross(i, j)
    if check and np.linalg.norm(a) < 0.001:
        return None
    return a


EX = np.array((1.0, 0.0, 0.0))
EY = np.array((0.0, 1.0, 0.0))
EZ = np.array((0.0, 0.0, 1.0))

# calculate quaternions from a rotation matrix (three orthogonal unit vectors)
def rotmat2quat0(i, j, k):
    a = op(i - EX, j - EY)
    if a is None:
        a = op(i - EX, k - EZ)
        if a is None:
            a = op(k - EZ, j - EY)
            if a is None:
                return 1.0, 0.0, 0.0, 0.0
    a /= np.linalg.norm(a)
    x0 = EX - a[0] * a
    i0 = i - a[0] * a
    if np.linalg.norm(i0) < 0.1:
        i0 = j - a[1] * a
        x0 = EY - a[1] * a

    i0 /= np.linalg.norm(i0)
    x0 /= np.linalg.norm(x0)
    cosine = np.dot(i0, x0)
    if cosine < -1.0 or cosine > 1.0:
        cosh = 0.0
        sinh = 1.0
    else:
        cosh = sqrt((1.0 + cosine) * 0.5)
        sinh = sqrt(1.0 - cosh * cosh)
    o = op(i0, x0, False)
    if np.dot(o, a) < 0.0:
        sinh = -sinh

    return np.array((cosh, -sinh * a[0], +sinh * a[1], -sinh * a[2]))


def rotmat2quat_orig(m):
    n = m.transpose()
    return rotmat2quat0(np.array(n[0]), np.array(n[1]), np.array(n[2]))

def rotmat2quat(m):
    """
    In mathematics, vectors are regarded as vertical matrices, so rotation matrices are applied from the front, but in a computer, vectors are easier to handle as horizontal matrices. In that case, the rotation matrix is multiplied from right by transposing the normal one. 

    Return the rotation matrix in posterior form to quaternion.
    """
    return rotmat2quat0(np.array(m[0]), np.array(m[1]), np.array(m[2]))



def quat2euler(q):
    """
    Qaternion to Euler angle
    """
    e = np.zeros(3)
    if q[0] == 1.0:
        return e
    p = 2 * (q[0]**2 + q[3]**2) - 1
    if p > 1.0:
        p = 1.0
    if p < -1.0:
        p = -1.0
    e[0] = acos(p)
    thh = e[0] / 2.0
    sinthh = sin(thh)
    costhh = cos(thh)
    p = q[0] / costhh
    if p > 1.0:
        p = 1.0
    if p < -1.0:
        p = -1.0
    p = acos(p)
    if sinthh == 0.0:
        s = 1.0
    else:
        s = q[1] / sinthh
    if s > 1.0:
        s = 1.0
    if s < -1.0:
        s = -1.0
    s = acos(s)
    if q[3] < 0.0:
        p = 2 * pi - p
    if q[2] > 0:
        e[2] = p + s
        e[1] = p - s
    else:
        e[2] = p - s
        e[1] = p + s
    e[1:3] %= (2 * pi)
    return e




def quat2rotmat(q):
    """
    Quaternion to Rotation matrix
    """
    a, b, c, d = q
    sp11 = (a * a + b * b - (c * c + d * d))
    sp12 = -2.0 * (a * d + b * c)
    sp13 = 2.0 * (b * d - a * c)
    sp21 = 2.0 * (a * d - b * c)
    sp22 = a * a + c * c - (b * b + d * d)
    sp23 = -2.0 * (a * b + c * d)
    sp31 = 2.0 * (a * c + b * d)
    sp32 = 2.0 * (a * b - c * d)
    sp33 = a * a + d * d - (b * b + c * c)
    return np.array([[sp11, sp12, sp13], [sp21, sp22, sp23], [sp31, sp32, sp33]]).T




def euler2quat(e):
    """
    Euler angle to quaternion
    """
    ea, eb, ec = e
    a = cos(ea / 2) * cos((ec + eb) / 2)
    b = sin(ea / 2) * cos((ec - eb) / 2)
    c = sin(ea / 2) * sin((ec - eb) / 2)
    d = cos(ea / 2) * sin((ec + eb) / 2)
    return np.array((a, b, c, d))




def accum2(comeus, cellmat, sources=None, targets=None, maxdist=13.0, exclude=None):
    """
    相互作用を近い順に加算するのみ。
    """
    water   = tip4picesites()
    charges = tip4picecharges()
    AA,BB   = tip4piceLJ()
    celli   = np.linalg.inv(cellmat)
    unitcell = cellmat / 16

    if sources is None:
        sources = np.arange(len(comeus), dtype=int)
    if targets is None:
        targets = np.arange(len(comeus), dtype=int)

    if exclude is not None:
        rcom = comeus[:,:3] @ np.linalg.inv(cellmat)

    Nmol = comeus.shape[0]
    atoms = np.zeros([Nmol,4,3])

    d_e = []
    Ro = np.zeros([Nmol,3,3])

    for i in range(Nmol):
        e = comeus[i,3:6]
        R = quat2rotmat(euler2quat(e))
        Ro[i] = regulate_ori(R)
        atoms[i] = (water @ R) + comeus[i,:3]

    for i in sources:
        if exclude is not None:
            # targetsから、特定の格子点にいる分子のみ除外する。これはかなり面倒だが、必要。
            newtargets = []
            for t in targets:
                pos = rcom[t] - rcom[i]
                pos -= np.floor(pos+0.5)
                pos = pos @ cellmat @ np.linalg.inv(Ro[i])
                # print(pos @ R45 @ np.linalg.inv(unitcell))
                grid = np.floor(pos @ R45 @ np.linalg.inv(unitcell)+0.5).astype(int)
                assert (grid[0]%2==0 and
                        grid[1]%2==0 and
                        grid[2]%2==0) or (grid[0]%2==1 and
                                          grid[1]%2==1 and
                                          grid[2]%2==1), grid
                tgrid = tuple(grid)
                if tgrid in exclude:
                    pass
                else:
                    newtargets.append(t)
            newtargets = np.array(newtargets)
        else:
            newtargets = targets

        # displacement vectors
        #print(targets, "targets")
        disp = comeus[newtargets,:3] - comeus[i,:3]
        # 簡潔さと読み易さが両立しないのでループで書く。
        wrap = np.zeros_like(disp)
        for j in range(disp.shape[0]):
            wrap[j] = np.floor(disp[j] @ celli + 0.5) @ cellmat
        dist = np.linalg.norm(disp-wrap, axis=1)
        js = np.argsort(dist)
        last = np.count_nonzero(dist < maxdist)
        js = js[:last]

        Coulombs, LJs = interactions4px(atoms[i], atoms[newtargets[js]], wrap[js], charges, AA, BB)

        sorteddist = dist[js]
        sortedEpot = Coulombs + LJs
        if np.isnan(sortedEpot[0]):
            sortedEpot[0] = 0 # avoid NaN
        sortedcumsum = np.cumsum(sortedEpot)
        d_e.append([sorteddist, sortedcumsum])
    return d_e


def accum0(comeus, cellmat, targets=None, maxdist=13.0, atoms=None, LJ=tip4piceLJ()):
    """
    相互作用を近い順に加算するのみ。

    2020-11 atomsが与えられた時もcomeusは重心位置として参照する。
    """
    charges = tip4picecharges()
    AA,BB   = LJ
    celli   = np.linalg.inv(cellmat)


    if atoms is None:
        water   = tip4picesites()
        atoms = np.zeros([comeus.shape[0],4,3])
        for i in range(comeus.shape[0]):
            e = comeus[i,3:6]
            R = quat2rotmat(euler2quat(e))
            atoms[i] = (water @ R) + comeus[i,:3]

    Nmol = atoms.shape[0]
    if targets is None:
        targets = range(Nmol)

    d_e = []

    for i in targets:
        # displacement vectors
        disp = comeus[:,:3] - comeus[i,:3]
        wrap = np.zeros_like(disp)
        for j in range(disp.shape[0]):
            wrap[j] = np.floor(disp[j] @ celli + 0.5) @ cellmat
        dist = np.linalg.norm(disp-wrap, axis=1)
        js = np.argsort(dist)
        last = np.count_nonzero(dist < maxdist)
        js = js[:last]

        Coulombs, LJs = interactions4p(0, atoms[js], wrap[js], charges, AA, BB)

        sorteddist = dist[js]
        sortedEpot = Coulombs + LJs

        sortedEpot[0] = 0 # avoid NaN
        sortedcumsum = np.cumsum(sortedEpot)
        d_e.append([sorteddist, sortedcumsum])
    return d_e


def accum_sd(comeus, cellmat, targets, r1=3.5, lw=1, tag="", style='-'):
    """
    相互作用を近い順に加算して、距離ごとの、累積対相互作用の標準偏差を計算する。
    """
    d_e = accum0(comeus, cellmat, targets)

    equalspacing = np.linspace(0, 13.0, 250)
    css = np.zeros_like(equalspacing)
    csss = np.zeros_like(equalspacing)
    cnt = 0

    for sorteddist, sortedcumsum in d_e:
        cs = np.zeros_like(equalspacing)
        for i in range(equalspacing.shape[0]):
            cut = sortedcumsum[sorteddist < equalspacing[i]]
            if cut.shape[0] == 0:
                cs[i] = 0
            else:
                cs[i] = cut[-1] # last value
        css += cs
        csss += cs**2
        cnt += 1
    #avg
    css /= cnt
    # sq avg
    csss /= cnt
    # SD
    sd = (csss - css**2)**0.5
    plt.plot(equalspacing, sd, style, label="{1} water {0}".format(i,tag), lw=lw)
    print(np.interp(r1, equalspacing, sd) / np.interp(13, equalspacing, sd), " convergence")

def accum_sd2(comeus, cellmat, sources, targets, r1=3.5, maxdist=13.0, **kwarg):
    """
    相互作用を近い順に加算する。
    距離ごとの、積算対相互作用の標準偏差を計算する。
    """
    d_e = accum2(comeus, cellmat, sources=sources, targets=targets, maxdist=maxdist)

    equalspacing = np.linspace(0, 13.0, 1000)
    css = np.zeros_like(equalspacing)
    csss = np.zeros_like(equalspacing)
    cnt = 0

    for sorteddist, sortedcumsum in d_e:
        # 距離を等間隔にとりなおす。少し難しい。
        cs = np.zeros_like(equalspacing)
        for i in range(equalspacing.shape[0]):
            cut = sortedcumsum[sorteddist < equalspacing[i]]
            if cut.shape[0] == 0:
                cs[i] = 0
            else:
                cs[i] = cut[-1] # last value
        css += cs
        csss += cs**2
        cnt += 1
    #avg
    css /= cnt
    # sq avg
    csss /= cnt
    # SD
    sd = (csss - css**2)**0.5
    plt.plot(equalspacing, sd, **kwarg)
    print(np.interp(r1, equalspacing, sd) / np.interp(13, equalspacing, sd), " convergence")



# ネットワークを分ける。結晶の格子点を参考にするのが一番確実。
# ただし、位置だけでは判別できない。まじめに水素結合ネットワークを分けるのか、それとも?

def primseco(comeus, cellmat):
    assert len(comeus) == 1024
    unitcell = cellmat / 16
    pos = comeus[:,:3]
    grid = np.floor(pos @ np.linalg.inv(unitcell)+0.5).astype(int)
    grid %= 16
    isprim = [None for i in range(len(comeus))]
    for i,g in enumerate(grid):
        x,y,z = g
        if x%2 == 1: #1つが奇数なら全部奇数
            x-=1
            y-=1
            z-=1
        isprim[i] = (x+y+z) % 4 == 0
    return np.array(isprim, dtype=bool)



def nodes2rotmat(c,i,j,rpos,cellmat):
    h1 = rpos[i] - rpos[c]
    h1 -= np.floor(h1+0.5)
    h2 = rpos[j] - rpos[c]
    h2 -= np.floor(h2+0.5)
    h1 = h1 @ cellmat
    h2 = h2 @ cellmat
    z = h1+h2
    y = h1-h2
    x = np.cross(y,z)
    x /= np.linalg.norm(x)
    y /= np.linalg.norm(y)
    z /= np.linalg.norm(z)
    return np.array([x,y,z])

def hbn(rcom, cellmat, R, water, icetest=True):
    """
    分子配置から有向グラフを構成。
    属性として、どっちの水素かを記録しておく。
    相対ベクトルを保存する。
    """
    dg = nx.DiGraph()
    for i, j in pl.pairs_iter(rcom, 3.2, cellmat, distance=False):
        rij = rcom[j] - rcom[i]
        rij -= np.floor(rij + 0.5)
        dij = rij @ cellmat
        wi = water @ R[i]
        wj = water @ R[j] + dij
        dmin = 2.3**2 # ice 16のために長めにした。
        hb   = None
        rij2 = rij.copy()
        for k in (0,1):
            dOH = wi[k] - wj[2]
            L2 = dOH @ dOH
            if L2 < dmin:
                dmin = L2
                hb   = (i, j, k) # i (hydrogen is k) to j
            dOH = wj[k] - wi[2]
            L2 = dOH @ dOH
            if L2 < dmin:
                dmin = L2
                hb   = (j, i, k) # j (hydrogen is k) to i
                rij2 = -rij
        if hb is not None:
            dg.add_edge(hb[0], hb[1], H=hb[2], vec=rij2)
    if icetest:
        for x in dg.out_degree(dg):
            assert x[1]==2, x
    return dg


def cycle_energy(cycle, center, rcom, cellmat, waters, dg):
    # 酸素の座標は一括抽出できる。
    # 酸素の電荷はO位置ではなくM位置(3番目)
    oxygens = waters[cycle, 3]

    # cycleの中にcenterを含む場合はoxygensとhydrogensからそれらを除去しておく必要がある。
    if center in cycle:
        oxygens = waters[cycle[cycle!=center], 3]

    # 水素はどっちか判別必要。
    which = []
    N = len(cycle)
    #print(cycle)
    for i, w in enumerate(cycle):
        nxt = cycle[(i+1) % N]
        which.append(dg[w][nxt]["H"])
    hydrogens = waters[cycle, which]

    if center in cycle:
        which = np.array(which)
        hydrogens = waters[cycle[cycle!=center], which[cycle!=center]]

    # 分子の相対位置
    d = rcom[cycle] - rcom[center]
    d -= np.floor(d + 0.5)
    d = d @ cellmat

    morder = cycle[np.argsort(np.linalg.norm(d, axis=1))]
    dmin = np.min(np.linalg.norm(d, axis=1))
    # print(np.linalg.norm(d, axis=1))
    # coulomb相互作用の和
    ec = 0
    charges = [1.0, 1.0, 0.0, -2.0] # 中心分子の電荷。中心分子は半分ではない。
    if center in cycle:
        d = d[cycle!=center]
    for i in (0,1,3):
        # 中心分子の原子位置
        ac = waters[center, i]

        # cycle上の半酸素位置との差
        rel = oxygens - ac + d
        L   = np.linalg.norm(rel, axis=1)
        ec -= np.sum(1/L) * charges[i]

        # cycle上の水素位置との差
        rel = hydrogens - ac + d
        L   = np.linalg.norm(rel, axis=1)
        ec += np.sum(1/L) * charges[i]

    unitcharge = tip4picecharges()[0]
    return ec * unitcharge**2 * qe0, dmin, list(morder)



from scipy.stats import norm

# 正規分布をたたみこんでスパイクをなくす。
def smooth(sd):
    Nf = 7
    filter = norm.pdf(range(-Nf,Nf+1),0,Nf)
    sd2= np.convolve(sd, filter)
    N = len(sd)
    sd2 = sd2[Nf:N+Nf]
    # print(len(sd2), len(sd))
    return sd2
