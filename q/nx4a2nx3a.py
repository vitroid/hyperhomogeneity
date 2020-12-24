#!/usr/bin/env python


from math import *

def quat2euler(q):
    e = [0.0] * 3
    if q[0] == 1.0:
        return e
    p = 2*(q[0]**2+q[3]**2) - 1
    if p>1.0:
        p=1.0
    if p<-1.0:
        p=-1.0
    e[0] = acos(p)
    thh = e[0] / 2.0
    sinthh = sin(thh)
    costhh = cos(thh)
    p = q[0]/costhh
    if p>1.0:
        p=1.0
    if p<-1.0:
        p=-1.0
    p = acos(p)
    if sinthh == 0.0:
        s = 1.0
    else:
        s = q[1]/sinthh
    if s>1.0:
        s=1.0
    if s<-1.0:
        s=-1.0
    s=acos(s)
    if q[3]<0.0:
        p=2*pi-p
    if q[2]>0:
        e[2] = p+s
        e[1] = p-s
    else:
        e[2] = p-s
        e[1] = p+s
    return e


import sys

while True:
    line = sys.stdin.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if columns[0] == '@BOX3':
        print("@BOX3")
        line = sys.stdin.readline()
        print(line, end="")
    if columns[0] == '@NX4A':
        print("@NX3A")
        line = sys.stdin.readline()
        nmol = int(line)
        print(nmol)
        for i in range(nmol):
            line = sys.stdin.readline()
            columns = line.split()
            for i in columns[0:3]:
                print(i, end=" ")
            for i in quat2euler([float(x) for x in columns[3:7]]):
                print(i, end=" ")
            print()
