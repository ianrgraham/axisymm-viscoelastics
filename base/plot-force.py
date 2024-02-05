#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import glob
# import pandas as pd

def read_wi(path):
    with open(path) as f:
        for line in f:
            if "Wi:" in line:
                return float(line.split()[1])
    return None
        

def read_force_data(path):
    visc = []
    press = []
    poly = []
    time = []
    with open(path) as f:
        for line in f.readlines()[3:]:
            data = line.split()
            # print(last_line)
            data = [d.strip("()") for d in data]

            time.append(data[0])

            data = [data[i:i+3] for i in range(1, len(data), 3)]

            press.append(float(data[0][0]))
            visc.append(float(data[1][0]))
            poly.append(float(data[2][0]))
    return np.array(time), np.array(press), np.array(visc), np.array(poly)

def read_force_data_normal(path):
    visc = []
    press = []
    # poly = []
    time = []
    with open(path) as f:
        for line in f.readlines()[3:]:
            data = line.split()
            # print(last_line)
            data = [d.strip("()") for d in data]

            time.append(data[0])

            data = [data[i:i+3] for i in range(1, len(data), 3)]

            press.append(float(data[0][0]))
            visc.append(float(data[1][0]))
            # poly.append(float(data[2][0]))
    return np.array(time), np.array(press), np.array(visc)
time, rbc1_press, rbc1_visc, rbc1_poly = read_force_data("postProcessing/forces0/0/rheoForces.dat")

rbc1_force = (rbc1_press + rbc1_visc + rbc1_poly)
plt.plot(rbc1_force)
print(rbc1_force[-1])

try:
    wi = read_wi("log.postProcess")
    print(wi)
    plt.text(0.5, 0.5, "Wi = {:.3f}".format(wi), transform=plt.gca().transAxes)
except:
    pass


# plt.xscale("log")
plt.ylabel("Force")
plt.xlabel("Time")

plt.savefig("test.png", dpi=300)

