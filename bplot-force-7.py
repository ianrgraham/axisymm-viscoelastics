#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import glob
import os
import argparse
# import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--test", action="store_true")
parser.add_argument("--namespace", default="")
parser.add_argument("--namespaces", default="", nargs="+")
parser.add_argument("--wi", default="1.0")
parser.add_argument("--model", default="giesekus")

args = parser.parse_args()

TEST = args.test
NAMESPACE = args.namespace
WI = args.wi
if NAMESPACE != "":
    NAMESPACE = f"{NAMESPACE}/"
print(NAMESPACE)

NAMESPACES = args.namespaces
if len(NAMESPACES) > 0:
    assert NAMESPACE == ""
    NAMESPACES = [f"{ns}/" for ns in NAMESPACES]
else:
    NAMESPACES = [NAMESPACE]
print(NAMESPACES)

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

markers = ["s", "o", "x",]
# wi_num = WI
# model = args.model

plt.figure(figsize=(12, 8))

for NAMESPACE in NAMESPACES:
    m = markers.pop()
    if TEST:
        path = f"test/{NAMESPACE}model/*/wi/*/wallv/*/dist"
    else:
        path = f"prod/{NAMESPACE}model/*/wi/*/wallv/*/dist"

    runs = glob.glob(path)
    cmap = cm.magma
    norm = colors.LogNorm(vmin=0.01, vmax=20)

    print(len(runs))

    for run in runs:

        model = run.split("model/")[-1].split("/")[0]

        if "fene-p" in run:
            ls = ":"
        else:
            ls = "-"

        wi = float(run.split("/")[-4])

        if wi > 6.0:
            continue
            
        xs = []
        ys = []
        yas = []
        ybs = []
        ycs = []
        for run in glob.glob(f"{run}/*"):
            try:
                # if not os.path.exists(f"{run}/done"):
                    # continue
                wallv = run.split("/")[-1]
                wallv = float(wallv)
                # print(wallv)
                c = "k" # cmap(norm(wallv))
                time, rbc1_press, rbc1_visc, rbc1_poly = read_force_data(f"{run}/postProcessing/forces0/0/rheoForces.dat")

                rbc1_force = (rbc1_press - rbc1_visc - rbc1_poly)
                # r1d = (rbc1_press - rbc1_visc)

                time, rbc2_press, rbc2_visc, rbc2_poly = read_force_data(f"{run}/postProcessing/forces1/0/rheoForces.dat")

                rbc2_force = (rbc2_press - rbc2_visc - rbc2_poly)
                # r2d = (rbc1_press - rbc1_visc)
                diff = (rbc2_force - rbc1_force)
                # print(wallv, len(diff))
                # plt.plot(wallv, diff[-1], "o", color=c)
                # plt.plot(wallv, rbc1_part_force[-1], "s", color=c)
                # plt.plot(rbc1_force, color=c)
                # print(rbc1_force[-1])
                xs.append(wallv)
                ys.append(diff[-1])
                yas.append(rbc2_press[-1] - rbc1_press[-1])
                ybs.append(rbc1_visc[-1] - rbc2_visc[-1])
                ycs.append(rbc1_poly[-1] - rbc2_poly[-1])

                    
            except:
                pass

        if len(xs) > 0:
            print()
            xs, ys, yas, ybs, ycs = zip(*sorted(zip(xs, ys, yas, ybs, ycs)))
            xs = np.array(xs)
            ys = np.array(ys)
            yas = np.array(yas)
            ybs = np.array(ybs)
            ycs = np.array(ycs)
            print(len(xs), len(ys))
            print(xs)
            print(ys)

            # plt.plot(xs, ys, ls, color=cmap(norm(wi)), label=f"{model} {wi}")
            plt.plot(xs, yas, ls, marker="D", color=cmap(norm(wi)), label=f"{model} {wi} press")
            plt.plot(xs, ybs, ls, marker="x", color=cmap(norm(wi)), label=f"{model} {wi} visc")
            plt.plot(xs, ycs, ls, marker=".", color=cmap(norm(wi)), label=f"{model} {wi} poly")
            plt.plot(xs, ys, ls, marker="s", color=cmap(norm(wi)), label=f"{model} {wi} total")
        



plt.legend(title=r"model $\mathrm{Wi}$", loc="upper right")
# plt.legend(title="sample", loc="upper right")
# plt.legend(["Press", "Visc", "Poly", "Total", "Press + Visc"])
# plt.xscale("log")
plt.ylabel(r"$F_x$", size="x-large")
plt.xlabel(r"$\Delta x$", size="x-large")

# plt.vlines(0.0, *plt.ylim(), linestyle="--", color="k")
# plt.hlines(0.0, *plt.xlim(), linestyle="--", color="k")

# plt.ylim(-2e-7, 2e-7)

plt.savefig("test-7.png", dpi=300)

