#!/usr/bin/env python3

import numpy as np
import os
import sys
import shutil
import glob
import argparse
import subprocess
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--test", action="store_true")
parser.add_argument("--namespace", default="")

# KEYWORDS
# $$MODEL$$ - constituitiveProperties model and parameters
# $$TIME$$ - end time of inlet turn-on
# $$UMAX$$ - maximum inlet velocity
# $$WALLV$$ - wall velocity
# $$RADIUS$$ - radius of the channel
# $$NBLOCKS$$ - number of blocks in the channel for parallelization

args = parser.parse_args()

TEST = args.test
OVERWRITE = args.overwrite
NAMESPACE = args.namespace
if NAMESPACE != "":
    NAMESPACE = f"{NAMESPACE}/"

if TEST:
    path = f"test/{NAMESPACE}model/*/wi/*/wallv/*/dist/*"
else:
    path = f"prod/{NAMESPACE}model/*/wi/*/wallv/*/dist/*"

runs = glob.glob(path)
runs = [os.path.abspath(run) for run in runs]

for run in runs:
    os.chdir(run)

    marker = pathlib.Path("done")
    if marker.exists() and not OVERWRITE:
        print("Skipping", run)
        continue

    print(run)

    command = "./Allclean & ./Allrun"
    rc = os.system(command)

    idx = 0
    while rc != 0:
        if idx == 5:
            break
        process = os.system(command)
        idx += 1
    if idx == 5:
        continue

    # touch file
    marker.touch()
    
