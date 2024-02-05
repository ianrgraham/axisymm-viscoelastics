#!/usr/bin/env python3
#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500
#SBATCH --ntasks=1
#SBATCH --job-name=rbve-final2
#SBATCH --output="slurm-%A_%a.out"
#SBATCH --error="slurm-%A_%a.err"
#SBATCH --array=0-271
#SBATCH --partition=p_rrig

import numpy as np
import os
import sys
import shutil
import glob
import argparse
import subprocess
import pathlib

# parser = argparse.ArgumentParser()
# parser.add_argument("--overwrite", action="store_true")
# parser.add_argument("--test", action="store_true")
# parser.add_argument("--namespace", default="")

# KEYWORDS
# $$MODEL$$ - constituitiveProperties model and parameters
# $$TIME$$ - end time of inlet turn-on
# $$UMAX$$ - maximum inlet velocity
# $$WALLV$$ - wall velocity
# $$RADIUS$$ - radius of the channel
# $$NBLOCKS$$ - number of blocks in the channel for parallelization

# args = parser.parse_args()

# TEST = args.test
# OVERWRITE = args.overwrite
# NAMESPACE = args.namespace
TEST = False
OVERWRITE = False
NAMESPACE = "final2"
if NAMESPACE != "":
    NAMESPACE = f"{NAMESPACE}/"

if TEST:
    path = f"test/{NAMESPACE}model/*/wi/*/wallv/*/dist/*"
else:
    path = f"prod/{NAMESPACE}model/*/wi/*/wallv/*/dist/*"

runs = glob.glob(path)
runs = sorted([os.path.abspath(run) for run in runs])

print(len(runs))
print(runs[0])
task_id = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(task_id)

run = runs[task_id]
# for run in runs:
os.chdir(run)

marker = pathlib.Path("done")
running = pathlib.Path("running")
if marker.exists() and not OVERWRITE:
    print("Skipping", run)
    sys.exit()

print(run)

command = "module load gcc/7.3.0 ; source $OF9DIR/etc/bashrc ; ./Allclean ; ./Allrun"
rc = os.system(command)

idx = 0
while rc != 0:
    if idx == 10:
        break
    process = os.system("./Allclean ; ./Allrun")
    idx += 1
if idx >= 10:
    sys.exit()

# touch file
marker.touch()

