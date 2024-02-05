#!/usr/bin/env python3

import numpy as np
import os
import sys
import glob
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--lambda", type=float, default=0.1, dest="_lambda")
parser.add_argument("--poly-visc", type=float, default=4e-9) # 0.3e-6
parser.add_argument("--beta", type=float, default=0.2)
parser.add_argument("--radius", type=float, default=3.7)
parser.add_argument("--ramp-time", type=float, default=1.0)
parser.add_argument("--total-time", type=float, default=6.0)
parser.add_argument("--nblocks", type=int, default=7)
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--test", action="store_true")
parser.add_argument("--namespace", default="")
parser.add_argument("--rel-wallv", type=float, default=None)

# KEYWORDS
# $$MODEL$$ - constituitiveProperties model and parameters
# $$TIME$$ - end time of inlet turn-on
# $$UMAX$$ - maximum inlet velocity
# $$WALLV$$ - wall velocity
# $$RADIUS$$ - radius of the channel
# $$NBLOCKS$$ - number of blocks in the channel for parallelization

args = parser.parse_args()

LAMBDA = args._lambda
RADIUS = args.radius # um
RAMPTIME = args.ramp_time # s
TOTALTIME = args.total_time # s
NBLOCKS = args.nblocks
POLYVISC = args.poly_visc 
BETA = args.beta
SOLVVISC = BETA * POLYVISC / (1.0 - BETA)
OVERWRITE = args.overwrite
TEST = args.test
NAMESPACE = args.namespace
REL_WALLV = args.rel_wallv
if NAMESPACE != "":
    NAMESPACE = f"{NAMESPACE}/"

def replace_keywords(path, keys_and_subs):
    with open(path, "r") as f:
        lines = f.readlines()
    with open(path, "w") as f:
        for line in lines:
            for (keyword, sub) in keys_and_subs:
                if keyword in line:
                    line = line.replace(keyword, sub)
            f.write(line)

models = [
("fene-p", f"""
    type             FENE-PLog;
    
    rho              rho [1 -3 0 0 0 0 0] 1.0e-18;
    etaS             etaS [1 -1 -1 0 0 0 0] {SOLVVISC:e};
    etaP             etaP [1 -1 -1 0 0 0 0] {POLYVISC:e};            
    lambda           lambda [0 0 1 0 0 0 0] {LAMBDA};
    L2               L2 [0 0 0 0 0 0 0] 5;
    solveInTau       false;
    
    stabilization    coupling;
"""),
("giesekus", f"""
    type             GiesekusLog;
    
    rho              rho [1 -3 0 0 0 0 0] 1.0e-18;
    etaS             etaS [1 -1 -1 0 0 0 0] {SOLVVISC:e};
    etaP             etaP [1 -1 -1 0 0 0 0] {POLYVISC:e};             
    lambda           lambda [0 0 1 0 0 0 0] {LAMBDA};
    alpha            alpha [0 0 0 0 0 0 0] 0.4;
    
    stabilization    coupling;
""")
]

# ,
# ("oldroyd-b", f"""
#     type             Oldroyd-BLog;
    
#     rho              rho [1 -3 0 0 0 0 0] 1.0e-18;
#     etaS             etaS [1 -1 -1 0 0 0 0] {SOLVVISC:e};
#     etaP             etaP [1 -1 -1 0 0 0 0] {POLYVISC:e};            
#     lambda           lambda [0 0 1 0 0 0 0] {LAMBDA};
    
#     stabilization    coupling;
# """)


wallvs = {
    ("giesekus", 0.1): -3.1396827418827233,
    ("giesekus", 1.0): -29.803194236265306,
    ("giesekus", 5.0): -150.9941630111803,
    ("fene-p", 0.1): -3.144287387208624,
    ("fene-p", 1.0): -31.01366875007078,
    ("fene-p", 5.0): -150.57528254493394
}

mshs = sorted(glob.glob("msh/*.msh"), key=lambda x: float(x.split("dist-")[-1].split(".msh")[0]))

if TEST:
    mshs = mshs[:1]

if TEST:
    models = models[:1] # [1:] # just giesekus
    wis = [5.0]
else:
    wis = [0.1, 1.0, 5.0, 10.0]


for (model, model_sub) in models:
    for wi in wis:
        umax = wi * RADIUS / LAMBDA

        for msh in mshs:
            if REL_WALLV is not None:
                assert REL_WALLV >= 0.0
                wallv = -REL_WALLV * umax
            else:
                wallv = 0.0

            dist = msh.split("dist-")[-1].split(".msh")[0]

            wall_str = f"{wallv:.4f}".replace("-", "n")

            # copy base into new folder
            if TEST:
                path = f"test/{NAMESPACE}model/{model}/wi/{wi:.1f}/wallv/{wall_str}/dist/{dist}"
            else:
                path = f"prod/{NAMESPACE}model/{model}/wi/{wi:.1f}/wallv/{wall_str}/dist/{dist}"

            if os.path.exists(path):
                if OVERWRITE:
                    shutil.rmtree(path)
                else:
                    continue

            shutil.copytree("base", path)

            # copy msh file
            shutil.copy(msh, f"{path}/msh/")

            # modify system/constitutiveProperties
            keys_and_subs = [("$$MODEL$$", model_sub)]
            replace_keywords(f"{path}/constant/constitutiveProperties", keys_and_subs)

            # modify 0/U
            keys_and_subs = [
                ("$$TIME$$", str(RAMPTIME)),
                ("$$UMAX$$", str(umax)),
                ("$$WALLV$$", str(wallv)),
                ("$$RADIUS$$", str(RADIUS))
            ]
            replace_keywords(f"{path}/0/U", keys_and_subs)

            # modify system/controlDict
            keys_and_subs = [
                ("$$TOTALTIME$$", str(TOTALTIME))
            ]
            replace_keywords(f"{path}/system/controlDict", keys_and_subs)


            # modify system/decomposeParDict
            keys_and_subs = [
                ("$$NBLOCKS$$", str(NBLOCKS))
            ]
            replace_keywords(f"{path}/system/decomposeParDict", keys_and_subs)
