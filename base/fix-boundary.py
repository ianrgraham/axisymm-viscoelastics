#!/usr/bin/env python3

import sys
import os
import shutil
import subprocess
import glob

file = "constant/polyMesh/boundary"

msh_file = glob.glob("msh/2rbc-channel_scale-1.4*.msh")[0]  # error if none found
subprocess.run(["gmshToFoam", msh_file])

shutil.rmtree("constant/polyMesh.inc", ignore_errors=True)

wedges = ["front", "back"]
walls = ["rbc", "rbc2", "wall"]

# Read the input file
with open('constant/polyMesh/boundary', 'r') as f:
    lines = f.readlines()

state = None

with open('constant/polyMesh/boundary', 'w') as f:
    for line in lines:

        strip_line = line.strip()

        if state is None:
            if strip_line in wedges:
                state = "wedge"
            elif strip_line in walls:
                state = "wall"
        elif state == "wedge" and "type" in strip_line:
            line = line.replace('patch', 'wedge')
        elif state == "wall" and "type" in strip_line:
            line = line.replace('patch', 'wall')

        if "}" in strip_line:
            state = None

        f.write(line)

shutil.copytree("constant/polyMesh", "constant/polyMesh.inc")