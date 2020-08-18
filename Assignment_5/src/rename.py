#!/usr/bin/python3

import sys
import shutil

if len(sys.argv) == 1:
    print("This script sets up the mini demo for the rasterizer.")
    print("The following demos are supported: \nbase, lines, attributes, depth, view, blend, animation")
else:
    shutil.copyfile("extra/attributes_" + str(sys.argv[1]) + ".h","attributes.h")
    shutil.copyfile("extra/main_" + str(sys.argv[1]) + ".cpp","main.cpp")