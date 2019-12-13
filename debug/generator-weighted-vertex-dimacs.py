import sys
import os
import time
from os import listdir
from os.path import isfile, join
from pathlib import Path
from datetime import datetime
import random

file_path = list()
current_path = str(Path(__file__).parent.absolute())

path = current_path + "/input/DIMACS/"

path_v = current_path + "/input/DIMACS_V/"
if not os.path.exists(path_v):
    os.mkdir(path_v)

for r, d, f in os.walk(path):
    for file in f:
        if '.mtx' in file:
            file_path.append([file, os.path.join(r,file)])

file_path.sort()

T = [(2,"/2-clique/")]

for (t,d) in T:
    for (x, y) in file_path:
        p = os.path.basename(os.path.dirname(y))
        if not os.path.exists(path_v + p + "/"):
            os.mkdir(path_v + p + "/")
        if not os.path.exists(path_v + p + d):
            os.mkdir(path_v + p + d)
        f = open(y, 'r', buffering=1)
        s = open(path_v + p + d + x + str(-t), 'w', buffering=1)
        vertex_count, edge_count = [int(z) for z in next(f).split()]
        edges = []
        edges = [[int(z) for z in line.split()] for line in f]
        s.write(str(vertex_count) + " " + str(edge_count) + "\n")
        for i in range(vertex_count):
            s.write(str((i%200)+1) + '\n')
        for edge_list in edges:
            for i in range(len(edge_list)):
                s.write(str(edge_list[i]) + " ")
            s.write("\n")
        f.close()
        s.close()

