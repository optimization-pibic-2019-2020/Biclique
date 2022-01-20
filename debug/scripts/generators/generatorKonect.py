import sys
import os
import time
from os import listdir
from os.path import isfile, join
from pathlib import Path
from datetime import datetime

files = list()

absolute = str(Path(__file__).parent.absolute())

folders = [(absolute+"../input/KONECT/", absolute+"../input/KONECT-V/")]

for (inputPath,outputPath) in folders:
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    files.clear()
    for r, d, f in os.walk(inputPath):
        for file in f:
            files.append([file, os.path.join(r,file)])
    files.sort(key=lambda x: (len(x[0]), x[0]))
    for (x, y) in files:
        p = os.path.basename(os.path.dirname(y))
        if not os.path.exists(outputPath + p + "/"):
            os.mkdir(outputPath + p + "/")
        f = open(y, 'r', buffering=1)
        s = open(outputPath + p + "/" + x, 'w', buffering=1)
        partitionA_count, partitionB_count, edge_count = [int(z) for z in next(f).split()]
        edges = []
        edges = [[int(z) for z in line.split(",")] for line in f]
        s.write(str(2) + " " + str(partitionA_count) + " " + str(partitionB_count) + " " + str(edge_count) + "\n")
        for i in range(partitionA_count):
            s.write(str((i%200)+1) + '\n')
        for j in range(partitionB_count):
            s.write(str((j%200)+1) + '\n')
        for edge_list in edges:
            for i in range(len(edge_list)):
                if i == 0:
                    s.write(str(edge_list[i]) + " ")
                else:
                    s.write(str(edge_list[i] + partitionA_count))

            s.write("\n")
        f.close()
        s.close()





