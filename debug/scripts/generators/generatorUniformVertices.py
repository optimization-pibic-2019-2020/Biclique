import sys
import os
import time
from os import listdir
from os.path import isfile, join
from pathlib import Path
from datetime import datetime

files = list()

absolute = str(Path(__file__).parent.absolute())

folders = [(absolute+"../input/DIMACS/", absolute+"../input/DIMACS-U/")] 
#folders += [(absolute+"../input/BHOSLIB/", absolute+"../input/BHOSLIB-U/")]
#folders += [(absolute+"../input/Barabasi-Albert/", absolute+"../input/Barabasi-Albert-U/")]
#folders += [(absolute+"../input/Erdos-Renyi/", absolute+"../input/Erdos-Renyi-U/")]

for (inputPath,outputPath) in folders:
	if not os.path.exists(outputPath):
		os.mkdir(outputPath)
	files.clear()
	for r, d, f in os.walk(inputPath):
		for file in f:
			files.append([file, os.path.join(r,file)])
	files.sort(key=lambda x: (len(x[0]), x[0]))
	folders.sort(key=lambda x: (len(x[0]), x[0]))
	for (x, y) in files:
		p = os.path.basename(os.path.dirname(y))
		if not os.path.exists(outputPath + p + "/"):
			os.mkdir(outputPath + p + "/")
		f = open(y, 'r', buffering=1)
		s = open(outputPath + p + "/" + x, 'w', buffering=1)
		vertex_count, edge_count = [int(z) for z in next(f).split()]
		edges = []
		edges = [[int(z) for z in line.split()] for line in f]
		s.write(str(vertex_count) + " " + str(edge_count) + "\n")
		for i in range(vertex_count):
			s.write(str(1) + '\n')
		for edge_list in edges:
			for i in range(len(edge_list)):
				s.write(str(edge_list[i]) + " ")
			s.write("\n")
		f.close()
		s.close()
	

