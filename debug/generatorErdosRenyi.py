import os
from pathlib import Path
import random
import numpy

current_path = str(Path(__file__).parent.absolute())

path = current_path + "/input/Erdos-Renyi"
if not os.path.exists(path):
    os.mkdir(path)

T = list()
for i in numpy.arange(0.01, 0.99, 0.01):
    T.append(['/{}/'.format(round(i, 2)), i]);

E = list()
for i in range(128,129):
    E += [i];
R = [4]

for (t,p) in T:
    if not os.path.exists(path + t):
        os.mkdir(path + t)
    for e in E:
        s = open(path + t + str(e), 'w', buffering=1)
        edges = list()
        for u in range(e):
            for v in range(u+1, e):
                if random.uniform(0, 1) <= p:
                    edges.append([u,v])
        vertex_count = e
        edge_count = len(edges)
        s.write(str(vertex_count) + " " + str(edge_count)  + "\n")
        for pair in edges:
            s.write('{} {}'.format(pair[0]+1, pair[1]+1) + "\n")
        s.close()


