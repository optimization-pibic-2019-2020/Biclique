import os
from pathlib import Path
import random

current_path = str(Path(__file__).parent.absolute())

path = current_path + "../input/Barabasi-Albert/"
if not os.path.exists(path):
    os.mkdir(path)

E = list()
for i in range(128,256):
    E += [i]
R = [(2, "/2-clique/")]

for e in E:
    s = open(path + str(e), 'w', buffering=1)
    vertices = list()
    degree = [0] * e
    edges = list()
    vertices += [0, 1]
    edges.append([0, 1])
    degree[0] += 1
    degree[1] += 1
    created = 2
    for u in range(created, e):
        summation = 0
        for v in range(0, created):
            summation += degree[v]
        for v in range(0, created):
            p = degree[v] / summation
            if random.uniform(0, 1) < p + 1e-6:
                edges.append([u, v])
                degree[u] += 1
                degree[v] += 1
        created += 1
    vertex_count = e
    edge_count = len(edges)
    s.write(str(vertex_count) + " " + str(edge_count) + "\n")
    for pair in edges:
        s.write('{} {}'.format(pair[0]+1, pair[1]+1) + "\n")
    s.close()

