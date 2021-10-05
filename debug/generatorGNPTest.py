import os
from pathlib import Path
import random

current_path = str(Path(__file__).parent.absolute())

path = current_path + "/input/GNPTest-V/"
if not os.path.exists(path):
    os.mkdir(path)


p = 0.050
seed = 1
while(p <= 1.000):
    index = 1
    while(index <= 10):
        s = open(path + format(p, '.3f') + '-' + str(index), 'w', buffering=1)
        random.seed(seed)
        vertices = list(range(0, 1000))
        edges = list()
        for v in vertices:
            for u in range(v + 1, 1000):
                probability_res = random.choices([0, 1], [p, 100.0 - p])
                if probability_res[0] == 0 and v != u:
                    if not ([u, v] or [v, u]) in edges:
                        edges.append([v, u])
        vertex_count = len(vertices)
        edge_count = len(edges)
        s.write("1 " + str(vertex_count) + " " + str(edge_count) + "\n")
        for v in vertices:
            s.write('{}'.format((v % 200) + 1) + "\n")
        for pair in edges:
            s.write('{} {}'.format(pair[0]+1, pair[1]+1) + "\n")
        s.close()
        seed += 1
        index += 1
    p += 0.001
    p = round(p, 3)
    

