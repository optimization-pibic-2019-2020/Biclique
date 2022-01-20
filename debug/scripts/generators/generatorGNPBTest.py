import os
from pathlib import Path
import random

current_path = str(Path(__file__).parent.absolute())

path = current_path + "../input/GNPBTest-V/"
if not os.path.exists(path):
    os.mkdir(path)


p = 0.050
seed = 1

while(p <= 1.800):
    index = 1
    while(index <= 10):
        s = open(path + format(p, '.3f') + '-' + str(index), 'w', buffering=1)
        random.seed(seed)
        partition_A_vertices = list(range(0, 500))
        partition_B_vertices = list(range(500, 1000))
        edges = list()
        for v in partition_A_vertices:
            for u in partition_B_vertices:
                probability_res = random.choices([0, 1], [p, 100.0 - p])
                if probability_res[0] == 0 and v != u:
                    if not ([u, v] or [v, u]) in edges:
                        edges.append([v, u])
        partition_A_count = len(partition_A_vertices)
        partition_B_count = len(partition_B_vertices)
        vertices = partition_A_count + partition_B_count
        edge_count = len(edges)
        s.write("2 " + str(partition_A_count) + " " + str(partition_B_count) + " " + str(edge_count) + "\n")
        for v in range(0, vertices):
            s.write('{}'.format((v % 200) + 1) + "\n")
        for pair in edges:
            s.write('{} {}'.format(pair[0]+1, pair[1]+1) + "\n")
        s.close()
        seed += 1
        index += 1
    p += 0.001
    p = round(p, 3)
    

