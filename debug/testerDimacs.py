from os.path import isfile, join
import os
from pathlib import Path
import subprocess
import time

file_path = list()
current_path = str(Path(__file__).parent.absolute()) 

#input_path = current_path + "/input/DIMACS/"
input_path = current_path + "/input/DIMACS-V/"
#input_path = current_path + "/input/DIMACS-U/"

if not os.path.exists(input_path):
    os.mkdir(input_path)

output_path = current_path + "/output/"

if not os.path.exists(output_path):
    os.mkdir(output_path)
else:
    os.remove(output_path + "/DIMACS")


for r, d, f in os.walk(input_path):
    for file in f:
        file_path.append([file, os.path.join(r,file)])

file_path.sort(key=lambda x: (len(x[0]), x[0]))

for (x, y) in file_path:
    print("Testing: " + x)
    f = open(output_path + "/DIMACS", 'a+', buffering=1)
    cmd = str(current_path + "/main" +  " < " +  y)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #Send it to run
    pStartTime = time.time()
    (output, err) = process.communicate()
    pEndTime = time.time()
    #Now store values
    f.write(x + "\t" + str(output.decode('UTF-8')))
    f.flush()
    f.close()
