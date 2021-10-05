from os.path import isfile, join, basename
from datetime import datetime
import os
from pathlib import Path
import subprocess
import time

file_path = list()
current_path = str(Path(__file__).parent.absolute()) 

# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

input_path = current_path + "/input/GNP-V/"

if not os.path.exists(input_path):
    os.mkdir(input_path)

output_path = current_path + "/logs/GNP/"

if not os.path.exists(output_path):
    os.mkdir(output_path)

for r, d, f in os.walk(input_path):
    for file in f:
        file_path.append([file, os.path.join(r,file)])

file_path.sort(key=lambda x: (len(x[0]), x[0]))

for (x, y) in file_path:
    z = basename(Path(y).parent)
    if not os.path.exists(output_path + z):
        os.mkdir(output_path + z)

    f = open(output_path + z + '/' + x, 'a+', buffering=1)
    print("Testing: " + x)
    cmd = str("timeout 3600s " + current_path + "/main" +  " < " +  y)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #Send it to run
    pStartTime = time.time()
    (output, err) = process.communicate()
    pEndTime = time.time()
    #Now store values
    f.write("Testing: " + x + " in " + dt_string + "\n") 
    f.write(str(output.decode('UTF-8')))
    f.flush()
    f.close()
