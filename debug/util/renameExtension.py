import os
from pathlib import Path

current_path = str(Path(__file__).parent.absolute())

input_dir = current_path + "/input/"

for path, subdirs, files in os.walk(input_dir):
        for name in files:
            filePath = os.path.join(path, name)
            target_filePath = ''.join(filePath.split('.')[:-1])+""
            os.rename(filePath, target_filePath)
