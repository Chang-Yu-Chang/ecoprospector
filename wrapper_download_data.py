#!/usr/bin/python

import os
import subprocess

for i in range(4, 20):
    try:
        os.mkdir("/Users/chang-yu/Desktop/Lab/community-selection/data/species_pool_" + str(i) + "/")
    except FileExistsError:
        print("Folder already exists")

    # Replace the number at placeholder
    fin = open("template_download.sh", "rt")
    fout = open("temp_download.sh", "wt")

    for line in fin:
        line = line.replace("arg1", str(i)) 
        fout.write(line)
    
    fin.close(); fout.close()
  
        # Run the code
    subprocess.run("source temp_download.sh", shell=True)
  
        # Remove temp.sh
#        os.remove("temp.sh")
