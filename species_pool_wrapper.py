#!/usr/bin/python

import os
import subprocess

for i in range(5):
  # Replace the number at placeholder
  fin = open("template.sh", "rt")
  fout = open("temp.sh", "wt")
  for line in fin:
    fout.write(line.replace("arg1", str(i)))
  fin.close(); fout.close()
  
  # Submit a job to cluster
  subprocess.os.run("sbatch temp.sh")
  
  # Remove temp.sh
  os.remove("temp.sh")
