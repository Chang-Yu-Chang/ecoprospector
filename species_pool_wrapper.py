#!/usr/bin/python

import os
import subprocess

# Read the phenotype list
list_phenotypes = ["f1_additive", "f2_interaction", "f3_additive_binary", "f4_interaction_binary", "f5_invader_growth", "f6_resident_growth"]

#for i in range(3, 21):
#for i in [6, 7, 11, 13, 15, 18, 19, 20]:
#for i in [1, 17]:
#for i in range(21, 41):
for i in range(41, 101):
    for k in range(len(list_phenotypes)):
        # Replace the number at placeholder
        fin = open("template.sh", "rt")
        fout = open("temp.sh", "wt")
        
        for line in fin:
            line = line.replace("arg1", str(i)) 
            line = line.replace("arg2", list_phenotypes[k])
            fout.write(line)


        fin.close(); fout.close()
  
        # Submit a job to cluster
        subprocess.run("sbatch temp.sh", shell=True)
  
        # Remove temp.sh
#        os.remove("temp.sh")
