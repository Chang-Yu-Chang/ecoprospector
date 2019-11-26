#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script tests whether the submitted migration_fnuction has the default format which:
    1. Input an array of community function with length of wells (n)
    2. Ouput a 1-D array with length of n
"""
# Import modules
import time
from tqdm import tqdm
from IPython.display import Image
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf as bpdf
colors = sns.color_palette()

# Call essential functions for community simulator
exec(open("community_simulator-00A-simulation_algorithm.py").read()) # Single algorithm for simulation `simulate_community()`
exec(open("community_simulator-01-regional_pool.py").read())
exec(open("community_simulator-02-community_function.py").read())
exec(open("community_simulator-02A-community_function_student.py").read()) # Student submitted algorithms
exec(open("community_simulator-03-community_selection.py").read())
exec(open("community_simulator-03A-community_selection_student.py").read()) # Student submitted algorithms
exec(open("community_simulator-04-migration.py").read()) 
exec(open("community_simulator-04A-migration_student.py").read()) # Student submitted algorithms

# Call the datafame of student migration algorithms
list_algorithm = pd.read_csv("data/list_algorithm.csv")
migration_function_list = list_algorithm.loc[list_algorithm.Algorithm.isin(['migration_function'])]

# Test migration algorithms in one propagation
print("Test migration algorithms")
test_function_report =  ["" for x in range(migration_function_list.shape[0])] # Empty string array
for j in tqdm(range(migration_function_list.shape[0])):
    migration_algorithm = migration_function_list.iloc[j]['AlgorithmName']
    assumptions = a_default.copy()
    start = 
    
    try:
        test_result = test_migration_function(migration_algorithm, assumptions) # Test whether the selection function passes
        print(migration_algorithm + " passed!")
        test_function_report[j] = migration_algorithm + " passed!"
    except:
        print(migration_algorithm + " does not pass. Something wrong!")
        test_function_report[j] = migration_algorithm + " did not pass. Something wrong!"




## Write the test result into csv
df = pd.DataFrame({"TestResult": test_function_report})
df.to_csv("data/test_migration_function_report.csv", index=False)










