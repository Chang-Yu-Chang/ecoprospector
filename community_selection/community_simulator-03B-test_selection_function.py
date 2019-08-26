#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script tests whether the submitted selection_function has the default format which:
    1. Input an array of community function with length of wells (n)
    2. Ouput a transfer matrix with dimension n by n
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

# Call the datafame of student selection functions
list_algorithm = pd.read_csv("data/list_algorithm.csv")
selection_function_list = list_algorithm.loc[list_algorithm.Algorithm.isin(['selection_function'])]

# Test selection functions in one propagation
print("Test selection functions")
test_function_report =  ["" for x in range(selection_function_list.shape[0])] # Empty string array
for j in tqdm(range(selection_function_list.shape[0])):
    selection_algorithm = selection_function_list.iloc[j]['AlgorithmName']
    assumptions = a_default.copy()
    
    try:
        test_result = test_selection_function(selection_algorithm, assumptions) # Test whether the selection function passes
        print(selection_algorithm + " passed!")
        test_function_report[j] = selection_algorithm + " passed!"
    except:
        print(selection_algorithm + " does not pass. Something wrong!")
        test_function_report[j] = selection_algorithm + " did not pass. Something wrong!"

## Write the test result into csv
df = pd.DataFrame({"TestResult": test_function_report})
df.to_csv("data/test_selection_function_report.csv", index=False)










