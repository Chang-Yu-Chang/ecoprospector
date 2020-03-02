#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This python script run the community simulation based on the different combinations of protocols of:
    1. community_function (phenotype)
    2. selection_function
    3. migration_function
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

# Call essential algorithms for community simulation
exec(open("community_simulator-00A-simulation_algorithm.py").read()) # Single algorithm for simulation `simulate_community()`
exec(open("community_simulator-01-regional_pool.py").read())
exec(open("community_simulator-02-community_function.py").read())
exec(open("community_simulator-02A-community_function_student.py").read()) # Student submitted algorithms
exec(open("community_simulator-03-community_selection.py").read())
exec(open("community_simulator-03A-community_selection_student.py").read()) # Student submitted algorithms
exec(open("community_simulator-04-migration.py").read()) 
exec(open("community_simulator-04A-migration_student.py").read()) # Student submitted algorithms


# Call the list of algorithms
list_algorithm = pd.read_csv("data/list_algorithm.csv")
community_phenotype_list = list_algorithm.loc[list_algorithm.Algorithm.isin(['community_phenotype'])]
selection_function_list = list_algorithm.loc[list_algorithm.Algorithm.isin(['selection_function'])]
migration_function_list = list_algorithm.loc[list_algorithm.Algorithm.isin(['migration_function'])]


# Community phenotype of interest in this case is the resource_distance_community_function 
community_phenotype_list = community_phenotype_list.loc[community_phenotype_list.AlgorithmName.isin(['resource_distance_community_function'])]

#gapminder[~gapminder.continent.isin(continents)]

# Migration functions of interest 
#migration_function_list = migration_function_list.loc[community_phenotype_list.AlgorithmName.isin(['migrate_by_fraction', 'migration_into_winner'])]

# Selection functions of interest 
selection_function_list = selection_function_list.loc[selection_function_list.AlgorithmName.isin(['select_best_n'])]


# Simulation
print("Run community simulation")

## Paramaters
rr = 1              # pre-saved index for replicate
n_propagation = 20  # Propogration time
n_wells = 96        # Number of wells
dilution = 1/1000   # Dilution rate
assumptions = a_default.copy() # Assumptions required for community-simulator to work; tournament needs are modified in `simulate_community()`

## 
for index_CommunityPhenotype in range(community_phenotype_list.shape[0]):         # Community phenotype
    for index_SelectionFunction in tqdm(range(selection_function_list.shape[0])): # Selection function
        for index_MigrationFunction in range(migration_function_list.shape[0]):   # Migration function 
            # Pick community phenotype algorithm  
            phenotype_algorithm = community_phenotype_list.iloc[index_CommunityPhenotype]['AlgorithmName']
            
            # Pick selection algorithm
            selection_algorithm = selection_function_list.iloc[index_SelectionFunction]['AlgorithmName']
        
            # Pick migration algrothms
            migration_algorithm = migration_function_list.iloc[index_MigrationFunction]['AlgorithmName']
            
            # Run simulation
            print("\nphenotype algorithm: " + phenotype_algorithm + "\nselection algorithm: " + selection_algorithm + "\nmigration algorithm: " + migration_algorithm)
            try:
                function_time_df = simulate_community(assumptions, phenotype_algorithm, selection_algorithm, migration_algorithm, n = n_wells, n_propagation = n_propagation,  dilution=dilution, replicate_index = rr)
            except:
                print("Something wrong!")
            
            # Convert commnity_function into a df and melt the df
            df = pd.DataFrame(function_time_df)
            df.insert(0, "CommunityPhenotype", np.repeat(phenotype_algorithm, n_propagation))
            df.insert(0, "SelectionFunction", np.repeat(selection_algorithm, n_propagation))
            df.insert(0, "MigrationFuntion", np.repeat(migration_algorithm, n_propagation))
            df.insert(0, "Time", range(n_propagation))
            df = pd.melt(df, id_vars=["Time", "CommunityPhenotype", "SelectionFunction", "MigrationFuntion"], var_name="Well", value_name="PhenotypeValue")
            
        
            # Write the result into csv
            ## File name structure: XX-phenotype_algorithm-XX-selection_algorithm-XX-migration_algorithm-YY.csv
            ## XX means the algorithm ID in the code, where yy means the replicate number 
            df.to_csv("data/" + "{:02d}".format(index_CommunityPhenotype) + "-" + phenotype_algorithm + "-" +"{:02d}".format(index_SelectionFunction) + "-" + selection_algorithm + "-" +"{:02d}".format(index_MigrationFunction) + "-" + migration_algorithm + "-" + "{:02d}".format(rr) + ".csv", index=False)
        


