#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""

"""
Python functions for simulation in self-assembly, monoculture, and pairwise competition. 
"""

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *

# Import the algorithms
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_migration_algorithms import *


# Functions for experimental setup
def make_regional_pool(assumptions):
    """
    Create a regional species pool that has each species' relative abundance
    
    assumptions = dictionary of metaparameters from community-simulator

    """
    # Total number of species (specialist + generalist)
    S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 

    # Assign drawn values based on power-law distribution
    pool = np.random.power(1, size  = S_tot)
    
    # Relative species abundance in regional pool
    pool = pool/np.sum(pool)
    
    return pool
    

def draw_species_function(assumptions):
    """
    Draw species-specific functions
    
    assumptions = dictionary of metaparameters from community-simulator
    
    Return:
    function_species, function_interaction
    """
    # Total number of species in the pool
    S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 
    
    # Species-specific function, 1-D array
    function_species = np.random.normal(0, assumptions["sigma"], size = S_tot)
    
    # Interaction-specific function, 2-D n by n array
    function_interaction = np.array(np.random.normal(0, assumptions["sigma"] * assumptions["alpha"], size = S_tot * S_tot)).reshape(S_tot, S_tot)

    return function_species, function_interaction


def prepare_experiment(assumptions):
    """
    Prepare the experimental setting shared by all assembly experiments
    
    assumptions = dictionary of metaparameters from community-simulator
    
    Return:
    species_pool
    """
    np.random.seed(0) 
    
    # Make parameters
    params = MakeParams(assumptions) 
    
    # Generate a species pool with the species abundance
    species_pool = make_regional_pool(assumptions) 
    
    # Generate species function
    function_species, function_interaction = draw_species_function(assumptions)
    
    return params, species_pool, function_species, function_interaction


def sample_from_pool(plate_N, pool, scale=10**6, inocula=10**6):
    """
    Sample communities from regional species pool

    plate_N = consumer data.frame
    pool = 1-D array that defines the species relative abundances in the pool
    """
    S_tot = plate_N.shape[0] # Total number of species in the pool
    N0 = np.zeros((plate_N.shape)) # Make empty plate
    consumer_index = plate_N.index
    well_names = plate_N.columns

    # Sample initial community for each well
    # The part is from Jean's code
    for k in range(plate_N.shape[1]):
        pool = np.random.power(1, size  = S_tot) #* PA_vector
        pool = pool/np.sum(pool)
        consumer_list = np.random.choice(len(pool), size=inocula, replace=True, p=pool)
        my_tab = pd.crosstab(index=consumer_list, columns="count") # Calculate the cell count
        N0[my_tab.index.values,k] = np.ravel(my_tab.values / scale) # Scale to biomass
    
    # for k in range(plate_N.shape[1]):
    #     species_list = np.random.choice(len(pool), size=len(pool), replace=True, p=pool)
    #     my_tab = pd.crosstab(index=species_list, columns="count") # Calculate the biomass count
    #     N0[my_tab.index.values,k] = np.ravel(my_tab.values / np.sum(my_tab.values) * inocula / scale) # Scale to sum

    # Make data.frame
    N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

    return N0

    

# Make initial state
## Make monos
def make_synthetic_mono(assumptions):
    """
    Make the synthetic community inocula
    
    assumptions = dictionary of metaparameters
    
    Return:
    N0 = initial consumer populations
    """
    # Extract parameters from assumption
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    F = len(assumptions['SA'])
    
    # Construct lists of names of resources, consumers, resource types, consumer families and wells:
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S_tot)]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    well_names = ['W'+str(k) for k in range(S_tot)]
    
    # Make data.frame for community-simulator input
    N0 = np.eye(S_tot)
    N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)
    
    return N0

## Make synthetic pairs
def make_synthetic_community(species_list, assumptions, number_species = 2, initial_frequency = [[0.5, 0.5], [0.95, 0.05]]):
    """
    Make synthetic community inocula of all pairwise combination
    
    species_list = consumer ID that will be mixed 
    assumptions = dictionary of metaparameters
    number_species = number of species in a community. Set to 2 for pairs, 3 for trios, etc
    initial_frequency = list of pair frequencies

    Return: N0 = initial consumer populations
    """
    # Stopifnot
    #assert max(species_list) <= len(species_pool), "Some species in the list are not in the pool."
    assert len(species_list) >= number_species, "Cannot make pair from one species."
    assert any(list((sum(x) == 1 for x in initial_frequency))), "Sum of initial frequencies is not equal to 1."
    assert any(list((len(x) == number_species for x in initial_frequency))), "Length of initial frequencies is not equal to number of species."
    

    # All possible combinations of species for given number of species added to a well
    from itertools import combinations
    consumer_pairs = list(combinations(species_list, number_species))
    
    # Extract parameters from assumption
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    F = len(assumptions['SA'])
    
    # Construct lists of names of resources, consumers, resource types, consumer families and wells
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S_tot)]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    well_names = ['W'+str(k) for k in range(len(consumer_pairs) * len(initial_frequency))]
    

    # Create empty plate
    N0 = np.zeros(shape = [S_tot, len(consumer_pairs) * len(initial_frequency)])
    
    # Fill the plate with species pairs
    for k in range(len(initial_frequency)):
        for i in range(len(consumer_pairs)):
            N0[consumer_pairs[i], k * len(consumer_pairs) + i] = 1 * initial_frequency[k]

    # Make data.frame for community-simulator input
    N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

    return N0


## Make initial state. Function added from Bobby's code 
def make_initial_state(assumptions, N0):
    """
    Construct initial state, at unperturbed resource fixed point.
    
    assumptions = dictionary of metaparameters
        'SA' = number of species in each family
        'MA' = number of resources of each type
        'Sgen' = number of generalist species
        'n_wells' = number of independent wells in the experiment
        'S' = initial number of species per well
        'food' = index of supplied "food" resource
        'R0_food' = unperturbed fixed point for supplied food resource
    N0 = consumer populations made from make_synthetic_*()
    
    Returns:
    N0 = initial consumer populations
    R0 = initial resource concentrations
    """

    #PREPARE VARIABLES
    #Force number of species to be an array:
    if isinstance(assumptions['MA'],numbers.Number):
        assumptions['MA'] = [assumptions['MA']]
    if isinstance(assumptions['SA'],numbers.Number):
        assumptions['SA'] = [assumptions['SA']]
    #Force numbers of species to be integers:
    assumptions['MA'] = np.asarray(assumptions['MA'],dtype=int)
    assumptions['SA'] = np.asarray(assumptions['SA'],dtype=int)
    assumptions['Sgen'] = int(assumptions['Sgen'])

    #Extract total numbers of resources, consumers, resource types, and consumer families:
    M = int(np.sum(assumptions['MA']))
    T = len(assumptions['MA'])
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    F = len(assumptions['SA'])
    #Construct lists of names of resources, consumers, resource types, consumer families and wells:
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S_tot)]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                      resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    well_names = ['W'+str(k) for k in range(N0.shape[1])] # Modify the number of wells in R0 by using the well number of N0

    R0 = np.zeros((M,N0.shape[1]))

    if not isinstance(assumptions['food'],int):
        assert len(assumptions['food']) == N0.shape[1], 'Length of food vector must equal n_wells.'
        food_list = assumptions['food']
    else:
        food_list = np.ones(N0.shape[1],dtype=int)*assumptions['food']

    if not isinstance(assumptions['R0_food'],int):
        assert len(assumptions['R0_food']) == N0.shape[1], 'Length of food vector must equal n_wells.'
        R0_food_list = assumptions['R0_food']
    else:
        R0_food_list = np.ones(N0.shape[1],dtype=int)*assumptions['R0_food']

    for k in range(N0.shape[1]):
        R0[food_list[k],k] = R0_food_list[k]

    R0 = pd.DataFrame(R0,index=resource_index,columns=well_names)

    return N0, R0

# Simulate community dynamics
## Simulation parameters
params_simulation = {
    "n_propagation": 12, # Lenght of propagation, or hours within a growth cycle
    "n_transfer": 10, # Number of transfer, or number of passage
    "dilution": 1/125, # Dilution factor for transfer
}

## Main function for simulating community
def simulate_community(
    plate, 
    params_simulation, 
    params_algorithm = {"community_phenotype": "community_function_additive", "selection_algorithm": "no_selection", "migration_algorithm": "no_migration"}, 
    file_name = "data/self_assembly-community", 
    assembly_type = "self_assembly",
    write_composition = False):

    """
    Simulate community dynamics by given experimental regimes
    
    plate = plate
    params_simulation = dictionary of parameters for running experiment
    params_algorithm = dictionary of algorithms that determine the selection regime, migration regime, and community pheotypes of interest
    
    Return
    community_composition = concatenated, melted panda dataframe of community and resource composition in each transfer
    community_function = melted panda dataframe of community function
    """
    # Empty list for saving data
    plate_data_list = list()
    community_function_list = list()

    # Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = 0, assembly_type = assembly_type) # Initial state
    plate_data_list.append(plate_data)

    
    # Output the file if write_composition set True
    if write_composition == True:
        plate_data.to_csv(file_name + "-T" + "{:02d}".format(0) + ".txt", index = False)

    # Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        
        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
    
        # Append the composition to a list
        plate_data = reshape_plate_data(plate, transfer_loop_index = i, assembly_type = assembly_type) # Transfer = 0 means that it's before selection regime works upon
        plate_data_list.append(plate_data)

        ## Output the file if write_composition set True
        if write_composition == True:
            plate_data.to_csv(file_name + "-T" + "{:02d}".format(i) + ".txt", index = False) # Transfer = 0 means that it's before selection regime works upon

        # Community phenotype
        community_function = globals()[params_algorithm["community_phenotype"]](plate)
        community_function_list.append(reshape_function_data(community_function, transfer_loop_index = i)) # Transfer = 0 means that it's before selection regime works upon

        # Passage and tranfer matrix
        transfer_matrix = globals()[params_algorithm["selection_algorithm"]](community_function)
        plate.Passage(transfer_matrix * params_simulation["dilution"])
        
        # Migration
        m = globals()[params_algorithm["migration_algorithm"]](community_function) 
        plate.N = migrate_from_pool(plate, params_simulation["pool"], m)
        
        # Print the propagation progress
        print("Transfer " + str(i+1) + " done") 
        
    
    # Concatenate datafrom from different transfers
    plate_data_con = pd.concat(plate_data_list)
    community_function_con = pd.concat(community_function_list)
    
    return plate_data_con, community_function_con


## Reshape the plate resource and consumer matrix for saving into a txt file
def reshape_plate_data(plate, transfer_loop_index, assembly_type):
    # Temporary function for adding variables to and melting df
    def melt_df(plate_df, data_type = "consumer"):
        # Consumers
        temp_df = pd.DataFrame(plate_df)
        total_number = temp_df.shape[0]
        
        ## Add variables
        temp_df["Type"] = np.repeat(data_type, total_number)
        temp_df["ID"] = range(total_number)
        temp_df["Transfer"] = np.repeat(str(transfer_loop_index), total_number)
        temp_df["Assembly"] = np.repeat(assembly_type, total_number)
        
        ## Melt the df
        temp_df = pd.melt(temp_df, id_vars = ["Transfer", "Assembly", "Type", "ID"], var_name = "Well", value_name = "Abundance")
        temp_df = temp_df[["Assembly", "Well", "Transfer", "Type", "ID", "Abundance"]]
        temp_df = temp_df[temp_df.Abundance != 0] # Remove zero abundances
        return temp_df
        
    # Melt the df
    temp_plate = plate.copy() # Copy the original plate 
    df_N = melt_df(temp_plate.N, data_type = "consumer")
    df_R = melt_df(temp_plate.R, data_type = "resource")
    
    # Concatenate dataframes
    merged_df = pd.concat([df_N, df_R]) 
    merged_df["Index"] = list(range(0, merged_df.shape[0]))
    merged_df.set_index("Index", inplace = True)

    return merged_df # Return concatenated dataframe

def reshape_function_data(community_function, transfer_loop_index):
    temp_vector = community_function.copy()
    # Number of wells
    number_well = len(community_function)
    # Make data.frame
    temp_df = pd.DataFrame({"Transfer": np.repeat(str(transfer_loop_index), number_well), "Well": ["W" + str(i) for i in range(number_well)], "CommunityPhenotype": temp_vector})
    # Turn the transfer columns as numeric
    temp_df[["Transfer"]] = temp_df[["Transfer"]].apply(pd.to_numeric)
    
    return temp_df 
 

# Make library of algorithms
def make_algorithm_library():
    import re
    import pandas as pd
    
    # Find directory of community_selection modultes
    import community_selection
    module_dir = community_selection.__file__
    module_dir = re.sub("__init__.py", "", module_dir) 
    
    # 
    algorithm_types = ["community_phenotypes", "selection_algorithms", "migration_algorithms"]
    algorithms = list()
    
    for i in range(len(algorithm_types)):
    
        # Open files
        file_algorithm_phenotype = open(module_dir + ["B", "C", "D"][i] + "_" + algorithm_types[i] + ".py", "r")
        
        # Read lines
        line_list = list()
        line = file_algorithm_phenotype.readline()
        cnt = 1
        
        while line:
            line = file_algorithm_phenotype.readline()
            # Only count non-commented-out functions
            if line.startswith("#") == False: 
                line_list.append(line.strip())
            cnt += 1
        
        # Regular expression
        algorithm_names = re.findall("def \w+", " ".join(line_list))
        
        list_algorithm = [re.sub("^def ", "", x) for x in algorithm_names]
        
        # Write the files
        algorithms.append(pd.DataFrame({"AlgorithmType": re.sub("s$", "", algorithm_types[i]), "AlgorithmName": list_algorithm}))
     
    return pd.concat(algorithms)


## Migrate from species pool to the plate 
def migrate_from_pool(plate, pool, migration_factor = 1):
    # Migration plate
    migration_plate = sample_from_pool(plate.N, pool) * migration_factor
    
    # Migration
    plate_migrated = plate.N + migration_plate 

    return plate_migrated

## Plot community_function
def plot_community_function(function_df):
    function_df.plot.scatter(x = "Transfer", y = "CommunityPhenotype")



















