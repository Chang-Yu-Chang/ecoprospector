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
import itertools

from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *

# Import the algorithms
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_migration_algorithms import *


# Functions for experiment

def make_regional_pool(assumptions):
    """
    Create a regional species pool that has each species' relative abundance    
    
    assumptions = dictionary of metaparameters from community-simulator
    """
    # Total number of species (specialist + generalist)
    S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 
    
    # Assign drawn values based on power-law distribution
    pool = np.random.power(1, size  = S_tot)  

    # Relative speceis abundance in regional pool
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


def prepare_experiment(assumptions, seed = 1):
    """
    Prepare the experimental setting shared by all assembly experiments
    
    assumptions = dictionary of metaparameters from community-simulator
    
    Return:
    species_pool
    """
    np.random.seed(seed) 
    
    # Make parameters
    params = MakeParams(assumptions) 
    
    # Generate a species pool with the species abundance
    species_pool = make_regional_pool(assumptions) 
    
    # Generate species function
    function_species, function_interaction = draw_species_function(assumptions)
    
    # Subset parameters for simulation
    temp_p = ["n_propagation", "n_transfer", "n_transfer_selection", "dilution", "n_inoc", "selected_function"]
    params_simulation = {key: value for key, value in assumptions.items() if key in temp_p}
    params_simulation.update({
        "pool": species_pool, 
        "species_function": function_species,
        "interaction_function": function_interaction
    })
    
    return params, params_simulation

def sample_from_pool(plate_N, scale = 10**6, inocula = 10**6,migration=False):
    """
    Sample communities from regional species pool.
    In order to create variability in the pool, split the species pool into two pools, one for initial inocula and one migration.
    For initial inocula, use the even number species, whereas for initla pool, use odds number species.

    plate_N = consumer data.frame
    pool = 1-D array that defines the species relative abundances in the pool
    """
    S_tot = plate_N.shape[0] # Total number of species in the pool
    N0 = np.zeros((plate_N.shape)) # Make empty plate
    consumer_index = plate_N.index
    well_names = plate_N.columns

    # Sample initial community for each well
    for k in range(plate_N.shape[1]):
        # For each well, sample community from different microbiome sample
        if migration==False:
            np.random.seed(k + 1) 
        pool = np.random.power(0.01, size = S_tot) # Power-law distribution
        pool = pool/np.sum(pool) # Normalize the pool
        consumer_list = np.random.choice(S_tot, size = inocula, replace = True, p = pool) # Draw from the pool
        my_tab = pd.crosstab(index = consumer_list, columns = "count") # Calculate the cell count
        N0[my_tab.index.values,k] = np.ravel(my_tab.values / scale) # Scale to biomass
        
    # Make data.frame
    N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

    return N0

# Migrate from species pool to the plate 
def migrate_from_pool(plate, pool, migration_factor, scale, inocula):
    if np.sum(migration_factor) !=0:
        print('Migration')
        migration_plate = sample_from_pool(plate.N, scale = scale, inocula = inocula, migration = True) * migration_factor # Migration factor is a list determined by migration algorithms and community function
        plate_migrated = plate.N + migration_plate 
        
    else:
        plate_migrated = plate.N
        
    return plate_migrated


# Make rich medium
def make_rich_medium(plate_R, assumptions):
    """
    Make rich media for each well

    plate_R = resource dataframe
    """
    np.random.seed(2)
    
    # Total number of resource in this universe
    R_tot = plate_R.shape[0] 
    
    # Make empty plate
    R0 = np.zeros((plate_R.shape)) # Make empty plate
    
    # Resource index
    resource_index = plate_R.index 
    
    # Well index
    well_names = plate_R.columns
    
    resource_pool = np.random.uniform(0, 1, size = R_tot) # Uniform distribution
    resource_pool = resource_pool/np.sum(resource_pool)
    resource_list = np.random.choice(R_tot, size = assumptions["R0_food"], replace = True, p = resource_pool) # Draw from the pool
    my_tab = pd.crosstab(index = resource_list, columns = "count")
    food_compostion = np.ravel(my_tab.values)
    for i in range(plate_R.shape[1]):
        R0[my_tab.index.values,i] = food_compostion

    R0 = pd.DataFrame(R0, index = resource_index, columns = well_names)
    
    return R0


# Main function for simulating community

def simulate_community( 
    assumptions,
    params,
    dynamics,
    params_simulation,
    params_algorithm,
    file_name = "data/",
    assembly_type = "self_assembly",
    write_composition = False):
    """
    Simulate community dynamics by given experimental regimes
    
    params_simulation = dictionary of parameters for running experiment
    params_algorithm = dictionary of algorithms that determine the selection regime, migration regime, and community pheotypes
    
    Return:
    community_composition = concatenated, melted panda dataframe of community and resource composition in each transfer
    community_function = melted panda dataframe of community function
    """

    # Print out the algorithms
    print("\nAlgorithm: "+ params_algorithm["algorithm_name"][0])
    print("\n")
    print(params_algorithm[["transfer", "community_phenotype", "selection_algorithm", "migration_algorithm"]].to_string(index = False))

    # Set seeds
    np.random.seed(2)
    
    # Make initial state
    init_state = MakeInitialState(assumptions)

    # Make plate
    plate = Community(init_state, dynamics, params, scale = assumptions["scale"], parallel = True) 
    
    # Update the community composition by sampling from the pool
    print("\nGenerating initial plate")
    plate.N = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])
    
    # Update the supplied resource if assumptions["rich_medium"]
    if assumptions["rich_medium"]:
        plate.R = make_rich_medium(plate.R, assumptions)
        plate.R0 = make_rich_medium(plate.R, assumptions) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on

    # Add the attributes that are essential to the function measurement to the plate objects 
    print("\nAdding attributes that are essential to the community function to the plate object")
    plate = add_community_function(plate, dynamics, assumptions, params_simulation)
    
    # Empty list for saving data
    plate_data_list = list() # Plate composition
    community_function_list = list() # Community function

    # Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = 0, assembly_type = assembly_type, community_function_name = params_algorithm["community_phenotype"][0]) # Initial state
    plate_data_list.append(plate_data)
    
    # Save the initial function
    community_function = globals()[params_algorithm["community_phenotype"][0]](plate, assumptions = assumptions) # Community phenotype
    richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
    biomass = list(np.sum(plate.N, axis = 0)) # Biomass
    function_data = reshape_function_data(community_function_name = params_algorithm["community_phenotype"][0], community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = 0, assembly_type = assembly_type)        
    community_function_list.append(function_data) # Transfer = 0 means that it's before selection regime works upon


    # Output the plate composition and community functions if write_composition set True
    if write_composition == True:
        plate_data.to_csv(file_name + "-" + params_algorithm["community_phenotype"][0] + "-T" + "{:02d}".format(0) + "-composition.txt", index = False)
        function_data.to_csv(file_name + "-" + params_algorithm["community_phenotype"][0] + "-T" + "{:02d}".format(0) + "-function.txt", index = False)

    print("\nStart propogation")
    # Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        # Algorithms used in this transfer
        phenotype_algorithm = params_algorithm["community_phenotype"][i]
        selection_algorithm = params_algorithm["selection_algorithm"][i]
        migration_algorithm = params_algorithm["migration_algorithm"][i]

        # Print the propagation progress
        if (i % 5) == 0:
            print("Transfer " + str(i+1))

        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
    
        # Append the composition to a list
        plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1, assembly_type = assembly_type, community_function_name = phenotype_algorithm) # Transfer = 0 means that it's before selection regime works upon
        plate_data_list.append(plate_data)

        # Community phenotype, richness, and biomass
        community_function = globals()[phenotype_algorithm](plate, assumptions = assumptions) # Community phenotype
        richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
        biomass = list(np.sum(plate.N, axis = 0)) # Biomass
        function_data = reshape_function_data(community_function_name = phenotype_algorithm, community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = i + 1 , assembly_type = assembly_type)        
        community_function_list.append(function_data) # Transfer = 0 means that it's before selection regime works upon

        # Output the plate composition and community functions if write_composition set True
        if write_composition == True:
            plate_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-composition.txt", index = False) # Transfer = 0 means that it's before selection regime works upon
            function_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-function.txt", index = False)

        # Passage and tranfer matrix if is selection experiment

        # Passage and tranfer matrix
        transfer_matrix = globals()[selection_algorithm](community_function)
        plate.Passage(transfer_matrix * params_simulation["dilution"])
        
        # Migration
        m = globals()[migration_algorithm](community_function) 
        plate.N = migrate_from_pool(plate, pool = params_simulation["pool"], migration_factor = m, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])

        # Perturbation
        if params_algorithm["algorithm_name"][0] == 'knock_in' and selection_algorithm == 'select_top':
            winning_index = np.where(community_function >= np.max(community_function))[0][0] 
            for k in plate.N.columns:
                if k != plate.N.columns[winning_index]:
                    print(k)
                    s_id = np.random.choice(np.where(plate.N[k]==0)[0])
                    plate.N[k][s_id]= 1/params_simulation["dilution"] * 1/assumptions["scale"]
        if params_algorithm["algorithm_name"][0] == 'knock_out' and selection_algorithm == 'select_top':
            winning_index = np.where(community_function >= np.max(community_function))[0][0]
            for k in plate.N.columns:
                if k != plate.N.columns[winning_index]:
                    s_id = np.random.choice(np.where(plate.N[k]>0)[0])
                    plate.N[k][s_id]=0        
        if params_algorithm["algorithm_name"][0] == 'bottleneck' and selection_algorithm == 'select_top':
            winning_index = np.where(community_function >= np.max(community_function))[0][0] 
            dilution_matrix = np.eye(assumptions['n_wells'])
            dilution_matrix[winning_index,winning_index] = 0
            plate.Passage(np.eye(assumptions['n_wells'])* params_simulation["dilution"])
        if 'resource' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
            winning_index = np.where(community_function >= np.max(community_function))[0][0] 
			#Remove fresh environment that was added by passage
            plate.R = plate.R - plate.R0
			#change default fresh renvironment so that all subsequent rounds use R0
            for k in plate.R0.columns:
                if k != plate.R0.columns[winning_index]: 
                    r_id = np.random.choice(np.where(plate.R0[k]>=0)[0])
                    plate.R0[k][r_id] = assumptions['R0_food']/10
            plate.R0 = plate.R0/np.sum(plate.R0)*assumptions['R0_food']
            ##add new fresh environment (so that this round uses R0
            plate.R = plate.R + plate.R0
            
    print("\nAlgorithm "+ params_algorithm["algorithm_name"][0] + " finished")

    # Concatenate data from from different transfers
    plate_data_con = pd.concat(plate_data_list)
    community_function_con = pd.concat(community_function_list)

    return plate_data_con, community_function_con


## Reshape the plate resource and consumer matrix for saving into a txt file
def reshape_plate_data(plate, transfer_loop_index, assembly_type, community_function_name):
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
        temp_df["CommunityPhenotypeName"] = np.repeat(community_function_name, total_number)

        ## Melt the df
        temp_df = pd.melt(temp_df, id_vars = ["Transfer", "CommunityPhenotypeName", "Assembly", "Type", "ID"], var_name = "Well", value_name = "Abundance")
        temp_df = temp_df[["Assembly", "CommunityPhenotypeName", "Well", "Transfer", "Type", "ID", "Abundance"]]
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


def reshape_function_data(community_function_name, community_function, richness, biomass, transfer_loop_index, assembly_type):
    temp_vector1 = community_function.copy()
    temp_vector2 = richness.copy()
    temp_vector3 = biomass.copy()
    
    # Number of wells
    number_well = len(richness)
    # Make data.frame
    temp_df = pd.DataFrame({"Assembly": np.repeat(assembly_type, number_well),
    "CommunityPhenotypeName": np.repeat(community_function_name, number_well),
    "Well": ["W" + str(i) for i in range(number_well)], 
    "Transfer": np.repeat(str(transfer_loop_index), number_well), 
    "CommunityPhenotype": temp_vector1,
    "Richness": temp_vector2,
    "Biomass": temp_vector3})
    
    # Turn the transfer columns as numeric
    temp_df[["Transfer"]] = temp_df[["Transfer"]].apply(pd.to_numeric)
    
    return temp_df 

 

# Make library of algorithms
def make_algorithm_library():
    """
    Show the table of algorithms in this package
    """
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
            line_list.append(line.strip())
            cnt += 1
        
        # Regular expression
        algorithm_names = re.findall("def \w+", " ".join(line_list))
        list_algorithm = [re.sub("^def ", "", x) for x in algorithm_names]
        
        # Write the files
        algorithms.append(pd.DataFrame({"AlgorithmType": re.sub("s$", "", algorithm_types[i]), "AlgorithmName": list_algorithm}))
     
    return pd.concat(algorithms)


# Plot community function
def plot_community_function(function_df):
    function_df.plot.scatter(x = "Transfer", y = "CommunityPhenotype")

# Plot transfer matrix
def plot_transfer_matrix(transfer_matrix):
    import seaborn as sns
    fig,ax=plt.subplots()
    sns.heatmap(transfer_matrix,ax=ax)
    ax.set_xlabel('Old well',fontsize=14)
    ax.set_ylabel('New well',fontsize=14)
    ax.set_title(r'Transfer Matrix',fontsize=14)
    plt.show()
    
    return fig
    


def make_algorithms(params_simulation):
    # Algorithms
    ## Simple screening
    simple_screening = pd.DataFrame({
        "algorithm_name": "simple_screening",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": "no_selection",
        "migration_algorithm": "no_migration"
    })

    ## Direction selection
    directed_selection_migration = pd.DataFrame({
        "algorithm_name": "directed_selection_migration",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"], 
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["directed_selection_select"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": ["no_migration" for i in range(params_simulation["n_transfer_selection"]-1)] + ["directed_selection_migrate"] + ["no_migration" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])]
    })

    ## Select top 25%
    select_top25 = pd.DataFrame({
        "algorithm_name": "select_top25",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top25percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  

    ## Select top 10%
    select_top10 = pd.DataFrame({
        "algorithm_name": "select_top10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top10percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    ## Pair top communities
    pair_top_communities = pd.DataFrame({
        "algorithm_name": "pair_top_communities",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pair_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    ## Multiple pair top
    multiple_pair_top = pd.DataFrame({
        "algorithm_name": "multiple_pair_top",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pair_top" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Arora2019
    Arora2019 = pd.DataFrame({
        "algorithm_name": "Arora2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top33percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Blouin2015
    Blouin2015 = pd.DataFrame({
        "algorithm_name": "Blouin2015",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Blouin2015 control
    Blouin2015_control = pd.DataFrame({
        "algorithm_name": "Blouin2015_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Jochum2019
    Jochum2019 = pd.DataFrame({
        "algorithm_name": "Jochum2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Mueller2019
    Mueller2019 = pd.DataFrame({
        "algorithm_name": "Mueller2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Panke-Buisse2015
    Panke_Buisse2015 = pd.DataFrame({
        "algorithm_name": "Panke_Buisse2015",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000a
    Swenson2000a = pd.DataFrame({
        "algorithm_name": "Swenson2000a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000a control
    Swenson2000a_control = pd.DataFrame({
        "algorithm_name": "Swenson2000a_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    # Swenson2000b
    Swenson2000b = pd.DataFrame({
        "algorithm_name": "Swenson2000b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Williams2007a
    Williams2007a = pd.DataFrame({
        "algorithm_name": "Williams2007a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Williams2007b
    Williams2007b = pd.DataFrame({
        "algorithm_name": "Williams2007b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007b" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
        
    # Wright2019
    Wright2019 = pd.DataFrame({
        "algorithm_name": "Wright2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    #ctrl_pertubation
    ctrl = pd.DataFrame({
        "algorithm_name": "ctrl",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    
    #Knock_in_pertubation
    knock_in = pd.DataFrame({
        "algorithm_name": "knock_in",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    #Knock_out_pertubation
    knock_out = pd.DataFrame({
        "algorithm_name": "knock_out",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    #Bottleneck_Pertubation
    bottleneck = pd.DataFrame({
        "algorithm_name": "bottleneck",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })      
    
    #Coalescence_Pertubation
    coalescence = pd.DataFrame({
        "algorithm_name": "coalescence",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["coalescence"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })     
    
    #Migration Pertubation
    migration = pd.DataFrame({
        "algorithm_name": "migration",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": ["no_migration" for i in range(params_simulation["n_transfer_selection"]-1)] + ["parent_migration"] + ["no_migration" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])]
    })     
    
    #Resource pertubation
    resource = pd.DataFrame({
        "algorithm_name": "resource",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    #long experiments
    iterative_ctrl = pd.DataFrame({
        "algorithm_name": "iterative_ctrl",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm": "no_migration"
    })  
    
    iterative_resource = pd.DataFrame({
        "algorithm_name": "iterative_resource",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm": "no_migration"
    })  

    iterative_resource_migration = pd.DataFrame({
        "algorithm_name": "iterative_resource_migration",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm":  np.tile(["no_migration" for i in range(params_simulation["n_transfer_selection"]-1)] + ["parent_migration"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist()
    })  

    iterative_migration = pd.DataFrame({
        "algorithm_name": "iterative_migration",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm":  np.tile(["no_migration" for i in range(params_simulation["n_transfer_selection"]-1)] + ["parent_migration"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist()
    })  

    iterative_coalescence = pd.DataFrame({
        "algorithm_name": "iterative_coalescence",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["coalescence"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm": "no_migration"
    })   
    
    # Save the algorithms
    algorithms = pd.concat([
        simple_screening, directed_selection_migration, select_top25, select_top10, pair_top_communities, multiple_pair_top,
        Arora2019, Blouin2015, Blouin2015_control, Mueller2019, Panke_Buisse2015, Swenson2000a, Swenson2000a_control,
        Swenson2000b, Williams2007a, Williams2007b, Wright2019,
        ctrl, coalescence,migration,resource,bottleneck,knock_out,knock_in,
        iterative_ctrl,iterative_resource,iterative_migration,iterative_resource_migration,iterative_coalescence])
    
    return algorithms


def add_community_function(plate, dynamics, assumptions, params_simulation):
    """
    Add the function attribute to the community
    
    For f1 and f3, add species_function 
    For f2 and f4, add interaction_function
    For f5, add invasion_plate_t0 and invasion_plate_t1
    For f6, add XX
    
    """
    np.random.seed(2)
    # Species function, f1 and f3
    setattr(plate, "species_function", params_simulation["species_function"]) # Species function for additive community function

    # Interactive functions, f2 and f4
    setattr(plate, "interaction_function", params_simulation["interaction_function"]) # Interactive function for interactive community function

    # Invasion function f5 and resident function f6
    if (params_simulation["selected_function"] == "f5_invader_growth") or (params_simulation["selected_function"] == "f6_resident_growth"):
    # Make 96 communities and pick the best grown one as the focal invader community
        assumptions_invasion = assumptions.copy()
        
        if params_simulation["selected_function"] == "f5_invader_growth":
            assumptions_invasion.update({"n_wells": 96, "n_inoc": 1})
        elif params_simulation["selected_function"] == "f6_resident_growth":
            assumptions_invasion.update({"n_wells": 96})
        
        params_invasion = MakeParams(assumptions_invasion)
        init_state_invasion = MakeInitialState(assumptions_invasion)
        plate_invasion = Community(init_state_invasion, dynamics, params_invasion, scale = assumptions_invasion["scale"], parallel = True)
        plate_invasion.N = sample_from_pool(plate_invasion.N, scale = assumptions_invasion["scale"], inocula = assumptions_invasion["n_inoc"])
        plate_invasion.N = plate_invasion.N / np.sum(plate_invasion.N, axis = 0) # Rescale biomass to to one 


        if assumptions["rich_medium"]:
            plate_invasion.R = make_rich_medium(plate_invasion.R, assumptions_invasion)
            plate_invasion.R0 = make_rich_medium(plate_invasion.R, assumptions_invasion) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on 

        # Save the t0 plate
        temp = assumptions_invasion["n_transfer"] - assumptions_invasion["n_transfer_selection"]
        print("\nStabilizing the invader (or resident) community. Passage for " + str(temp) + " transfer. The plate has ", str(assumptions_invasion["n_wells"]), " wells.")
        plate_invasion_t0 = plate_invasion.N.copy()
        plate_invasion.Propagate(assumptions_invasion["n_propagation"])
        print("Passaging invader community. Transfer 1")

        # Grow the invader plate 
        for i in range(0, temp - 1):
            plate_invasion.Passage(np.eye(assumptions_invasion["n_wells"]) * assumptions_invasion["dilution"])
            plate_invasion.Propagate(assumptions_invasion["n_propagation"])
            print("Passaging invader community. Transfer " + str(i + 2))


        # Save the t1 plate
        plate_invasion_t1 = plate_invasion.N.copy()
        invasion_plate_growth = np.sum(plate_invasion.N, axis = 0)
        temp_index = np.where(invasion_plate_growth == np.max(invasion_plate_growth))[0][0] # Find the well with the highest biomass
        temp_column_t0 = plate_invasion_t0["W" + str(temp_index)]
        temp_column_t1 = plate_invasion_t1["W" + str(temp_index)] # Pick the best growth community as the invader community
    
        # Dupliate the chosen community to the whole plate
        temp_df_t0 = pd.DataFrame()
        temp_df_t1 = pd.DataFrame()
        for i in range(assumptions["n_wells"]):
            temp_df_t0["W" + str(i)] = temp_column_t0
            temp_df_t1["W" + str(i)] = temp_column_t1

        print("Finished passaging the invader (or resident) community.") 
        print("The community has " + str(np.sum(temp_column_t1>0)) + " species")
        print("The community has initial biomass " + str(int(np.sum(temp_column_t0))) + " and reaches total biomass", str(np.sum(temp_column_t1)))
        
    
        if params_simulation["selected_function"] == "f5_invader_growth":
            # Add the invasion plate to the attr of community
            setattr(plate, "invasion_plate_t0", temp_df_t0)
            setattr(plate, "invasion_plate_t1", temp_df_t1)
        elif params_simulation["selected_function"] == "f6_resident_growth":
            # Add the resident plate to the attr of community
            setattr(plate, "resident_plate_t0", temp_df_t0)
            setattr(plate, "resident_plate_t1", temp_df_t1)
    
    return plate



