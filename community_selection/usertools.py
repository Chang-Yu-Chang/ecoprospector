#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 09 2020
@author: changyuchang
"""

import numpy as np
import scipy as sp

from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
from community_selection.A_experiment_functions import *
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_migration_algorithms import *


#DEFAULT PARAMETERS FOR CONSUMER AND METABOLIC MATRICES, AND INITIAL STATE
a_default = {'sampling':'Binary', #{'Gaussian','Binary','Gamma'} specifies choice of sampling algorithm
          'SA': 60*np.ones(3), #Number of species in each specialist family (here, 3 families of 60 species)
          'MA': 30*np.ones(3), #Number of resources in each class 
          'Sgen': 30, #Number of generalist species (unbiased sampling over alll resource classes)
          'muc': 10, #Mean sum of consumption rates (used in all models)
          'sigc': 3, #Standard deviation of sum of consumption rates for Gaussian and Gamma models
          'q': 0.0, #Preference strength of specialist families (0 for generalist and 1 for specialist)
          'c0':0.0, #Sum of background consumption rates in binary model
          'c1':1., #Specific consumption rate in binary model
          'l':0.8, #Leakage fraction
          'fs':0.45, #Fraction of secretion flux with same resource type
          'fw':0.45, #Fraction of secretion flux to 'waste' resource
          'sparsity':0.2, #Effective sparsity of metabolic matrix (between 0 and 1)
          'n_wells':10, #Number of independent wells
          'S':100, #Number of species per well (randomly sampled from the pool of size Stot = sum(SA) + Sgen)
          'food':0, #index of food source (when a single resource is supplied externally)
          'R0_food':1000, #unperturbed fixed point for supplied food
          'regulation':'independent', #metabolic regulation (see dRdt)
          'response':'type I', #functional response (see dRdt)
          'supply':'off' #resource supply (see dRdt)
         }

assumptions = a_default.copy() # Start with default parameters

## Update parameters for community-selection
assumptions.update({
    'SA': 600*np.ones(3), #Number of species in each specialist family (here, 3 families of 60 species)
    'MA': 30*np.ones(3), #Number of resources in each class
    'Sgen': 300, #Number of generalist species (unbiased sampling over alll resource classes)
    "n_wells": 10,
    "m": 0, # Mortality
    "scale": 10**6,  #scale is a conversion factor specifying the number of individual microbial cells present when N = 1.
    "sigma" : 1, # Standard deviation for drawing specifc speices/interaction function
    "alpha": 1, # Scaling factor between species- and interaction-specific function variances
    "l": 0, # Set leakage function to 0 to switch off cross-feeding
    "response": "type III",
    "sigma_max": 5,
    'R0_food': 1000, # Total amount of supplied food
    "rich_medium": True, # Number of food types passed to R0
    "binary_threshold": 1,  
    # The parameters below will be passed to params_simulation
    "n_propagation": 1, # Length of propagation, or hours within a growth cycle
    "n_transfer": 10, # Number of total transfer, or number of passage
    "n_transfer_selection": 5, # Number of transfer implementing seleciton regimes
    "dilution": 1/1000, # Dilution factor at every transfer
    "n_inoc": 10**6,  # Number of cells sampled from the regional species at start
    "n_migration": 10**6, # Number of cells to be migrated in the migration perturbation algorithm
    "R_percent": 0, # Fracion of new resources to be spiked in to the media in the resource perturbation algorithm
    "selected_function": "f1_additive"
})

## Prepare experiment setup in this universe
#params, params_simulation = prepare_experiment(assumptions, seed = 1)

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
    function_interaction = np.random.normal(0, assumptions["sigma"] * assumptions["alpha"], size = S_tot * S_tot).reshape(S_tot, S_tot)
    function_interaction_p25 = np.random.binomial(1, 0.25, S_tot**2).reshape(S_tot, S_tot) * np.array(np.random.normal(0, assumptions["sigma"] * assumptions["alpha"], size = S_tot * S_tot)).reshape(S_tot, S_tot)

    return function_species, function_interaction, function_interaction_p25

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
    #species_pool = make_regional_pool(assumptions) 
    
    # Generate species function
    function_species, function_interaction, function_interaction_p25 = draw_species_function(assumptions)
    
    # Subset parameters for simulation
    temp_p = ["n_propagation", "n_transfer", "n_transfer_selection", "dilution", "n_inoc", "selected_function", "n_migration", "R_percent"]
    params_simulation = {key: value for key, value in assumptions.items() if key in temp_p}
    params_simulation.update({
        #"pool": species_pool, 
        "species_function": function_species,
        "interaction_function": function_interaction,
        "interaction_function_p25": function_interaction_p25
    })
    
    return params, params_simulation


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
    
    positive_control = pd.DataFrame({
        "algorithm_name": "positive_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": "no_selection",
        "migration_algorithm": "no_migration"
    })
    
    ## Monoculture
    monoculture = pd.DataFrame({
        "algorithm_name": "monoculture",
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

    ## Pool top 25%
    pool_top25 = pd.DataFrame({
        "algorithm_name": "pool_top25",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top25percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  

    ## Pool top 10%
    pool_top10 = pd.DataFrame({
        "algorithm_name": "pool_top10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top10percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
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
        "selection_algorithm": ["pool_top28percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Penn2004
    Penn2004 = pd.DataFrame({
        "algorithm_name": "Penn2004",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Raynaud2019a
    Raynaud2019a = pd.DataFrame({
        "algorithm_name": "Raynaud2019a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Raynaud2019b
    Raynaud2019b = pd.DataFrame({
        "algorithm_name": "Raynaud2019b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
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
    
    # Swenson2000b_control
    Swenson2000b_control = pd.DataFrame({
        "algorithm_name": "Swenson2000b_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000c
    Swenson2000c = pd.DataFrame({
        "algorithm_name": "Swenson2000c",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
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

    # Xie2019a
    Xie2019a = pd.DataFrame({
        "algorithm_name": "Xie2019a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top_dog" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Xie2019b
    Xie2019b = pd.DataFrame({
        "algorithm_name": "Xie2019b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Arora2019
    Arora2019_V2 = pd.DataFrame({
        "algorithm_name": "Arora2019_V2",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Arora2019" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })    
    Arora2019_V2_control = pd.DataFrame({
        "algorithm_name": "Arora2019_V2_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Arora2019_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })       
    # Raynaud2019a
    Raynaud2019a_V2	= pd.DataFrame({
        "algorithm_name": "Raynaud2019a_V2",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    Raynaud2019a_V2_control	= pd.DataFrame({
        "algorithm_name": "Raynaud2019a_V2_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019a_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    }) 
    # Raynaud2019b
    Raynaud2019b_V2 = pd.DataFrame({
        "algorithm_name": "Raynaud2019b_V2",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019b" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    Raynaud2019b_V2_control = pd.DataFrame({
        "algorithm_name": "Raynaud2019b_V2_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019b_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
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
    
    #Knock_in perutbation using high performing isolates
    knock_in_isolates = pd.DataFrame({
        "algorithm_name": "knock_in_isolates",
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
    #Resource pertubation algorithms
    resource = pd.DataFrame({
        "algorithm_name": "resource",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    resource_old = pd.DataFrame({
        "algorithm_name": "resource_old",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    resource_add = pd.DataFrame({
        "algorithm_name": "resource_add",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    resource_remove = pd.DataFrame({
        "algorithm_name": "resource_remove",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    resource_rescale_add = pd.DataFrame({
        "algorithm_name": "resource_rescale_add",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    resource_rescale_remove = pd.DataFrame({
        "algorithm_name": "resource_rescale_remove",
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
        # Control
        simple_screening, positive_control, monoculture, select_top25, select_top10, pool_top25, pool_top10,
        # Directed selection
        pair_top_communities, multiple_pair_top, directed_selection_migration,
        # Literature
        Arora2019, Blouin2015, Blouin2015_control, Jochum2019, Mueller2019, Panke_Buisse2015, Penn2004,
        Raynaud2019a, Raynaud2019b, Swenson2000a, Swenson2000a_control, Swenson2000b, Swenson2000b_control, Swenson2000c,
        Williams2007a, Williams2007b, Wright2019, Xie2019a, Xie2019b,
        Arora2019_V2, Arora2019_V2_control, Raynaud2019a_V2, Raynaud2019a_V2_control, Raynaud2019b_V2, Raynaud2019b_V2_control,
        # Perturbation
        ctrl, coalescence, migration, bottleneck, knock_out, knock_in, knock_in_isolates,
        resource_old, resource, resource_add, resource_remove, resource_rescale_add, resource_rescale_remove,
        # Iterative perturabation
        iterative_ctrl,iterative_resource,iterative_migration,iterative_resource_migration,iterative_coalescence])
    
    return algorithms
