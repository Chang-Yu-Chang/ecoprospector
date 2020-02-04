"""
Execution file for bulk simulation
"""


# Community simulator package
from IPython.display import Image
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf as bpdf
import numpy as np
import scipy as sp

# Community selection package
from community_selection import *
from community_selection.A_experiment_functions import *
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_migration_algorithms import *

# Time tracking packages
import datetime
import time
from tqdm import tqdm
print("\nRun community simulation")
print(datetime.datetime.now()) # print time


# Make dynanmics by default we will use the microbial consumer resource model
def dNdt(N,R,params):
    return MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(assumptions)(N,R,params)
dynamics = [dNdt,dRdt]


def select_community(params_algorithm = {"community_phenotype": "f1_community_function_additive", "selection_algorithm": "select_best_n", "migration_algorithm": "no_migration"}):
    """
    Wrapper function for running control(screening) and selection given the selection algorithms and communtiy phenotype
    """
    # Control plate
    np.random.seed(2)
    # Make initial state
    init_state = MakeInitialState(assumptions)
    # Make plate
    plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
    setattr(plate, "species_function", species_function)
    setattr(plate, "interaction_function", interaction_function)
    # Simulation
    ctrl_plate_df, ctrl_function_df = simulate_community(plate, assumptions = assumptions, params_simulation = params_simulation, params_algorithm = {"community_phenotype": params_algorithm["community_phenotype"], "selection_algorithm": "no_selection", "migration_algorithm": params_algorithm["migration_algorithm"]}, write_composition = False)

    # Selection plate
    np.random.seed(2)
    # Make initial state
    init_state = MakeInitialState(assumptions)
    # Make plate
    plate = Community(init_state, dynamics, params, scale = 10**6, parallel = True) 
    setattr(plate, "species_function", species_function)
    setattr(plate, "interaction_function", interaction_function)
    # Simulation
    selection_plate_df, selection_function_df  = simulate_community(plate, assumptions = assumptions, params_simulation = params_simulation, params_algorithm = params_algorithm, write_composition = False)

    return ctrl_function_df, selection_function_df

def concat_result(f_ctrl, f_selc, species_pool_seed, phenotype_algorithm = "f6_resident_growth", selection_algorithm = "select_top25percent", migration_algorithm = "no_migration"):
    """
    Concatenate the result data.frame from control (no-selection) and selection experiment
    """
    # 
    f_ctrl["Experiment"] = "control"
    f_selc["Experiment"] = "selection"
    
    #
    f_result = pd.concat([f_ctrl, f_selc])
    f_result["SelectionAlgorithmName"] = selection_algorithm
    f_result["CommunityPhenotypeName"] = phenotype_algorithm
    f_result["MigrationAlgorithmName"] = phenotype_algorithm
    f_result["SpeciesPool"] = species_pool_seed

    # Write the result
    f_result.to_csv("data/" + selection_algorithm + "-" + phenotype_algorithm + "-" + migration_algorithm + "-pool" + str(species_pool_seed) + ".txt", index = False)


# Simulation




## List of algorithms
list_phenotype = make_algorithm_library().AlgorithmName[make_algorithm_library().AlgorithmType == "community_phenotype"]
#list_phenotype = list_phenotype[list_phenotype == "f3_community_function_additive_binary"]
list_phenotype = list_phenotype[list_phenotype != "f6_resident_growth"] # Remove f6 for now since it has bug
list_phenotype.index = range(len(list_phenotype))

list_selection = make_algorithm_library().AlgorithmName[make_algorithm_library().AlgorithmType == "selection_algorithm"]
list_selection = list_selection[list_selection.isin(["direct_selection_select"])]
#list_selection = list_selection[list_selection.isin(["pair_top", "select_top25percent", "select_top10percent", "pool_top25percent"])]
#list_selection = list_selection[list_selection != "no_selection"] # Remove no-selection from being chosen as selection algorithm.  
list_selection.index = range(len(list_selection))

list_migration = make_algorithm_library().AlgorithmName[make_algorithm_library().AlgorithmType == "migration_algorithm"]
list_migration = list_migration[list_migration.isin(["direct_selection_migrate"])]
list_migration.index = range(len(list_migration))

## Loop for each species pool
for pool_index in range(1):
    # Universal parameters
    assumptions = a_default.copy() # Start with default parameters
    assumptions = {
        'sampling':'Binary', #{'Gaussian','Binary','Gamma'} specifies choice of sampling algorithm    
        'SA': 60*np.ones(3), #Number of species in each specialist family (here, 3 families of 60 species)
        'MA': 30*np.ones(3), #Number of resources in each class 
        'Sgen': 30, #Number of generalist species (unbiased sampling over alll resource classes)
        'muc': 10, #Mean sum of consumption rates (used in all models)
        'sigc': 3, #Standard deviation of sum of consumption rates for Gaussian and Gamma models
        'q': 0.0, #Preference strength of specialist families (0 for generalist and 1 for specialist)
        'c0':0.0, #Sum of background consumption rates in binary model
        'c1':1, #Specific consumption rate in binary model
        'l':0.8, #Leakage fraction
        'fs':0.45, #Fraction of secretion flux with same resource type
        'fw':0.45, #Fraction of secretion flux to 'waste' resource
        'sparsity':0.2, #Effective sparsity of metabolic matrix (between 0 and 1)
        'n_wells':24, #Number of independent wells
        'S':100, #Number of species per well (randomly sampled from the pool of size Stot = sum(SA) + Sgen)
        'food':0, #index of food source (when a single resource is supplied externally)
        'R0_food':1000, #unperturbed fixed point for supplied food
        'regulation':'independent', #metabolic regulation (see dRdt)
        'response':'type III', #functional response (see dRdt)
        'supply':'off' #resource supply (see dRdt)
    }
    
    ## Parameters not included in the community-simulator package
    assumptions.update({
        "m": 0, # Mortality
        "sigma" : 1, # Standard deviation for drawing specifc speices/interaction function
        "alpha": 1 # Scaling factor between species- and interaction-specific function variances
    })
    
    # Prepare experiment setup in this universe
    params, species_pool, species_function, interaction_function = prepare_experiment(assumptions, seed = pool_index + 1)
    
    # Simulation parameters
    params_simulation = {
        "n_propagation": 24, # Length of propagation, or hours within a growth cycle
        "n_transfer": 20, # Number of transfer, or number of passage
        "dilution": 1/125, # Dilution factor at every transfer
        "pool": species_pool 
    }

    print("\nspecies pool seed: " + str(pool_index+1))
    
    # Simulation
    
    for phenotype_index in range(len(list_phenotype)):
        for selection_index in range(len(list_selection)):
            for migration_index in range(len(list_migration)):
                # Algorithms
                phenotype_algorithm = list_phenotype[phenotype_index]
                selection_algorithm = list_selection[selection_index]
                migration_algorithm = list_migration[migration_index]
                print("\nphenotype algorithm: " + phenotype_algorithm + "\nselection algorithm: " + selection_algorithm + "\nmigration algorithm: " + migration_algorithm)
                
                # Output 
                f_ctrl, f_selc = select_community(params_algorithm = {"community_phenotype": phenotype_algorithm, "selection_algorithm": selection_algorithm, "migration_algorithm": migration_algorithm})
                concat_result(f_ctrl, f_selc, pool_index+1, phenotype_algorithm = phenotype_algorithm, selection_algorithm = selection_algorithm)
            

    print(datetime.datetime.now()) # print time
    
print(datetime.datetime.now()) # print time



