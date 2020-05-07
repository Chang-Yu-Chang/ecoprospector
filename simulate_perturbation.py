#!/usr/bin/python

# Read the arguments from command line 
import sys
import os
import pickle

# External arguments
seed_temp = int(sys.argv[1]) # Species pool number
tested_phenotype = str(sys.argv[2]) # Tested communtiy function. Default to f1_additive 
target_well = str(sys.argv[3]) # Target well for perturbation. Default is the number of rank of community
perturbation = int(sys.argv[4]) # Perturbation treatment
migration_strength = int(sys.argv[5]) # Migration strength n_migration. 1 = knock-in one species, 
resource_strength = float(sys.argv[6]) # Resource stregnth R_percent.

print("\nSpecies pool: " +  str(seed_temp))
print("Tested function: " + tested_phenotype)
print("Target well: " + target_well)
print("Perturbation treatment: " + str(perturbation))
print("Migration strength: " + str(migration_strength))
print("Resource strength: " + str(resource_strength))

# Raw data directory
data_directory = "data/test/"

# Read plate_screen layout
plate_screen = pickle.load(open(data_directory + "plate_screen.p", "rb"))

# Community simulator package
#from IPython.display import Image
from community_simulator import *
from community_simulator.usertools import *
from community_simulator.visualization import *
#import seaborn as sns
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
from community_selection.usertools import *


# Make dynanmics by default we will use the microbial consumer resource model
def dNdt(N,R,params):
    return MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(assumptions)(N,R,params)
dynamics = [dNdt,dRdt]

# Global parameters
## Default parameters from community-simulator
## !!!Don't touch this dictionary!!!
assumptions = a_default.copy() # Start with default parameters

## Update parameters for community-selection
assumptions.update({
    'SA': 2100 * np.ones(1), #Number of species in each specialist family (here, 3 families of 60 species)
    'MA': 90*np.ones(1), #Number of resources in each class
    'Sgen': 0, #Number of generalist species (unbiased sampling over alll resource classes)
    "n_wells": 12,
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
    "n_inoc": 10**6,  #Number of cells sampled from the regional species at start
    "n_migration": 1000, # Number of cells to be migrated in the migration perturbation algorithm
    "s_migration": migration_strength, # Number of species to be migrated in the migration perturbation algorithm
    "R_percent": resource_strength, # Fracion of new resources to be spiked in to the media in the resource perturbation algorithm
    "selected_function": "f1_additive"
})

# Prepare experiment setup in this universe
params, params_simulation = prepare_experiment(assumptions, seed = seed_temp)
make_algorithms(params_simulation)

# Lists
list_phenotypes = [tested_phenotype]
list_algorithms = ["simple_screening"]

if resource_strength != 0:
    list_algorithms = ["resource_add"]
# elif migration_strength == 1:
#     list_algorithms = ["knock_in_isolates"]
elif migration_strength >= 1:
    list_algorithms = ["migration"]

# Parameters
i = 0; j = 0
params_simulation.update({"selected_function": list_phenotypes[j]}) # selected function
algorithms = make_algorithms(params_simulation) # Make algorithm list 
params_algorithm = algorithms[algorithms["algorithm_name"] == list_algorithms[i]] # Algorithm for this run
file_name = data_directory + "SP" + str(seed_temp) + "-" + target_well + "-P" + str(perturbation)
assembly_type = str(list_algorithms[i])

print("\nAlgorithm: "+ params_algorithm["algorithm_name"][0])
print("\n")
print(params_algorithm[["transfer", "community_phenotype", "selection_algorithm", "migration_algorithm"]].to_string(index = False))





# Simulation 
"""
Perturbation
"""
# Set seeds
np.random.seed(2)

# Make initial state
init_state = MakeInitialState(assumptions)

# Make plate
plate = Community(init_state, dynamics, params, scale = assumptions["scale"], parallel = True) 

# Update the community composition 
print("\nGenerating initial plate")
plate.N = plate_screen["plate_N"]

# Update the supplied resource if assumptions["rich_medium"]
plate.R = plate_screen["plate_R"]
if assumptions["rich_medium"]:
    plate.R0 = make_rich_medium(plate.R, assumptions) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on

# Add the attributes that are essential to the function measurement to the plate objects 
print("\nAdding attributes that are essential to the community function to the plate object")
plate = add_community_function(plate, dynamics, assumptions, params, params_simulation, params_algorithm)


# Algorithms used in this transfer
phenotype_algorithm = params_algorithm["community_phenotype"][params_simulation["n_transfer_selection"]-1]
selection_algorithm = params_algorithm["selection_algorithm"][params_simulation["n_transfer_selection"]-1]
migration_algorithm = params_algorithm["migration_algorithm"][params_simulation["n_transfer_selection"]-1]


# Refill the plate with media
community_function = globals()[params_algorithm["community_phenotype"][0]](plate, assumptions = assumptions) # Community phenotype
transfer_matrix = select_top_nth(community_function, int(target_well[2:len(target_well)])) # Select the top nth community as the targeted perturbed
plate.Passage(transfer_matrix * params_simulation["dilution"])

# Perturbation
## Migration
m = globals()[migration_algorithm](community_function)
plate.N = migrate_from_pool(plate, migration_factor = m, assumptions = assumptions, power_law = False, community_function = community_function) # By default, n_migration is the same as n_inoc

## Resource 
if 'resource' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
    winning_index = np.where(community_function >= np.max(community_function))[0][0]
    #Remove fresh environment that was added by passage
    plate.R = plate.R - plate.R0
    #change default fresh renvironment so that all subsequent rounds use R0
    for k in plate.R0.columns:
        if k != plate.R0.columns[winning_index]:
			#By default pick 2 resources at randomk
            r_id_remove = np.random.choice(np.where(plate.R0[k]>=0)[0])
            r_id_add = np.random.choice(np.where(plate.R0[k]>=0)[0])
            if 'add' in params_algorithm["algorithm_name"][0]: #remove from top and add to random
                r_id_remove = np.where(plate.R0[k]==np.max(plate.R0[k]))[0]
            if 'remove' in params_algorithm["algorithm_name"][0]: #remove from random and add to bottom
                r_id_add = np.where(plate.R0[k]==np.min(plate.R0[k]))[0]
            if 'rescale_add' in params_algorithm["algorithm_name"][0]:  # Increase Fraction of resource
                plate.R0[k][r_id_add] = plate.R0[k][r_id_add]*(1+params_simulation['R_percent']) #increase resource conc by fixed %
            elif 'rescale_remove' in params_algorithm["algorithm_name"][0]: # Decrease Fraction of resource
                plate.R0[k][r_id_remove] = plate.R0[k][r_id_remove]*(1-params_simulation['R_percent'])  #decrease resource
            elif 'resource_old' in params_algorithm["algorithm_name"][0]:
                plate.R0[k] = plate.R0[k] * (1-params_simulation['R_percent']) #Dilute old resources
                plate.R0[k][r_id_add] = plate.R0[k][r_id_add] + (assumptions['R0_food']*params_simulation['R_percent'])
            else: #Default to resource swap.
                plate.R0[k][r_id_add] = plate.R0[k][r_id_add] + (plate.R0[k][r_id_remove]*params_simulation['R_percent']) #add new resources
                plate.R0[k][r_id_remove] = plate.R0[k][r_id_remove]*(1-params_simulation['R_percent']) #remove new resources
    plate.R0 = plate.R0/np.sum(plate.R0)*assumptions['R0_food'] #Keep this to avoid floating point error and rescale when neeeded.
    #add new fresh environment (so that this round uses R0
    plate.R = plate.R + plate.R0


# Propagation
print("\nStart propogation")
# Run simulation. Starting from the first transfer of stablization
for i in range(params_simulation["n_transfer_selection"], params_simulation["n_transfer"]):
    # Algorithms used in this transfer
    phenotype_algorithm = params_algorithm["community_phenotype"][i]
    selection_algorithm = params_algorithm["selection_algorithm"][i]
    migration_algorithm = params_algorithm["migration_algorithm"][i]

    # Print the propagation progress
    print("Transfer " + str(i+1))

    # Propagation
    plate.Propagate(params_simulation["n_propagation"])

    # Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1, assembly_type = assembly_type, community_function_name = params_algorithm["community_phenotype"][0]) # Initial state
    plate_data["Perturbation"] = perturbation; plate_data["MigrationStrength"] = migration_strength; plate_data["ResourceStrength"] = resource_strength
    # Save the initial function
    community_function = globals()[params_algorithm["community_phenotype"][0]](plate, assumptions = assumptions) # Community phenotype
    richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
    biomass = list(np.sum(plate.N, axis = 0)) # Biomass
    function_data = reshape_function_data(community_function_name = params_algorithm["community_phenotype"][0], community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = i+1, assembly_type = assembly_type)        
    function_data["Perturbation"] = perturbation; function_data["MigrationStrength"] = migration_strength; function_data["ResourceStrength"] = resource_strength
    # Output the plate composition and community functions if write_composition set True
    if (i == assumptions["n_transfer_selection"]) or (i+1) == assumptions["n_transfer"]:
            plate_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-composition.txt", index = False) # Transfer = 0 means that it's before selection regime works upon
    function_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-function.txt", index = False)

    # Passage and tranfer matrix
    transfer_matrix = globals()["no_selection"](community_function)
    plate.Passage(transfer_matrix * params_simulation["dilution"])

    # Migration
    m = globals()["no_migration"](community_function) 
    plate.N = migrate_from_pool(plate, migration_factor = m, assumptions = assumptions, community_function = community_function) # By default, n_migration is the same as n_inoc

print("\nAlgorithm "+ params_algorithm["algorithm_name"][0] + " finished")


