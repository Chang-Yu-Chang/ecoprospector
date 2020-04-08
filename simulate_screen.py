#!/usr/bin/python

# Read the arguments from command line 
import sys
import os
import pickle

# External arguments
seed_temp = int(sys.argv[1]) # Species pool number
tested_phenotype = str(sys.argv[2])

print("Species pool: " +  str(seed_temp))
print(tested_phenotype)

# Raw data directory
data_directory = "data/test/"

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
    "n_transfer": 20, # Number of total transfer, or number of passage
    "n_transfer_selection": 10, # Number of transfer implementing seleciton regimes
    "dilution": 1/1000, # Dilution factor at every transfer
    "n_inoc": 10**6,  #Number of cells sampled from the regional species at start
    "n_migration": 10**6, # Number of cells to be migrated in the migration perturbation algorithm
    "R_percent": 0, # Fracion of new resources to be spiked in to the media in the resource perturbation algorithm
    "selected_function": "f1_additive"
})

# Prepare experiment setup in this universe
params, params_simulation = prepare_experiment(assumptions, seed = seed_temp)
make_algorithms(params_simulation)

# Lists
list_phenotypes = [tested_phenotype]
list_algorithms = ["simple_screening"]

# Parameters
i = 0; j = 0
params_simulation.update({"selected_function": list_phenotypes[j]}) # selected function
algorithms = make_algorithms(params_simulation) # Make algorithm list 
params_algorithm = algorithms[algorithms["algorithm_name"] == list_algorithms[i]]
file_name = data_directory + "SP" + str(seed_temp) + "-" + list_algorithms[i]
assembly_type = str(list_algorithms[i])

# Simulation
"""
simple screen
"""
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
plate = add_community_function(plate, dynamics, assumptions, params, params_simulation, params_algorithm)

# Save the inocula composition
plate_data = reshape_plate_data(plate, transfer_loop_index = 0, assembly_type = assembly_type, community_function_name = params_algorithm["community_phenotype"][0]) # Initial state
# Save the initial function
community_function = globals()[params_algorithm["community_phenotype"][0]](plate, assumptions = assumptions) # Community phenotype
richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
biomass = list(np.sum(plate.N, axis = 0)) # Biomass
function_data = reshape_function_data(community_function_name = params_algorithm["community_phenotype"][0], community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = 0, assembly_type = assembly_type)        
# Output the plate composition and community functions if write_composition set True
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
    if ((i+1) == assumptions["n_transfer_selection"]) or ((i+1) == assumptions["n_transfer"]):
            plate_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-composition.txt", index = False) # Transfer = 0 means that it's before selection regime works upon
    
    function_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-function.txt", index = False)

    # Save the plate object
    if (i+1) == assumptions["n_transfer_selection"]:
        with open(data_directory + "plate_screen.p", "wb") as f:
            pickle.dump({"plate_N":plate.N, "plate_R":plate.R}, f)

    # Passage and tranfer matrix
    transfer_matrix = globals()[selection_algorithm](community_function)
    plate.Passage(transfer_matrix * params_simulation["dilution"])

    # Migration
    m = globals()[migration_algorithm](community_function) 
    plate.N = migrate_from_pool(plate, migration_factor = m, scale = assumptions["scale"], inocula = params_simulation["n_migration"]) # By default, n_migration is the same as n_inoc

print("\nAlgorithm "+ params_algorithm["algorithm_name"][0] + " finished")


