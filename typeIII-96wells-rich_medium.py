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
    'SA': 600*np.ones(3), #Number of species in each specialist family (here, 3 families of 60 species)
    'MA': 30*np.ones(3), #Number of resources in each class 
    'Sgen': 300, #Number of generalist species (unbiased sampling over alll resource classes)
    "n_wells": 96,
    "c1": 1, #Rescale uptake rate part 1. This is needed to avoid numerical errors that slow down the simulations
    "muc": 10, # Rescale uptake part 2
    "m": 0, # Mortality
    "scale": 10**6,  #scale is a conversion factor specifying the number of individual microbial cells present when N = 1.
    "sigma" : 1, # Standard deviation for drawing specifc speices/interaction function
    "alpha": 1, # Scaling factor between species- and interaction-specific function variances
    "l": 0, # Set leakage function to 0 to switch off cross-feeding
    "response": "type III",
    "sigma_max": 100, 
    'R0_food': 1000, # Total amount of supplied food 
    "rich_medium": True, # Number of food types passed to R0
    # The parameters below will be passed to params_simulation
    "n_propagation": 8, # Length of propagation, or hours within a growth cycle
    "n_transfer": 40, # Number of total transfer, or number of passage
    "n_transfer_selection": 20, # Number of transfer implementing seleciton regimes 
    "dilution": 1/125, # Dilution factor at every transfer
    "n_inoc": 128,  #Number of cells sampled from the regional species at start
    "selected_function": "f1_additive"
})

# Prepare experiment setup in this universe
params, params_simulation = prepare_experiment(assumptions, seed = 1)

data_directory = "data/typeIII-96wells-rich_medium/"
list_phenotypes = ["f1_additive", "f2_interaction", "f3_additive_binary", "f4_interaction_binary", "f5_invader_growth", "f6_resident_growth"]
list_algorithms = ["simple_screening", "directed_selection_migration", "pair_top_communities", "multiple_pair_top", "Blouin2015", "Mueller2019", "Panke_Buisse2015", "Swenson2000a", "Swenson2000b", "Williams2007a", "Williams2007b", "Wright2019"]
#list_algorithms = ["simple_screening", "direct_selection", "pair_top_communities", "multiple_pair_top", "Blouin2015", "Mueller2019", "Panke_Buisse2015", "Swenson2000a", "Swenson2000b", "Wright2019"]

for j in range(len(list_phenotypes)):
    # Parameters for simulation
    params_simulation.update({"selected_function": list_phenotypes[j]}) # selected function

    # Make the list of algorithms
    algorithms = make_algorithms(params_simulation)
    
    # Simulation
    for i in range(len(list_algorithms)):
        simulate_community(
            assumptions = assumptions,
            params = params,
            dynamics = dynamics,
            params_simulation = params_simulation, 
            params_algorithm = algorithms[algorithms["algorithm_name"] == list_algorithms[i]], 
            write_composition = True,
            file_name = data_directory + list_algorithms[i],
            assembly_type = str(list_algorithms[i]),
        )
        
