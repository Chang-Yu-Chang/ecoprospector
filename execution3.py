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

# Universal parameters
## Thehse are the default functions from community-simulator
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
    'n_wells':96, #Number of independent wells
    'S':100, #Number of species per well (randomly sampled from the pool of size Stot = sum(SA) + Sgen)
    'food':0, #index of food source (when a single resource is supplied externally)
    'R0_food':1000, #unperturbed fixed point for supplied food
    'regulation':'independent', #metabolic regulation (see dRdt)
    'response':'type I', #functional response (see dRdt)
    'supply':'off' #resource supply (see dRdt)
}

## Update parameters for community-selection
assumptions.update({
    "n_wells": 24,
    "c1": 0.01, #Rescale uptake rate part 1. This is needed to avoid numerical errors that slow down the simulations
    "muc": 0.1, # Rescale uptake part 2
    "m": 0, # Mortality
    "scale": 10**9,  #scale is a conversion factor specifying the number of individual microbial cells present when N = 1.
    "sigma" : 1, # Standard deviation for drawing specifc speices/interaction function
    "alpha": 1, # Scaling factor between species- and interaction-specific function variances
    "response": "type I"
})



# Prepare experiment setup in this universe
params, species_pool, species_function, interaction_function = prepare_experiment(assumptions, seed = 1)

list_phenotype = ["f1_additive", "f2_interaction", "f3_additive_binary", "f4_interaction_binary"]

for i in range(len(list_phenotype)):
    plate_df, function_df = simulate_community(
        assumptions = assumptions,
        params = params,
        dynamics = dynamics,
        params_simulation = params_simulation, 
        params_algorithm = pair_top_communities, 
        write_composition = True,
        file_name = "data/pair_top",
        assembly_type = "pair_top",
    )



print(datetime.datetime.now()) # print time










