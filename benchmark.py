"""
Time the simulate_community
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
colors = sns.color_palette()

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
assumptions = a_default.copy() # Start with default parameters

## Update parameters for community-selection
assumptions.update({
    'SA': 600*np.ones(3), #Number of species in each specialist family (here, 3 families of 60 species)
    'MA': 30*np.ones(3), #Number of resources in each class
    'Sgen': 300, #Number of generalist species (unbiased sampling over alll resource classes)
    "n_wells": 50,
    "m": 0, # Mortality
    "scale": 10**6,  #scale is a conversion factor specifying the number of individual microbial cells present when N = 1.
    "sigma" : 1, # Standard deviation for drawing specifc speices/interaction function
    "alpha": 1, # Scaling factor between species- and interaction-specific function variances
    "l": 0, # Set leakage function to 0 to switch off cross-feeding
    "response": "type III",
    "sigma_max": 5,
    'R0_food': 1000, # Total amount of supplied food
    "rich_medium": True, # Number of food types passed to R0
    # The parameters below will be passed to params_simulation
    "n_propagation": 1, # Length of propagation, or hours within a growth cycle
    "n_transfer": 2, # Number of total transfer, or number of passage
    "n_transfer_selection": 1, # Number of transfer implementing seleciton regimes
    "dilution": 1/1000, # Dilution factor at every transfer
    "n_inoc": 10**6,  #Number of cells sampled from the regional species at start
    "selected_function": "f1_additive"
})
# Prepare experiment setup in this universe
params, params_simulation = prepare_experiment(assumptions, seed = 1)


#
data_directory = "data/test/"
list_algorithms = ["simple_screening"]
list_phenotypes = ["f1_additive"]
i = 0; j = 0
params_simulation.update({"selected_function": list_phenotypes[j]}) # selected function
algorithms = make_algorithms(params_simulation)

params_algorithm = algorithms[algorithms["algorithm_name"] == list_algorithms[i]]
write_composition = True
file_name = data_directory + list_algorithms[i]
assembly_type = str(list_algorithms[i])


import time
start_time = time.time()
temp_time_df = list()

# Print out the algorithms
print("\nAlgorithm: "+ params_algorithm["algorithm_name"][0])
print("\n")
print(params_algorithm[["transfer", "community_phenotype", "selection_algorithm", "migration_algorithm"]].to_string(index = False))

#----------------------------------------------------------------------------------------------------------------
prev_time = start_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("print out algorithm table\t---%s seconds ---" % (time.time() - start_time))
#----------------------------------------------------------------------------------------------------------------


# Set seeds
np.random.seed(2)

# Make initial state
init_state = MakeInitialState(assumptions)

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("initial state\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

# Make plate
plate = Community(init_state, dynamics, params, scale = assumptions["scale"], parallel = True) 
#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("make Community object\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

# Update the community composition by sampling from the pool
print("\nGenerating initial plate")
plate.N = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("make initial plate\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

plate_test = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])
#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("make test plate\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------


plate.N = plate.N + plate_test
#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("migration \t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------




# Update the supplied resource if assumptions["rich_medium"]
if assumptions["rich_medium"]:
    plate.R = make_rich_medium(plate.R, assumptions)
    plate.R0 = make_rich_medium(plate.R, assumptions) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on
    
#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("make rich medium\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

    
# Add the attributes that are essential to the function measurement to the plate objects 
print("\nAdding attributes that are essential to the community function to the plate object")
plate = add_community_function(plate, dynamics, assumptions, params_simulation)

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("add community function\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------


# Empty list for saving data
plate_data_list = list() # Plate composition
community_function_list = list() # Community function

# Save the inocula composition
plate_data = reshape_plate_data(plate, transfer_loop_index = 0, assembly_type = assembly_type, community_function_name = params_algorithm["community_phenotype"][0]) # Initial state
plate_data_list.append(plate_data)

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("reshape plate data\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------


# Save the initial function
community_function = globals()[params_algorithm["community_phenotype"][0]](plate, assumptions = assumptions) # Community phenotype
richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
biomass = list(np.sum(plate.N, axis = 0)) # Biomass
function_data = reshape_function_data(community_function_name = params_algorithm["community_phenotype"][0], community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = 0, assembly_type = assembly_type)        
community_function_list.append(function_data) # Transfer = 0 means that it's before selection regime works upon

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("reshape function data\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------



# Output the plate composition and community functions if write_composition set True
if write_composition == True:
    plate_data.to_csv(file_name + "-" + params_algorithm["community_phenotype"][0] + "-T" + "{:02d}".format(0) + "-composition.txt", index = False)
    function_data.to_csv(file_name + "-" + params_algorithm["community_phenotype"][0] + "-T" + "{:02d}".format(0) + "-function.txt", index = False)


#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("write plate and function data \t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

    
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

    
    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("propagation\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------

    # Append the composition to a list
    plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1, assembly_type = assembly_type, community_function_name = phenotype_algorithm) # Transfer = 0 means that it's before selection regime works upon
    plate_data_list.append(plate_data)

    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("reshape plate data\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------
    
    # Community phenotype, richness, and biomass
    community_function = globals()[phenotype_algorithm](plate, assumptions = assumptions) # Community phenotype
    richness = np.sum(plate.N >= 1/assumptions["scale"], axis = 0) # Richness
    biomass = list(np.sum(plate.N, axis = 0)) # Biomass
    function_data = reshape_function_data(community_function_name = phenotype_algorithm, community_function = community_function, richness = richness, biomass = biomass, transfer_loop_index = i + 1 , assembly_type = assembly_type)        
    community_function_list.append(function_data) # Transfer = 0 means that it's before selection regime works upon

    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("reshape function data\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------

    # Output the plate composition and community functions if write_composition set True
    if write_composition == True:
        plate_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-composition.txt", index = False) # Transfer = 0 means that it's before selection regime works upon
        function_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-function.txt", index = False)
    # Passage and tranfer matrix if is selection experiment
    
    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("write plate and function data\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------


    # Passage and tranfer matrix
    transfer_matrix = globals()[selection_algorithm](community_function)
    plate.Passage(transfer_matrix * params_simulation["dilution"])
    
    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("transfer\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------

    # Migration
    m = globals()[migration_algorithm](community_function) 
    
    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("migration\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------
    
    plate.N = migrate_from_pool(plate, pool = params_simulation["pool"], migration_factor = m, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])

    #----------------------------------------------------------------------------------------------------------------
    prev_time = current_time
    current_time = time.time()
    temp_time_df.append(current_time - prev_time)
    print("migration from pool\t---%s seconds ---" % (time.time() - prev_time))
    #----------------------------------------------------------------------------------------------------------------


    if params_algorithm["algorithm_name"][0] == 'knock_in' and selection_algorithm == 'select_top':
        for k in plate.N.columns:
            s_id = np.random.choice(np.where(plate.N[k]==0)[0])
            plate.N[k][s_id]= 1/params_simulation["dilution"] * 1/assumptions["scale"]
    if params_algorithm["algorithm_name"][0] == 'knock_out' and selection_algorithm == 'select_top':
        for k in plate.N.columns:
            s_id = np.random.choice(np.where(plate.N[k]>0)[0])
            plate.N[k][s_id]=0        
    if params_algorithm["algorithm_name"][0] == 'bottleneck' and selection_algorithm == 'select_top':
        plate.Passage(np.eye(assumptions['n_wells'])* params_simulation["dilution"])
    if params_algorithm["algorithm_name"][0] == 'resource' and selection_algorithm == 'select_top':
        #Remove fresh environment that was added by passage
        plate.R = plate.R - plate.R0
        #change default fresh renvironment so that all subsequent rounds use R0
        for k in plate.R0.columns:
            r_id = np.random.choice(np.where(plate.R0[k]>=0)[0])
            plate.R0[k][r_id] = assumptions['R0_food']/10
        plate.R0 = plate.R0/np.sum(plate.R0)*assumptions['R0_food']
        ##add new fresh environment (so that this round uses R0
        plate.R = plate.R + plate.R0
        
print("\nAlgorithm "+ params_algorithm["algorithm_name"][0] + " finished")


#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("finish loop\t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------


# Concatenate data from from different transfers
plate_data_con = pd.concat(plate_data_list)
community_function_con = pd.concat(community_function_list)


#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("concat plate and function \t---%s seconds ---" % (time.time() - prev_time))
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
prev_time = current_time
current_time = time.time()
temp_time_df.append(current_time - prev_time)
print("\t---%s seconds ---" % (time.time() - start_time))
#----------------------------------------------------------------------------------------------------------------
