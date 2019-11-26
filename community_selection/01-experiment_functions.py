#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""

"""
Python functions for simulation in self-assembly, monoculture, and pairwise competition. 
"""

# Make dynanmics by default we will use the microbial consumer resource model
def dNdt(N,R,params):
    return MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(assumptions)(N,R,params)
dynamics = [dNdt,dRdt]


def make_regional_pool(assumptions):
    """
    Create a regional species pool
    """
    # Total number of species (specialist + generalist)
    S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 
    # Assign drawn values based on power-law distribution
    pool = np.random.power(1, size  = S_tot)  # Creating an array with M rows an n_wells columns with 0 entries
    return pool/np.sum(pool)  # Relative species abundance in regional pool


def prepare_experiment(assumptions):
    """
    Prepare the experimental setting shared by self assembly and paiwise competition  
    
    assumptions = dictionary of metaparameters from community-simulator
    
    Return:
    params, species_pool
    """
    np.random.seed(0) 
    
    # Make parameters
    params = MakeParams(assumptions)
    
    # Generate a species pool
    species_pool = make_regional_pool(assumptions) 
    
    return params, species_pool


def sample_from_pool(plate_N, pool, scale=10**6, inocula=10**6):
    """
    Sample communities from regional species pool
    
    """
    N0 = np.zeros((plate_N.shape))
    consumer_index = plate_N.index
    well_names = plate_N.columns
    
    for k in range(plate_N.shape[1]):
        species_list = np.random.choice(len(pool), size=len(pool), replace=True, p=pool)
        my_tab = pd.crosstab(index=species_list, columns="count") # Calculate the biomass count
        N0[my_tab.index.values,k] = np.ravel(my_tab.values / np.sum(my_tab.values) * inocula / scale) # Scale to sum

    # Make data.frame
    N0 = pd.DataFrame(N0,index=consumer_index,columns=well_names)
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

    Return:
    N0 = initial consumer populations
    """
    # Stopifnot
    assert max(species_list) <= len(species_pool), "Some species in the list are not in the pool."
    assert len(species_list) >= number_species, "Cannot make pair from one species."
    assert any(list((sum(x) == 1 for x in initial_frequency))), "Sum of initial frequencies is not equal to 1."
    assert any(list((len(x) == number_species for x in initial_frequency))), "Length of initial frequencies is not equal to number of species."
    
    # All possible combinations of species for given number of species added to a well
    from itertools import combinations
    consumer_pairs = list(combinations(species_list, number_species))
    
    # Extract parameters from assumption
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    F = len(assumptions['SA'])
    
    # Construct lists of names of resources, consumers, resource types, consumer families and wells:
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

## Simulate community
def simulate_community(plate, assumptions, params_simulation, file_name = "data/self_assembly-community", write_composition = False):
    """
    Simulate community dynamics by given experimental regimes
    
    plate = plate
    assumptions = dictionary of metaparameters from community-simulator
    params_simulation = dictionary of parameters for running experiment
    
    return
    
    """
    np.random.seed(0) # Global random seed (i.e all participants)

    # Save the inocula composition
    plate_data = reshape_plate_data(plate, transfer_loop_index = 0) 
    plate_data.to_csv(file_name + "-T" + "{:02d}".format(0) + ".txt", index = False)

    # Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
        
        # Save composition to an empty data
        if write_composition == True:
            # Convert commnity_function into a df and melt the df
            plate_data = reshape_plate_data(plate, transfer_loop_index = i + 1)
            plate_data.to_csv(file_name + "-T" + "{:02d}".format(i+1) + ".txt", index = False)
            print("Species and resource composition saved")

        # Transfer/passage by usigg transfer matrix. For simple passage, the transfer matrix is an identity matrix
        transfer_matrix = np.eye(plate.N.shape[1]) * params_simulation["dilution"]
        plate.Passage(transfer_matrix)
        
        # Print the propagation progress
        print("propagation: " + str(i+1)) 
    

## Higher function that simlate pairs




## Reshape the plate resource and consumer matrix for saving into a txt file
def reshape_plate_data(plate, transfer_loop_index):
    # Temporary function for adding variables to and melting df
    def melt_df(plate_df, data_type = "consumer"):
        # Consumers
        temp_df = pd.DataFrame(plate_df)
        total_number = temp_df.shape[0]
        
        ## Add variables
        temp_df["Type"] = np.repeat(data_type, total_number)
        temp_df["ID"] = range(total_number)
        temp_df["Transfer"] = np.repeat(str(transfer_loop_index), total_number)
        temp_df["Assembly"] = np.repeat("self-assembly", total_number)
        
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

