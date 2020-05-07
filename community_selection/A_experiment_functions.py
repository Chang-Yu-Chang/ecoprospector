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


def sample_from_pool(plate_N, scale = 10**6, inocula = 10**6, migration = False, monoculture = False):
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
    
    # Draw community
    if monoculture == False:
        # Sample initial community for each well
        for k in range(plate_N.shape[1]):
            # For each well, sample community from different microbiome sample
            if migration == False:
                np.random.seed(k + 1) 
            pool = np.random.power(0.01, size = S_tot) # Power-law distribution
            pool = pool/np.sum(pool) # Normalize the pool
            consumer_list = np.random.choice(S_tot, size = inocula, replace = True, p = pool) # Draw from the pool
            my_tab = pd.crosstab(index = consumer_list, columns = "count") # Calculate the cell count
            N0[my_tab.index.values,k] = np.ravel(my_tab.values / scale) # Scale to biomass

        # Make data.frame
        N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

    # Monoculture plate
    elif monoculture == True:
        N0 = np.eye(plate_N.shape[0]) / scale
        N0 = pd.DataFrame(N0, index = consumer_index, columns = ["W" + str(i) for i in range(plate_N.shape[0])])


    return N0




# Migrate from species pool to the plate 
def migrate_from_pool(plate,migration_factor,assumptions,power_law = True,community_function = None):
    '''
    If power_law pool is true than sample n_migration cells from species pool following pair law distribution
    If power_law is false sample s_migration species from isolates with each species starting at a cell count to survive dilution
    If community function is passed in as n argument it will skip the best performing community.
    '''
    if power_law:
        if np.sum(migration_factor) !=0:
            migration_plate = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = assumptions["n_migration"], migration = True) * migration_factor # Migration factor is a list determined by migration algorithms and community function
            plate_migrated = plate.N + migration_plate 
        else:
            plate_migrated = plate.N
    else: 
        plate_migrated = plate.N.copy()
        if community_function is None:
            winning_index = 0
        else:
            winning_index = np.where(community_function >= np.max(community_function))[0][0]
        for k in plate.N.columns:
            if k != plate.N.columns[winning_index] or community_function is None:
                for j in range(0,assumptions['s_migration']):
                    s_id = np.random.choice(np.where(plate.N[k]==0)[0])
                    plate.N[k][s_id]= 1/assumptions["dilution"] * 1/assumptions["scale"] * 1/assumptions['s_migration']
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
    write_composition = False,
    write_function = False):
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
    
    if "monoculture" in params_algorithm["algorithm_name"][0]:
        assumptions.update({"n_wells": np.sum(assumptions["SA"])  + assumptions["Sgen"]})

    # Set seeds
    np.random.seed(2)
    
    # Make initial state
    init_state = MakeInitialState(assumptions)

    # Make plate
    plate = Community(init_state, dynamics, params, scale = assumptions["scale"], parallel = True) 
    
    # Update the community composition by sampling from the pool
    print("\nGenerating initial plate")
    if "monoculture" in params_algorithm["algorithm_name"][0]:
        plate.N = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = params_simulation["n_inoc"], monoculture = True) 
    else:
        plate.N = sample_from_pool(plate.N, scale = assumptions["scale"], inocula = params_simulation["n_inoc"])

    
    # Update the supplied resource if assumptions["rich_medium"]
    if assumptions["rich_medium"]:
        plate.R = make_rich_medium(plate.R, assumptions)
        plate.R0 = make_rich_medium(plate.R, assumptions) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on

    # Add the attributes that are essential to the function measurement to the plate objects 
    print("\nAdding attributes that are essential to the community function to the plate object")
    plate = add_community_function(plate, dynamics, assumptions, params, params_simulation, params_algorithm)

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
    if write_function == True:
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
            if ((i+1) == assumptions["n_transfer_selection"]) or ((i+1) == assumptions["n_transfer"]):
                plate_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-composition.txt", index = False) # Transfer = 0 means that it's before selection regime works upon
        if write_function == True:
            function_data.to_csv(file_name + "-" + phenotype_algorithm + "-T" + "{:02d}".format(i + 1) + "-function.txt", index = False)

        # Passage and tranfer matrix if is selection experiment

        # Passage and tranfer matrix
        transfer_matrix = globals()[selection_algorithm](community_function)
        if "monoculture" in params_algorithm["algorithm_name"][0]:
            plate = passage_monoculture(plate, params_simulation["dilution"])
        else:
            plate.Passage(transfer_matrix * params_simulation["dilution"])
        
        
        # Migration
        m = globals()[migration_algorithm](community_function) 
        plate.N = migrate_from_pool(plate, migration_factor = m, assumptions = assumptions) # By default, n_migration is the same as n_inoc

        # Perturbation section
        if 'bottleneck' in params_algorithm["algorithm_name"][0]  and selection_algorithm == 'select_top':
            winning_index = np.where(community_function >= np.max(community_function))[0][0] 
            dilution_matrix = np.eye(assumptions['n_wells'])
            dilution_matrix[winning_index,winning_index] = 1
            # For additional bottleneck. For example, bottleneck_10 means the total dilution factor is 1000*10 = 10000
            if params_algorithm["algorithm_name"][0] == "bottleneck_1":
                plate.Passage(dilution_matrix / 1)
            if params_algorithm["algorithm_name"][0] == "bottleneck_10":
                plate.Passage(dilution_matrix / 10)
            if params_algorithm["algorithm_name"][0] == "bottleneck_100":
                plate.Passage(dilution_matrix / 100)
        if  'knock_in' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
            selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice
            winning_index = np.where(community_function >= np.max(community_function))[0][0] 
            for k in plate.N.columns:
                if k != plate.N.columns[winning_index]:
                    sel = [x for x in np.where(plate.N[k]==0)[0] if x not in selected] #List to be selected from
                    if 'isolates' in params_algorithm["algorithm_name"][0]:
                        # Knock in an isolate that is 1) not in the resident community 2) have high monoculture function (above 90% isolates in the pool)  
                        temp = np.logical_and(np.array(plate.N[k]==0), plate.isolate_function >= np.percentile(plate.isolate_function, q = 90))
                        sel = [x for x in np.where(temp)[0] if x not in selected] #List to be selected from                
                    if len(sel) == 0:
                        continue
                    s_id = np.random.choice(sel)
                    selected.append(s_id)
                    plate.N[k][s_id]= 1/params_simulation["dilution"] * 1/assumptions["scale"]
        if 'knock_out' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
            selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice=
            winning_index = np.where(community_function >= np.max(community_function))[0][0]
            for k in plate.N.columns:
                if k != plate.N.columns[winning_index]:
                    sel = [x for x in np.where(plate.N[k]>0)[0] if x not in selected] #List to be selected from
                    if len(sel) == 0: #If you've already knocked out every species skip
                        continue
                    s_id = np.random.choice(sel)
                    selected.append(s_id)
                    plate.N[k][s_id]=0      
        if 'resource' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
            selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice
            winning_index = np.where(community_function >= np.max(community_function))[0][0]
            #Remove fresh environment that was added by passage
            plate.R = plate.R - plate.R0
            #change default fresh renvironment so that all subsequent rounds use R0
            for k in plate.R0.columns:
                if k != plate.R0.columns[winning_index]: 
                    #By default pick 2 resources at randomk
                    sel = [(x,y) for x in np.where(plate.R0[k]>=0)[0] for y in np.where(plate.R0[k]>=0)[0] if x !=y]
                    sel = [x for x in sel if x not in selected]
                    if 'add' in params_algorithm["algorithm_name"][0]: #remove from top and add to random
                        top = np.where(plate.R0[k]==np.max(plate.R0[k]))[0]
                        sel = [x for x in sel if x[0] ==top]
                    if 'remove' in params_algorithm["algorithm_name"][0]: #remove from random and add to bottom
                        bottom = np.where(plate.R0[k]==np.min(plate.R0[k]))[0]
                        sel = [x for x in sel if x[1] ==bottom]
                    if len(sel) == 0: #If you've done every possible resource pertubation skip.
                        continue
                    r_id = sel[np.random.choice(len(sel))]
                    selected.append(r_id)   
                    if 'rescale_add' in params_algorithm["algorithm_name"][0]:  # Increase Fraction of resource
                        plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]]*(1+params_simulation['R_percent']) #increase resource conc by fixed %
                    elif 'rescale_remove' in params_algorithm["algorithm_name"][0]: # Decrease Fraction of resource
                        plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-params_simulation['R_percent'])  #decrease resource 
                    elif 'resource_old' in params_algorithm["algorithm_name"][0]:
                        plate.R0[k] = plate.R0[k] * (1-params_simulation['R_percent']) #Dilute old resources
                        plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (assumptions['R0_food']*params_simulation['R_percent'])
                    else: #Default to resource swap.
                        # For resource_add with R_percent 0.1, 0.5, and 1
                        if "resource_add_10" in params_algorithm["algorithm_name"][0]:
                            plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*0.1) #add new resources
                            plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-0.1) #remove new resources
                        elif "resource_add_2" in params_algorithm["algorithm_name"][0]:
                            plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*0.5) #add new resources
                            plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-0.5) #remove new resources
                        elif "resource_add_1" in params_algorithm["algorithm_name"][0]:
                            plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*1) #add new resources
                            plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-1) #remove new resources
                        # Else default remain the same
                        else:
                            plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*params_simulation['R_percent']) #add new resources
                            plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-params_simulation['R_percent']) #remove new resources

            plate.R0 = plate.R0/np.sum(plate.R0)*assumptions['R0_food'] #Keep this to avoid floating point error and rescale when neeeded.
            #add new fresh environment (so that this round uses R0
            plate.R = plate.R + plate.R0
            
    print("\nAlgorithm "+ params_algorithm["algorithm_name"][0] + " finished")

    # Concatenate data from from different transfers
    if (write_composition == False) and (write_function == False):
        plate_data_con = pd.concat(plate_data_list)
        community_function_con = pd.concat(community_function_list)
        return plate_data_con, community_function_con
    else: 
        return None



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
    """
    if resource or migation is in the algorithm name, chage the assemly name
    """

    # Make data.frame
    temp_df = pd.DataFrame({
        "Assembly": np.repeat(assembly_type, number_well),
        "CommunityPhenotypeName": np.repeat(community_function_name, number_well),
        "Well": ["W" + str(i) for i in range(number_well)], 
        "Transfer": np.repeat(str(transfer_loop_index), number_well), 
        "CommunityPhenotype": temp_vector1,
        "Richness": temp_vector2,
        "Biomass": temp_vector3})
    
    # Turn the transfer columns as numeric
    temp_df[["Transfer"]] = temp_df[["Transfer"]].apply(pd.to_numeric)
    
    return temp_df 

 

def add_community_function(plate, dynamics, assumptions, params, params_simulation, params_algorithm):
    """
    Add the function attribute to the community
    
    For f1 and f3, add species_function 
    For f2 and f4, add interaction_function
    For f5, add invasion_plate_t0 and invasion_plate_t1
    For f6, f7, and f8, add resident_plate_t0_N, resident_plate_t1_N, resident_plate_t0_R, and resident_plate_t1_R
    
    if isolates calculate function for every isolate in monoculture.
    """
    np.random.seed(2)
    plate_func = plate.copy()

    # Species function, f1 and f3
    setattr(plate_func, "species_function", params_simulation["species_function"]) # Species function for additive community function

    # Interactive functions, f2 and f4
    setattr(plate_func, "interaction_function", params_simulation["interaction_function"]) # Interactive function for interactive community function
    setattr(plate_func, "interaction_function_p25", params_simulation["interaction_function_p25"])


    # Invasion function f5 and resident function f6, f7 and f8 as well
    if ("invader" in params_simulation["selected_function"]) or ("resident" in params_simulation["selected_function"]):
        # Keep the initial plate R0 for function f7 
        setattr(plate_func, "R0_initial", plate_func.R0)
        
        assumptions_invasion = assumptions.copy()
        params_invasion = params.copy()

        # For invader, only sample one species. For resident community, sample multiple species
        if "invader" in params_simulation["selected_function"]:
            assumptions_invasion.update({"n_wells": np.sum(assumptions["SA"])  + assumptions["Sgen"]})
        elif "resident" in params_simulation["selected_function"]:
            assumptions_invasion.update({"n_wells": 96})

        # Make plates
        init_state_invasion = MakeInitialState(assumptions_invasion)
        plate_invasion = Community(init_state_invasion, dynamics, params_invasion, scale = assumptions_invasion["scale"], parallel = True)
#        print(plate.N.shape)
        print("Sampling invader (resident) community")
        if "invader" in params_simulation["selected_function"]:
            plate_invasion.N = sample_from_pool(plate_invasion.N, scale = assumptions_invasion["scale"], monoculture = True, migration = False) 
        elif "resident" in params_simulation["selected_function"]:
            plate_invasion.N = sample_from_pool(plate_invasion.N, scale = assumptions_invasion["scale"], inocula = assumptions_invasion["n_inoc"]) 
        
        if assumptions["rich_medium"]:
            plate_invasion.R = make_rich_medium(plate_invasion.R, assumptions_invasion)
            plate_invasion.R0 = make_rich_medium(plate_invasion.R, assumptions_invasion) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on 

        # Save the t0 plate
        temp = assumptions_invasion["n_transfer"] - assumptions_invasion["n_transfer_selection"]
        print("\nStabilizing the invader (resident) community. Passage for " + str(temp) + " transfers. The plate has ", str(assumptions_invasion["n_wells"]), " wells")
        plate_invasion_t0 = plate_invasion.copy()
        plate_invasion.Propagate(assumptions_invasion["n_propagation"])
        print("Passaging invader (resident) community. Transfer 1")



        # Grow the invader plate 
        if "invader" in params_simulation["selected_function"]:
            for i in range(3):
                plate_invasion = passage_monoculture(plate_invasion, assumptions_invasion["dilution"])
                plate_invasion.Propagate(assumptions_invasion["n_propagation"])
                print("Passaging invader (resident) community. Transfer " + str(i + 2))

        elif "resident" in params_simulation["selected_function"]:
            for i in range(temp - 1):
                plate_invasion.Passage(np.eye(assumptions_invasion["n_wells"]) * assumptions_invasion["dilution"])
                plate_invasion.Propagate(assumptions_invasion["n_propagation"])
                print("Passaging invader (resident) community. Transfer " + str(i + 2))


        # Save the t1 plate
        plate_invasion_t1 = plate_invasion.copy()
        invasion_plate_growth = np.sum(plate_invasion.N, axis = 0)
        temp_index = np.where(invasion_plate_growth == np.max(invasion_plate_growth))[0][0] # Find the well with the highest biomass

        # Dupliate the chosen community to the whole plate
        temp_df_t0_N = pd.DataFrame(); temp_df_t1_N = pd.DataFrame()
        temp_df_t0_R = pd.DataFrame(); temp_df_t1_R = pd.DataFrame()
        
        for i in range(assumptions["n_wells"]):
            temp_df_t0_N["W" + str(i)] = plate_invasion_t0.N["W" + str(temp_index)]
            temp_df_t1_N["W" + str(i)] = plate_invasion_t1.N["W" + str(temp_index)]
            temp_df_t0_R["W" + str(i)] = plate_invasion_t0.R["W" + str(temp_index)]
            temp_df_t1_R["W" + str(i)] = plate_invasion_t1.R["W" + str(temp_index)]
            
        # Index of the most abundant species in the resident community
        dominant_index = list(np.where(temp_df_t1_N["W0"] == np.max(temp_df_t1_N["W0"]))[0]) 
        print("invader index = " + str(dominant_index))
 
        # Print out the characteristics of (invader) resident community 
        print("Finished passaging the invader (resident) community") 
        print("The community has " + str(np.sum(temp_df_t1_N["W0"] > 0)) + " species")
        print("The community has initial biomass " + str(np.sum(temp_df_t0_N["W0"])) + " and reaches total biomass", str(np.sum(temp_df_t1_N["W0"])))
        print("The invader species (or dominant species in the resident community) has the biomass ", temp_df_t1_N["W0"].iloc[dominant_index][0], " at equilibrium")
    
    
    
        if "invader" in params_simulation["selected_function"]:
            # Add the invasion plate to the attr of community
            setattr(plate_func, "invasion_plate_t0", temp_df_t0_N)
            setattr(plate_func, "invasion_plate_t1", temp_df_t1_N)
        elif "resident" in params_simulation["selected_function"]:
            # Add the resident plate to the attr of community
            setattr(plate_func, "resident_plate_t0_N", temp_df_t0_N)
            setattr(plate_func, "resident_plate_t1_N", temp_df_t1_N)
            setattr(plate_func, "resident_plate_t0_R", temp_df_t0_R)
            setattr(plate_func, "resident_plate_t1_R", temp_df_t1_R)
            ### Firstly are we preserving the paramaters for the isolate plat   
    
    
    # Single isolate functions 
    # Make a plate of all single isolates and measure function    
    if ('isolate' in params_algorithm["algorithm_name"][0]):# or (assumptions["n_migration"] == 1):
        assumptions_isolate = assumptions.copy()
        params_isolate = params.copy()
        
        # For invader, only sample one species. For resident community, sample multiple species
        assumptions_isolate.update({"n_wells": np.sum(assumptions["SA"])  + assumptions["Sgen"]})

        # Make plates
        print("Making monoculture plate")
        init_state_invasion = MakeInitialState(assumptions_isolate)
        plate_isolate = Community(init_state_invasion, dynamics, params_isolate, scale = assumptions_isolate["scale"], parallel = True)
        plate_isolate.N = sample_from_pool(plate_isolate.N, scale = assumptions_isolate["scale"], inocula = assumptions_isolate["n_inoc"], monoculture = True) 
        if assumptions["rich_medium"]:
            plate_isolate.R = make_rich_medium(plate_isolate.R, assumptions_isolate)
            plate_isolate.R0 = make_rich_medium(plate_isolate.R, assumptions_isolate) # R0 for refreshing media on passaging if "refresh_resoruce" is turned on 

        temp = assumptions_isolate["n_transfer"] - assumptions_isolate["n_transfer_selection"]
        print("\nStabilizing monoculture plate. Passage for " + str(temp) + " transfers. The plate has ", str(assumptions_isolate["n_wells"]), " wells")
        plate_isolate.Propagate(assumptions_isolate["n_propagation"])
        print("Passaging monoculture plate. Transfer 1")

        # Grow the monoculture plate 
        for i in range(temp - 1):
            plate_isolate = passage_monoculture(plate_isolate, assumptions_isolate["dilution"])
            plate_isolate.Propagate(assumptions_isolate["n_propagation"])
            print("Passaging monoculture plate. Transfer " + str(i + 2))

        # Save the isolate functions
        isolate_function = np.array(np.sum(plate_isolate.N, axis = 0)) 
        # print(isolate_function)
        # print(len(isolate_function))

        setattr(plate_func, "isolate_function", isolate_function)

    
    return plate_func



def passage_monoculture(plate_mono, f, scale = None, refresh_resource=True):
    """
    Reduced version of function Passage(), for passaging a large set of wells
    """
    self = plate_mono.copy()
    #HOUSEKEEPING
    if scale == None:
        scale = self.scale #Use scale from initialization by default
    self.N[self.N<0] = 0 #Remove any negative values that may have crept in
    self.R[self.R<0] = 0
    
    #DEFINE NEW VARIABLES
    N_tot = np.sum(self.N)
    R_tot = np.sum(self.R)
    N = np.zeros(np.shape(self.N))
    
    self.N = self.N * f 
    plate_mono.N[plate_mono.N <= 1./plate_mono.scale] = 0

    if refresh_resource:
        self.R = self.R * f
        self.R = self.R+self.R0
        
    #In continuous culture, it is useful to eliminate the resources that are
    #going extinct, to avoid numerical instability
    else:
        R_tot = np.sum(self.R)
        R = np.zeros(np.shape(self.R))
        for k in range(self.n_wells):
            if f[k,k] > 0 and R_tot[k] > 0:
                R[:,k] += np.random.multinomial(int(scale*R_tot[k]*f[k,k]),(self.R/R_tot).values[:,k])*1./scale
        self.R = pd.DataFrame(R, index = self.R.index, columns = self.R.keys())
    

    return self
