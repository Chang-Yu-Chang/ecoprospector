#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp
import pandas as pd
import random

from community_simulator import *
from community_simulator.usertools import *
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_perturbation_algorithms import *


def reshape_plate_data(plate, params_simulation,transfer_loop_index):
    """
    Reshape the plate resource and consumer matrices (wider form) into a melted data.frame (longer form)
    """
    # Temporary function for adding variables to and melting df
    def melt_df(plate_df, data_type = "consumer"):
        # Consumers
        temp_df = pd.DataFrame(plate_df)
        total_number = temp_df.shape[0]
        
        ## Add variables
        temp_df["Type"] = np.repeat(data_type, total_number)
        temp_df["ID"] = range(total_number)
        temp_df["Transfer"] = np.repeat(str(transfer_loop_index), total_number)
        temp_df["exp_id"] = np.repeat(params_simulation['exp_id'] , total_number)

        ## Melt the df
        temp_df = pd.melt(temp_df, id_vars = ["exp_id","Transfer", "Type", "ID"], var_name = "Well", value_name = "Abundance")
        temp_df = temp_df[temp_df.Abundance != 0] # Remove zero abundances
        return temp_df
        
    # Melt the df
    temp_plate = plate.copy() # Copy the original plate 
    df_N = melt_df(temp_plate.N, data_type = "consumer")
    df_R = melt_df(temp_plate.R, data_type = "resource")
    df_R0 = melt_df(temp_plate.R0,data_type = "R0")
    
    # Concatenate dataframes
    merged_df = pd.concat([df_N, df_R,df_R0]) 
    merged_df["Index"] = list(range(0, merged_df.shape[0]))
    merged_df.set_index("Index", inplace = True)

    return merged_df # Return concatenated dataframe


def reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index):
    """
    Reshape the community function, richness, biomass into a melted data.frame
    """
    temp_vector1 = community_function.copy()
    temp_vector2 = richness.copy()
    temp_vector3 = biomass.copy()
    
    # Number of wells
    number_well = len(richness)

    # Make data.frame
    temp_df = pd.DataFrame({
        "exp_id": np.repeat(params_simulation['exp_id'], number_well),
        "Well": ["W" + str(i) for i in range(number_well)], 
        "Transfer": np.repeat(str(transfer_loop_index), number_well), 
        "CommunityPhenotype": temp_vector1,
        "Richness": temp_vector2,
        "Biomass": temp_vector3})
    
    # Turn the transfer columns as numeric
    temp_df[["Transfer"]] = temp_df[["Transfer"]].apply(pd.to_numeric)
    
    return temp_df 
    

def migrate_from_pool(plate,migration_factor,params_simulation, power_law = True, n = None):
    """
    Migrate from species pool to the plate mainly for directed selection)
    If power_law pool is true than sample n cells from species pool following power law distribution (default is same as inoculum)
    If power_law is false sample s_migration species from isolates with each total number of cells equivalent to n
    """
    from community_selection.usertools import sample_from_pool
    if n is None:
        n = params_simulation['n_migration']
    if power_law:
        if np.sum(migration_factor) !=0:
            temp_params_simulation = params_simulation.copy() 
            migration_plate = sample_from_pool(plate.N, params_simulation,n=n) * migration_factor # Migration factor is a list determined by migration algorithms and community function
            plate_migrated = plate.N + migration_plate 
        else:
            plate_migrated = plate.N
    else: 
        if np.sum(migration_factor) !=0:
            migration_plate = plate.N.copy()
            migration_plate[:]  =0
            for k in plate.N.columns:
                if migration_factor[np.where(plate.N.columns == k)[0]]>0: 
                    for j in range(0,params_simulation['s_migration']):
                        s_id = np.random.choice(np.where(plate.N[k]==0)[0])
                        migration_plate[k][s_id]= n * 1/params_simulation["scale"] * 1/params_simulation['s_migration']
            plate_migrated = plate.N + migration_plate
        else:
            plate_migrated = plate.N
    return plate_migrated


def passage_monoculture(plate_mono, f, scale = None, refresh_resource=True):
    """
    Reduced version of Passage(), for passaging a large set of wells without multinomial sampling
    Most code adapted from community-simulator
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
    
    #Poisson sample cells
    self.N = self.N * f *scale
    self.N.applymap(np.random.poisson)   
    self.N = self.N/scale

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


def simulate_community(params, params_simulation, params_algorithm, plate):
    """
    Simulate community dynamics by given experimental regimes
    
    params = parameter passed from community-simulator
    params_simulation = dictionary of parameters for running experiment
    params_algorithm = dictionary of algorithms that determine the selection regime, migration regime, and community pheotypes
    plate = Plate object specified by community-simulator
    
    Return:
    community_composition = concatenated, melted panda dataframe of community and resource composition in each transfer
    community_function = melted panda dataframe of community function
    """

    # Test the community function
    try:
        community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
    except:
        print('/n Community phenotype test failed')
        raise SystemExit

    # Save the inocula composition
    if params_simulation['save_composition']:
        plate_data_list = list() # Plate composition
        plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=0)  # Initial state
        plate_data_list.append(plate_data)
        composition_filename = params_simulation['output_dir'] + params_simulation['exp_id'] + '_composition.txt'   
        
    # Save the initial community function + richness + biomass
    if params_simulation['save_function']:
        community_function_list = list() # Plate composition
        richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
        biomass = list(np.sum(plate.N, axis = 0)) # Biomass
        function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index =0)
        community_function_list.append(function_data)
        function_filename = params_simulation['output_dir'] + params_simulation['exp_id'] + '_function.txt'   


    print("\nStart propogation")
    # Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        # Algorithms used in this transfer
        phenotype_algorithm = params_algorithm["community_phenotype"][i]
        selection_algorithm = params_algorithm["selection_algorithm"][i]
#        migration_algorithm = params_algorithm["migration_algorithm"][i]
        print("Transfer " + str(i+1))

        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
    
        # Measure Community phenotype
        community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
        
        # Append the composition to a list
        if params_simulation['save_composition'] and ((i+1) % params_simulation['composition_lograte'] == 0):
            plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=i+1)  # Initial state
            plate_data_list.append(plate_data)

        if params_simulation['save_function'] and ((i+1) % params_simulation['function_lograte'] == 0):
            richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
            biomass = list(np.sum(plate.N, axis = 0)) # Biomass
            function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index =i+1)
            community_function_list.append(function_data)

		#Store prior state before passaging (For coalescence)
        setattr(plate, "prior_N", plate.N)
        setattr(plate, "prior_R", plate.R)
        setattr(plate, "prior_R0", plate.R0)

        # Passage and tranfer matrix
        transfer_matrix = globals()[selection_algorithm](community_function)
        if params_simulation['monoculture']:
            plate = passage_monoculture(plate, params_simulation["dilution"])
        else:
            plate.Passage(transfer_matrix * params_simulation["dilution"])
        
        # Migration
        # m = globals()[migration_algorithm](community_function) 
        # plate.N = migrate_from_pool(plate, migration_factor = m, params_simulation = params_simulation) # By default, n_migration is the same as n_inoc
        
        # Perturbation
        if (i+1) % params_simulation['n_transfer_selection'] == 0 and params_simulation['directed_selection'] and params_simulation['n_transfer'] != (i+1):
            keep = 0 # Always keep the first communtiy, by the design of selection transfer matrix, this is always the top of selected communities
            plate = perturb(plate,params_simulation, keep = keep)

    if params_simulation['save_composition']:
        pd.concat(plate_data_list).to_csv(composition_filename ,index=False)
    if params_simulation['save_function']:
        pd.concat(community_function_list).to_csv(function_filename ,index=False)
    print("\n"+ params_simulation["exp_id"]+ " finished")

