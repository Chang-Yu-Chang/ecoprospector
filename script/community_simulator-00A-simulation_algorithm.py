#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

# Algorithms for simulating community 
def simulate_community(assumptions, phenotype_algorithm, selection_algorithm, migration_algorithm, n = 96, n_propagation = 20, dilution=1/1000, replicate_index = 1, write_composition=False):
    np.random.seed(0) #you are all going to get the same species pool
 #   assumptions=a_default.copy() #Start with default parameters
    assumptions.update({'n_wells':n, 'c1' :.01, 'muc':0.1, 'm':0})  # default setting
    init_state = MakeInitialState(assumptions)
    params = MakeParams(assumptions)
    species_pool = RegionalSpeciesPool(assumptions) # Generate a species pool

    # Make plate
    np.random.seed(replicate_index) # Change to a unique number so that you each start the experiment with a slightly different species pool
    plate1 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community
    plate1.N = SampleFromPool(plate1.N,species_pool) # Populate the well by sampling from the species pool
    
    # Propagate communities
    function_df = list()
    plate_N_df = list()
    
    for k in range(n_propagation):
        # Propagation
        plate1.Propagate(24)
        
        if write_composition == True:
#            plate_N_df = list()
            plate1.N.to_csv("data/passages/" + "{:02d}".format(phenotype_algorithmID) + "-" + phenotype_algorithm + "-" +"{:02d}".format(selection_algorithmID) + "-" + selection_algorithm + "-" +"{:02d}".format(migration_algorithmID) + "-" + migration_algorithm + "-" + "{:02d}".format(rr) + "-T" + "{:02d}".format(k) + "-N" + ".csv", index=False)
            plate1.R.to_csv("data/passages/" + "{:02d}".format(phenotype_algorithmID) + "-" + phenotype_algorithm + "-" +"{:02d}".format(selection_algorithmID) + "-" + selection_algorithm + "-" +"{:02d}".format(migration_algorithmID) + "-" + migration_algorithm + "-" + "{:02d}".format(rr) + "-T" + "{:02d}".format(k) + "-R" + ".csv", index=False)
            print("Species and resource composition recorded in folder data/passage/")

        # Community phenotype
        community_function = globals()[phenotype_algorithm](plate1)
        function_df.append(community_function)
        
        # Passage and tranfer matrix
        if selection_algorithm in ["exponealing", "exponealing2"]: # Propagation time dependent function
            t = globals()[selection_algorithm](community_function, k)
        else:
            t = globals()[selection_algorithm](community_function)
        
        plate1.Passage(t*dilution)
        
        # Migration
        m = globals()[migration_algorithm](community_function) 
        plate1.N = migrate_from_pool(plate1, species_pool, m)
        
        # Print the propagation progress
        print("propagation: " + str(k)) 
    
    return function_df



# Test

def test_simulate_community(assumptions, phenotype_algorithm, selection_algorithm, migration_algorithm, n = 96, n_propagation = 2, dilution=1/1000, replicate_index = 1, write_composition=False):
    np.random.seed(0) #you are all going to get the same species pool
 #   assumptions=a_default.copy() #Start with default parameters
    assumptions.update({'n_wells':n, 'c1' :.01, 'muc':0.1, 'm':0})  # default setting
    init_state = MakeInitialState(assumptions)
    params = MakeParams(assumptions)
    species_pool = RegionalSpeciesPool(assumptions) # Generate a species pool

    # Make plate
    np.random.seed(replicate_index) # Change to a unique number so that you each start the experiment with a slightly different species pool
    plate1 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community
    plate1.N = SampleFromPool(plate1.N,species_pool) # Populate the well by sampling from the species pool
    
    # Propagate communities
    function_df = list()
    plate_N_df = list()
    
    for k in range(n_propagation):
        # Propagation
        plate1.Propagate(24)
        
        if write_composition == True:
#            plate_N_df = list()
#            plate1.N.to_csv("data/passages/" + "{:02d}".format(phenotype_algorithmID) + "-" + phenotype_algorithm + "-" +"{:02d}".format(selection_algorithmID) + "-" + selection_algorithm + "-" +"{:02d}".format(migration_algorithmID) + "-" + migration_algorithm + "-" + "{:02d}".format(rr) + "-T" + "{:02d}".format(k) + ".csv", index=False)
            print("data/passages/" + "{:02d}".format(phenotype_algorithmID) + "-" + phenotype_algorithm + "-" +"{:02d}".format(selection_algorithmID) + "-" + selection_algorithm + "-" +"{:02d}".format(migration_algorithmID) + "-" + migration_algorithm + "-" + "{:02d}".format(rr) + "-T" + "{:02d}".format(k) + ".csv")

        # Community phenotype
        community_function = globals()[phenotype_algorithm](plate1)
        function_df.append(community_function)
        
        # Passage and tranfer matrix
        if selection_algorithm in ["exponealing", "exponealing2"]: # Propagation time dependent function
            t = globals()[selection_algorithm](community_function, k)
        else:
            t = globals()[selection_algorithm](community_function)
        
        plate1.Passage(t*dilution)
        
        # Migration
        m = globals()[migration_algorithm](community_function) 
        plate1.N = migrate_from_pool(plate1, species_pool, m)
        
        # Print the propagation progress
        print("propagation: " + str(k)) 
    
    return plate1, function_df


