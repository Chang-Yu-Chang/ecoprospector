#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

# Design migration_factor (a sequence of binary factors) 
def migration_factor_from_function(community_function, migration_type = 'migrate_half'):
    # Number of wells
    n_wells = len(community_function)
    
    # Migration
    if migration_type == 'migrate_half': # Migrate half 
        migration_factor = [1, 0] * int(n_wells/2)
    elif migration_type == 'random': # Random migration. Each well has probability=0.5 being migrated
        migration_factor = np.random.binomial(1, 0.5, size = n_wells)
    else: # No migration 
        migration_factor = np.zeros(n_wells)

    return migration_factor

# Migrate from species pool to the plate 
def migrate_from_pool(plate, pool, migration_factor):
    # Migration plate
    migration_plate = SampleFromPool(plate.N, pool) * migration_factor
    
    # Migration
    plate_migrated = plate.N + migration_plate 

    return plate_migrated


# Test migration function algorithms
def test_migration_function(migration_algorithm, assumptions):
    n=96
    np.random.seed(0) #you are all going to get the same species pool
 #   assumptions=a_default.copy() #Start with default parameters
    assumptions.update({'n_wells':n, 'c1' :.01, 'muc':0.1, 'm':0}) 
    init_state = MakeInitialState(assumptions)
    params = MakeParams(assumptions)
    species_pool = RegionalSpeciesPool(assumptions) # Generate a species pool

    # Make plate
    np.random.seed(0) # Change to a unique number so that you each start the experiment with a slightly different species pool
    plate1 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community
    plate1.N = SampleFromPool(plate1.N,species_pool) # Populate the well by sampling from the species pool
    
    # Propagate
    plate1.Propagate(24)
    community_function = additive_community_function(plate1)
    m = globals()[migration_algorithm](community_function) # Call function name from selection_function_list
    
    return len(m) == n
