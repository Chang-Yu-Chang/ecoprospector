#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp


def f1_additive(plate, params_simulation):
    """
    Additive community function(F1)
    
    plate = plate object from package
    k = an 1-D array of saturation factors. set k = np.zeros(n) for binary function (species presence or absense)
    """
    
    community_function = np.sum(plate.N.values * plate.species_function[:,None], axis = 0)
    
    return community_function


def f2_interaction(plate, params_simulation):
    """
    Additive community function with interaction (F2)
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    """

    # Number of species in the pool 
    S_tot = plate.N.shape[0]
    
    # Additive term
    additive_term = np.sum(plate.N.values * plate.species_function[:,None], axis = 0)
    
    # Interaction term
    interaction_term = np.zeros(plate.N.shape[1])
    for i in range(plate.N.shape[1]): # For each community
        community_composition = np.array(plate.N.iloc[:,i]).reshape(S_tot, 1)
        community_composition_square = np.multiply(community_composition, community_composition.reshape(1, S_tot))
        interaction_term[i] = np.sum(community_composition_square * plate.interaction_function)

    return additive_term + interaction_term


def f2a_interaction(plate, params_simulation):
    """
    Additive community function with interaction (F2)
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    """

    # Number of species in the pool 
    S_tot = plate.N.shape[0]
    
    # Additive term
    additive_term = np.sum(plate.N.values * plate.species_function[:,None], axis = 0)
    
    # Interaction term
    interaction_term = np.zeros(plate.N.shape[1])
    for i in range(plate.N.shape[1]): # For each community
        community_composition = np.array(plate.N.iloc[:,i]).reshape(S_tot, 1)
        community_composition_square = np.multiply(community_composition, community_composition.reshape(1, S_tot))
        interaction_term[i] = np.sum(community_composition_square * plate.interaction_function_p25)

    return additive_term + interaction_term


def f3_additive_binary(plate, params_simulation):
    """
    Complex community function
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    """
    # Binary function using type III response
    plate_temp = plate.copy()
    n = 10; Sm = 1
    plate_temp.N = plate_temp.N / params_simulation["binary_threshold"]
    plate_temp.N = plate_temp.N**n / (1 + plate_temp.N**n/Sm) 
    community_function = np.sum(plate_temp.N.values * plate_temp.species_function[:,None], axis = 0)

    return community_function


def f4_interaction_binary(plate, params_simulation):
    """
    Complex community function
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    k = an 2-D array of saturation factors. set k = np.zeros([n, n]) for binary function (species presence or absense)
 
    """
    # Number of species in the pool 
    S_tot = plate.N.shape[0]

    # Binary function using type III response
    plate_temp = plate.copy()
    n = 10; Sm = 1
    plate_temp.N = plate_temp.N / params_simulation["binary_threshold"]
    plate_temp.N = plate_temp.N**n / (1 + plate_temp.N**n/Sm) 
    
    # Additive term
    additive_term = np.sum(plate_temp.N.values * plate_temp.species_function[:,None], axis = 0)
    
    # Interaction term
    interaction_term = np.zeros(plate_temp.N.shape[1])
    for i in range(plate_temp.N.shape[1]): # For each community
        community_composition = np.array(plate_temp.N.iloc[:,i]).reshape(S_tot, 1)
        community_composition_square = np.multiply(community_composition, community_composition.reshape(1, S_tot))
        interaction_term[i] = np.sum(community_composition_square * plate_temp.interaction_function)
    
    return additive_term + interaction_term


def f5_invader_growth(plate, params_simulation):
    """
    Community function in which an indentical alien community (single or multiple species) invades the selected resident communities.
    This community function is the ratio between the biomass when invader grows with the community and when invader grows alone.
    The biomass of invader growing alone (plate.invasion_plate_t1) should have been included in the plate object attribute.
    
    """
    # Number of species and community
    S_tot = plate.N.shape[0]
    n_wells = plate.N.shape[1]
    
    # Dilute the tested communities
    plate_test = plate.copy()
    
    # Coalesce the tested community with invasion community (or single invader) 
    plate_test.N = plate_test.N + plate.invader_N 
    
    plate_test.Passage(np.eye(n_wells) * params_simulation["dilution"])

    # Grow the coalesced communities
    plate_test.Propagate(params_simulation["n_propagation"])
    
    # Calculate the function by dividing the final x(t) with x(o) of pathogen 
    temp = plate.invasion_plate_t1["W0"]
    dominant_index = list(np.where(temp == np.max(temp))[0]) # Index of the most abundant speceis in the resident community
    invader_growth_along = np.sum(plate.invasion_plate_t1.iloc[dominant_index], axis = 0)
    invader_growth_together = np.sum(plate_test.N.iloc[dominant_index], axis = 0)
    
    #
    function_invader_suppressed_growth = invader_growth_along / invader_growth_together

    return function_invader_suppressed_growth


# Compute the distances from the target resource 
# This function is from Jean
def resource_distance_community_function(plate,R_target,sigma = 0.01): # Sigma is the measurement error
    R_tot = plate.R.shape[0]
    well_tot = plate.R.shape[1]
    relative_resource = np.array(plate.R) #Load plate resource data
    relative_resource[0,:]  = 0.0 #Set supplied resource to 0
    relative_resource = relative_resource/relative_resource.sum(0)  #Look at relative abundance of remaining resource
    R_dist = np.sqrt(np.sum(np.array((np.tile(R_target,(well_tot,1)) - relative_resource.T)**2)[:,1:],axis=1))
    return (np.array(R_dist.T)* -1) * (1+ np.random.normal(0,sigma,well_tot))#(so we select for positive community function)


