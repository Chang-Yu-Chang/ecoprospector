#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""

"""
Python functions for computing the community-level phenotypes
"""

import numpy as np
import scipy as sp

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



# Consumer community composition
def community_function_additive(plate, species_function = np.random.normal(0, 1, size = 210)):
    """
    Additive community function with saturation
    
    plate = plate object from package
    species_function = an 1-D array with the length of the size of species pool
    k = an 1-D array of saturation factors. set k = np.zeros(n) for binary function (species presence or absense)
    """
    assert len(species_function) == plate.N.shape[0], "Length of species_function does not match species number in plate."

    community_function = np.sum(plate.N.values * species_function[:,None], axis = 0)
    
    return community_function


def community_function_additive_saturation(plate, species_function = np.random.normal(0, 1, size = 210), k = np.zeros(210)):
    """
    Additive community function with saturation
    
    plate = plate object from package
    species_function = an 1-D array with the length of the size of species pool
    k = an 1-D array of saturation factors. set k = np.zeros(n) for binary function (species presence or absense)
    """
    assert len(species_function) == plate.N.shape[0], "Length of species_function does not match species number in plate."
    assert len(k) == plate.N.shape[0], "Length of k does not match species number in plate."
    
    community_function = np.nansum((plate.N.values * (species_function[:,None]) / (plate.N.values + k[:, None])), axis = 0)

    return community_function

def community_function_complex(plate, species_function):
    """
    Complex community function
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    """
    assert len(species_function) == plate.N.shape[0], "Length of species_function does not match species number in plate."
    
    # Number of species in the pool 
    S_tot = plate.N.shape[0]
    
    # Additive term; diagonal of the species function matrix
    additive_term = np.sum(plate.N.values * np.diag(species_function)[:,None], axis=0)

    # Interaction term; matrix multiplication
    community_composition = np.array(plate.N.iloc[:,0]).reshape(S_tot, 1)
    community_composition_square = np.multiply(community_composition, community_composition.reshape(1, S_tot))
    
    interaction_term = np.sum(community_composition_square * species_function)

    
    return additive_term + interaction_term


def community_function_complex_saturation(plate, species_function, k):
    """
    Complex community function
    
    plate = plate object from package
    species_function = a n by n 2-D array; n is the size of species pool
    k = an 2-D array of saturation factors. set k = np.zeros([n, n]) for binary function (species presence or absense)
 
    """
    assert species_function.shape[0] == plate.N.shape[0], "Length of species_function does not match species number in plate."
    
    # Number of species in the pool 
    S_tot = plate.N.shape[0]
    
    # Additive term; diagonal of the species function matrix
    additive_term = np.nansum((plate.N.values * np.diag(species_function)[:,None]) / (plate.N.values + np.diag(k)[:, None]), axis = 0)

    # Interaction term; matrix multiplication
    community_composition = np.array(plate.N.iloc[:,0]).reshape(S_tot, 1)
    community_composition_square = np.multiply(community_composition, community_composition.reshape(1, S_tot))

    interaction_term = np.nansum(community_composition_square * species_function / ( community_composition_square + k))

    
    return additive_term + interaction_term



# Resource function
def resource_additive(plate, resource_function):
    """
    Additive community function with saturation
    
    plate = plate object from package
    resource_function = an 1-D n-length array; n is the number of available resources
    k = an 1-D array of saturation factors. set k = np.zeros(n) for binary function (species presence or absense)
    """
    assert len(resource_function) == plate.R.shape[0], "Length of resource_function does not match species number in plate."

    community_function = np.sum(plate.R.values * resource_function[:,None], axis = 0)
    
    return community_function



def resource_additive_saturation(plate, resource_function, k):
    """
    Additive community function with saturation
    
    plate = plate object from package
    resource_function = an 1-D n-length array; n is the number of available resources
    k = an 1-D array of saturation factors. set k = np.zeros(n) for binary function (species presence or absense)
    """
    assert len(resource_function) == plate.R.shape[0], "Length of resource_function does not match species number in plate."
    assert len(k) == plate.N.shape[0], "Length of k does not match species number in plate."
    
    community_function = np.nansum((plate.R.values * (resource_function[:,None]) / (plate.R.values + k[:, None])), axis = 0)

    return community_function


def resource_additive(plate, resource_function):
    """
    Complex community function
    
    plate = plate object from package
    resource_function = a n by n 2-D array; n is the number of available resources
    """
    assert len(resource_function) == plate.R.shape[0], "Length of resource_function does not match species number in plate."
    
    # Number of species in the pool 
    R_tot = plate.R.shape[0]
    
    # Additive term; diagonal of the species function matrix
    additive_term = np.sum(plate.R.values * np.diag(resource_function)[:,None], axis=0)

    # Interaction term; matrix multiplication
    community_composition = np.array(plate.R.iloc[:,0]).reshape(R_tot, 1)
    community_composition_square = np.multiply(community_composition, community_composition.reshape(1, R_tot))
    
    interaction_term = np.sum(community_composition_square * resource_function)
    
    return additive_term + interaction_term



def resource_complex_saturation(plate, resource_function, k):
    """
    Complex community function
    
    plate = plate object from package
    resource_function = a n by n 2-D array; n is the number of available resources
    k = an 2-D array of saturation factors. set k = np.zeros([n, n]) for binary function (species presence or absense)
 
    """
    assert resource_function.shape[0] == plate.R.shape[0], "Length of resource_function does not match species number in plate."
    
    # Number of species in the pool 
    R_tot = plate.R.shape[0]
    
    # Additive term; diagonal of the species function matrix
    additive_term = np.nansum((plate.R.values * np.diag(resource_function)[:,None]) / (plate.R.values + np.diag(k)[:, None]), axis = 0)

    # Interaction term; matrix multiplication
    community_composition = np.array(plate.R.iloc[:,0]).reshape(R_tot, 1)
    community_composition_square = np.multiply(community_composition, community_composition.reshape(1, R_tot))

    interaction_term = np.nansum(community_composition_square * resource_function / ( community_composition_square + k))

    
    return additive_term + interaction_term



