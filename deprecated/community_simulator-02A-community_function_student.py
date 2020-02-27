#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This python script contains the community function (phenotype) algorithms written by students.

"""

# Chang-Yu Chang ----
#Sum up the community function for a well
def additive_community_function(plate, sigma = 0.001): #Sigma is the measurement error
    N_tot = plate.N.shape[1]
    
    # Assign additive neutral traits to each of species
    ## Number of species in pools
    S_tot = plate.N.shape[0]
    ## Assign additive traits to each species
    np.random.seed(0)
    traits = np.random.normal(0, 0.1, size=S_tot)

    return(np.sum(plate.N.values*traits[:,None],axis=0)*(1+ np.random.normal(0,sigma,N_tot)))


# Jean Villa ----

def resource_distance_community_function(plate,low=0,high=0.1,sigma = 0.01): #Sigma is the measurement error
   R_tot = plate.R.shape[0]
   well_tot = plate.R.shape[1]
   # Assign additive traits to each species
   if not hasattr(resource_distance_community_function, 'R_target'):
       np.random.seed(0)
       R_target = np.concatenate([np.zeros(1),np.random.uniform(low, high, size=R_tot-1)])
   R_dist = np.sqrt(np.sum((np.tile(R_target,(well_tot,1)) - plate.R.T)**2,axis=1))
   return np.array(R_dist.T)* -1 #(so we select for positive community function)


# Paul Gerald Sanchez ----

def Resource_Environment(plate):

    from scipy.spatial import distance
    
    n = len(plate.R) # number of different resources
    n_wells = len(plate.R.columns) # number of wells
    
    Rstar = np.zeros(n)
    Rstar[15] = 10
    Rstar[18] = 10
    
    distances = np.zeros(n_wells)
    
    for i in range(0,n_wells):
        
        Rwell = plate.R.iloc[:,i]
        distances[i] = -1 * distance.euclidean(Rstar,Rwell) # closer to zero means shorter distance = closer to Rstar
    
    return distances


def Efficiency_Community_Function(plate):
    
    community_function = additive_community_function(plate)
    
    #Read number of wells
    n = len(community_function)
    
    R0_left = plate.R.iloc[0,:]
    R0_left = R0_left.values # remaining R0 after propagation, before passage
    community_function_R0 = community_function * R0_left
    
    return community_function_R0


