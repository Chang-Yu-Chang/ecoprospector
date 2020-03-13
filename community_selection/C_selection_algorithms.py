#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 27 2019
@author: changyuchang
"""

"""
Python functions for selection algorithms during plate passage in 96-well plates.
"""

import numpy as np
import scipy as sp
from functools import partial
  
def no_selection(community_function):
    """
    Direct well-to-well transfer without selection
    """
    # Read number of wells 
    n_wells = len(community_function)
    return np.eye(n_wells)

 
def select_top(community_function):
    """
    Select the top community 
    """
    # Read number of wells 
    n_wells = len(community_function)
    
    # Winner wells
    winner_index = np.where(community_function >= np.max(community_function))[0][::-1] # Reverse the list so the higher 
    
    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * n_wells # Old wells
        
    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
  
    return transfer_matrix

def select_top_dog(community_function):
  """
  100 communities. Reproduce the best one to 60-70 newborns, and reproduce the second best to 30-40 newborns. 
  """
  n_wells = len(community_function)
  sorted_community_function = np.sort(community_function)
  cut_off = sorted_community_function[int(np.round(len(community_function)*0.5)) - 1]
  winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
  
  # Transfer matrix
  transfer_matrix = np.zeros((n_wells,n_wells))
  t_new = range(n_wells) # New wells
  # The best performed community
  t_old = [list(winner_index)[0]] * int(0.6 * n_wells) + [list(winner_index)[1]] * int(0.5 * n_wells) # Old wells
  
  # Fill in the transfer matrix
  for i in range(n_wells):
      transfer_matrix[t_new[i], t_old[i]] = 1
  
  return transfer_matrix

# Make selection algorithms with similar names 
## Select top n%
def temp_select_top(community_function, p):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix

for i in [10, 15, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['select_top%spercent' %i] = partial(temp_select_top, p = i/100)

## Select top n% control
def temp_select_top_control(community_function, p):
    n_wells = len(community_function)
    sorted_community_function = community_function[list(np.argsort(np.random.uniform(size = n_wells)))] # Randomize function
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix

for i in [10, 15, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['select_top%spercent_control' %i] = partial(temp_select_top_control, p = i/100)


## Pooling
def temp_pool_top(community_function, p):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
    return transfer_matrix

for i in [10, 15, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['pool_top%spercent' %i] = partial(temp_pool_top, p = i/100)

## Pooling control
def temp_pool_top_control(community_function, p):
    n_wells = len(community_function)
    sorted_community_function = community_function[list(np.argsort(np.random.uniform(size = n_wells)))] # Randomize function
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
    return transfer_matrix

for i in [10, 15, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['pool_top%spercent_control' %i] = partial(temp_pool_top_control, p = i/100)



# Other algorithms
def Williams2007a(community_function):
    """
    Williams2007a
    Select the top community and impose an bottleneck
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    winner_index = np.where(community_function == np.max(community_function))[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * n_wells # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 10**(-4) # An additional strong bottleneck
    return transfer_matrix


def Williams2007b(community_function, p = 0.2):
    """
    Williams2007b
    Select and pool the top 20% community and impose an bottleneck
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1]
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 10**(-4) # An additional strong bottleneck
    return transfer_matrix


# All directed selection algorithm: keep the top and perturb
def pair_top(community_function):
    """
    Pair the top communities. Each pairwise combination has roughly two replicates
    """
    import itertools

    # Read number of wells
    n_wells = len(community_function)

    # Compute the cutoff based on the number of wells
    cut_off_percent = (np.sqrt(n_wells))/n_wells

    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)

    # Community function value cutoff for selecting communities    
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-cut_off_percent)))]

    # Winner wells
    winner_index = np.where(community_function >= cut_off)[0] # Reverse the list so the higher 
    pairs_list = list(itertools.combinations(winner_index, 2)) # Pair list based on the winer wells

    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) + pairs_list * (int(np.round(1/cut_off_percent)) + 1) # Old wells
 
    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1

    return transfer_matrix

def coalescence(community_function):
    """
    Select the top community and coalesce with all other communities (including self)
    """
    # Read number of wells 
    n_wells = len(community_function)
    
    # Winner wells
    winner_index = np.where(community_function >= np.max(community_function))[0][::-1] # Reverse the list so the higher 
    
    # Transfer matrix
    transfer_matrix = np.eye(n_wells)
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * n_wells # Old wells
        
    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix


def directed_selection_select(community_function):
    """
    Select an keep the top communities, and coalesce top community replicates with migrant communities in the rest of the well
    
    This is only the transfer part of the selection algorithm. It has to work with the migration algorithm `direct_selection_migrate()`
    """
    # Number of cells
    n = len(community_function)
    
    # Compute the cutoff based on the number of wells
    cut_off_percent = (np.sqrt(n))/n

    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)

    # Community function value cutoff for selecting communities    
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-cut_off_percent)))]

    # Winner wells
    winner_index = np.where(community_function >= cut_off)[0][::-1]
    
    # Transfer matrix
    transfer_matrix = np.zeros((n,n))
    t_new = range(n) # New wells
    t_old = list(winner_index) + list(winner_index) * (int(np.round(1/cut_off_percent)) + 1) 
    
    # Fill in the transfer matrix
    for i in range(n):
        transfer_matrix[t_new[i], t_old[i]] = 1
    
    return transfer_matrix


