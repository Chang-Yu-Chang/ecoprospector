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

def select_top_nth(community_function, n):
    """
    Select the top nth single community. Designed for perturbation effect
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)[::-1]
    cut_off = sorted_community_function[n-1] # The top nth
    winner_index = np.where(community_function == cut_off)[0]

    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * n_wells # Old wells
        
    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
  
    return transfer_matrix

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
    randomized_community_function = community_function
    np.random.shuffle(randomized_community_function)
    sorted_community_function = np.sort(randomized_community_function)
    cut_off = sorted_community_function[int(np.round(len(randomized_community_function)*(1-p))) - 1]
    winner_index = np.where(randomized_community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
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
    randomized_community_function = community_function
    np.random.shuffle(randomized_community_function)
    sorted_community_function = np.sort(randomized_community_function)
    cut_off = sorted_community_function[int(np.round(len(randomized_community_function)*(1-p))) - 1]
    winner_index = np.where(randomized_community_function > cut_off)[0][::-1] # Reverse the list so the higher 
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

### Sub-lineage algorithms
def Arora2019(community_function, n_rep=3):
    """
    Arora2019
    Sub-divide wells of plate into lines where each 'line' consists of n_rep communities' 
    Each round the highest function member of the line is used to colonize the next three triplicate wells of that line
    """
    n_wells = len(community_function)
    n_lines = int(np.ceil(n_wells/n_rep)) #Number of lines
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix
def Arora2019_control(community_function, n_rep=3):
    """
  	Same as Arora2019 except the line member is selected at Random
    """
    n_wells = len(community_function)
    n_lines = int(np.ceil(n_wells/n_rep)) #Number of lines
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
  	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
  	    cut_off = np.max(sorted_community_function)
  	    winner_index = np.random.randint(0,n_rep)
  	    if winner_index+i*n_rep >= n_wells:
  	    	  corrected_n_rep  = n_wells % n_rep
  	    	  winner_index = np.random.randint(0,corrected_n_rep)
  	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix
def Raynaud2019a(community_function, n_lines=3):
    """
    Raynaud2019a
    Sub-divide wells of plate into n_lines' 
    Each round the highest function member of the line is used to colonize the  wells of that lineage
    """
    n_wells = len(community_function)
    n_rep  = int(np.ceil(n_wells/n_lines)) #Number of replicates per line
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix
def Raynaud2019a_control(community_function, n_lines=3):
    """
	Same as Raynaud2019a except the lineage member is selected at Random
    """
    n_wells = len(community_function)
    n_rep  = int(np.ceil(n_wells/n_lines)) #Number of replicates per line
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.random.randint(0,n_rep)
	    if winner_index+i*n_rep >= n_wells:
	    	  corrected_n_rep  = n_wells % n_rep
	    	  winner_index = np.random.randint(0,corrected_n_rep)
	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix
def Raynaud2019b(community_function, n_lines=3):
    """
    same as Raynaud2019a except top from each lineage is pooled
    """
    n_wells = len(community_function)
    n_rep  = int(np.ceil(n_wells/n_lines)) #Number of replicates per line
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[:, winner_index+i*n_rep] = 1
    return transfer_matrix
def Raynaud2019b_control(community_function, n_lines=3):
    """
	Same as Raynaud2019b except the lineage member is selected at Random
    """
    n_wells = len(community_function)
    n_rep  = int(np.ceil(n_wells/n_lines)) #Number of replicates per line
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.random.randint(0,n_rep)
	    if winner_index+i*n_rep >= n_wells:
	    	  corrected_n_rep  = n_wells % n_rep
	    	  winner_index = np.random.randint(0,corrected_n_rep)
	    transfer_matrix[:, winner_index+i*n_rep] = 1
    return transfer_matrix


# Bottleneck sizes

def select_top_25_bottleneck_10(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.1
def select_top_25_bottleneck_100(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.01
def select_top_25_bottleneck_1000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.001
def select_top_25_bottleneck_10000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.0001
def pool_top_25_bottleneck_10(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.1
def pool_top_25_bottleneck_100(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.01
def pool_top_25_bottleneck_1000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.001
def pool_top_25_bottleneck_10000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.0001


def pool_top_25_bottleneck_100000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.00001
def pool_top_25_bottleneck_1000000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1.0
    return transfer_matrix*0.000001
def select_top_25_bottleneck_100000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.00001
def select_top_25_bottleneck_1000000(community_function, p=0.25):
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1.0
    return transfer_matrix*0.000001

