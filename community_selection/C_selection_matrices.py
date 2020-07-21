#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 27 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp
from functools import partial


def no_selection(community_function):
    """
    Direct well-to-well transfer without selection
    """
    n_wells = len(community_function)
    return np.eye(n_wells)


def temp_select_top(community_function, p):
    """
    Select top p% of communities. Round up to the nearest integer.
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.floor(len(community_function)*(1-p)))]
    winner_index = np.where(community_function >= cut_off)[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells)
    t_old = list(winner_index) * (int(np.ceil(1/p) + 1))
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix

for i in [10, 15, 16, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['select_top%spercent' %i] = partial(temp_select_top, p = i/100)



def temp_select_top_control(community_function, p):
    """
    Select random p% of communities. Round up to the nearest integer.
    """
    n_wells = len(community_function)
    randomized_community_function = community_function.copy()
    np.random.shuffle(randomized_community_function)
    sorted_community_function = np.sort(randomized_community_function)
    cut_off = sorted_community_function[int(np.floor(len(randomized_community_function)*(1-p)))]
    winner_index = np.where(randomized_community_function >= cut_off)[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells)
    t_old = list(winner_index) * (int(np.ceil(1/p)+1))
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix

for i in [10, 15, 16, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['select_top%spercent_control' %i] = partial(temp_select_top_control, p = i/100)


def temp_pool_top(community_function, p):
    """
    Pool top p% of communities. Round up to the nearest integer.
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.floor(len(community_function)*(1-p)))]
    winner_index = np.where(community_function >= cut_off)[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, list(winner_index)] = 1
    return transfer_matrix

for i in [10, 15, 16, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['pool_top%spercent' %i] = partial(temp_pool_top, p = i/100)

def temp_pool_top_control(community_function, p):
    """
    Pool random p% of communities. Round up to the nearest integer.
    """
    n_wells = len(community_function)
    randomized_community_function = community_function.copy()
    np.random.shuffle(randomized_community_function)
    sorted_community_function = np.sort(randomized_community_function)
    cut_off = sorted_community_function[int(np.floor(len(randomized_community_function)*(1-p)))]
    winner_index = np.where(randomized_community_function >= cut_off)[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
    return transfer_matrix

for i in [10, 15, 16, 20, 25, 28, 30, 33, 40, 50, 60]:
    globals()['pool_top%spercent_control' %i] = partial(temp_pool_top_control, p = i/100)


def Arora2019(community_function, n_rep=3):
    """
    Sub-divide wells of plate into multiple lines, each consisting of n_rep=3 communities
    Each round the highest function member of the line is used to colonize the next three triplicate wells of that line
    """
    n_wells = len(community_function)
    n_lines = int(np.ceil(n_wells/n_rep))
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix
    
    
def Arora2019_control(community_function, n_rep=3):
    """
  	Same as Arora2019 except within a subline the community is selected at random
    """
    n_wells = len(community_function)
    n_lines = int(np.ceil(n_wells/n_rep))
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
    Sub-divide wells of plate into n_lines=3 sublines
    Each round the highest function member of the line is used to colonize the wells of that lineage
    """
    n_wells = len(community_function)
    n_rep = int(np.ceil(n_wells/n_lines))
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[i*n_rep:(i*n_rep)+n_rep, winner_index+i*n_rep] = 1
    return transfer_matrix


def Raynaud2019a_control(community_function, n_lines=3):
    """
	Same as Raynaud2019a except the lineage member is selected at random
    """
    n_wells = len(community_function)
    n_rep = int(np.ceil(n_wells/n_lines))
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
    Same as Raynaud2019a except the top from each lineage is pooled (pooling across sublines)
    """
    n_wells = len(community_function)
    n_rep = int(np.ceil(n_wells/n_lines))
    transfer_matrix = np.zeros((n_wells,n_wells))
    for i in range(n_lines):
	    sorted_community_function = np.sort(community_function[i*n_rep:(i*n_rep)+n_rep])
	    cut_off = np.max(sorted_community_function)
	    winner_index = np.where(community_function[i*n_rep:(i*n_rep)+n_rep] == cut_off)[0]
	    transfer_matrix[:, winner_index+i*n_rep] = 1
    return transfer_matrix


def Raynaud2019b_control(community_function, n_lines=3):
    """
	Same as Raynaud2019b except the lineage member is selected at random
    """
    n_wells = len(community_function)
    n_rep = int(np.ceil(n_wells/n_lines))
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


def select_top(community_function):
    """
    Select the top community 
    """
    n_wells = len(community_function)
    winner_index = np.where(community_function >= np.max(community_function))[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells)
    t_old = list(winner_index) * n_wells
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
  winner_index = np.where(community_function >= cut_off)[0][::-1] 
  transfer_matrix = np.zeros((n_wells,n_wells))
  t_new = range(n_wells)
  t_old = [list(winner_index)[0]] * int(0.6 * n_wells) + [list(winner_index)[1]] * int(0.5 * n_wells)
  for i in range(n_wells):
      transfer_matrix[t_new[i], t_old[i]] = 1
  return transfer_matrix


def Williams2007a(community_function):
    """
    Williams2007a
    Select the top community and impose an bottleneck
    """
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    winner_index = np.where(community_function == np.max(community_function))[0][::-1] 
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells)
    t_old = list(winner_index) * n_wells
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 10**(-4) 
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
    transfer_matrix[:, winner_index] = 10**(-4)
    return transfer_matrix


def pair_top(community_function):
    """
    Pair the top communities. Each pairwise combination has roughly two replicates
    """
    import itertools
    n_wells = len(community_function)
    cut_off_percent = (np.sqrt(n_wells))/n_wells
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-cut_off_percent)))]
    winner_index = np.where(community_function >= cut_off)[0]
    pairs_list = list(itertools.combinations(winner_index, 2))
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells)
    t_old = list(winner_index) + pairs_list * (int(np.round(1/cut_off_percent)) + 1)
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
    return transfer_matrix



