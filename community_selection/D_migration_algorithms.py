#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""

"""
This python script contains the migration algorithms written by students
"""

import numpy as np
import scipy as sp


# Chang-Yu Chang ----

# Design migration_factor (a sequence of binary factors) 
def no_migration(community_function, migration_type = 'no_migration'):
    # Number of wells
    n_wells = len(community_function)
    
    # No migration 
    migration_factor = np.zeros(n_wells)

    return migration_factor
  
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
  
 
 
# # Xinwen Zhu ----
# 
# def migration_factor_for_pairwise_XZ(community_function):
#     # Number of wells
#     n_wells = len(community_function)
#     migration_factor=np.zeros(n_wells)
#     migration_factor[0:44] = 1
#     return migration_factor 
# 
# 
# 
# # Rachel Waymack ----
# 
# def migrate_by_fraction(community_function, migration_type = 'migrate_by_fraction', n_repeats=5):
#     # Number of wells
#     n_wells = len(community_function)
#     
#     # Migration
#     if migration_type == 'migrate_half': # Migrate half 
#         migration_factor = [1, 0] * int(n_wells/2)
#     elif migration_type == 'migrate_by_faction':
#         migration_factor = [1]
#         for i in range(n_repeats-1):
#             migration_factor = np.append(migration_factor,0)
#     elif migration_type == 'random': # Random migration. Each well has probability=0.5 being migrated
#         migration_factor = np.random.binomial(1, 0.5, size = n_wells)
#     else: # No migration 
#         migration_factor = np.zeros(n_wells)
#     return migration_factor
# 
# 
# def migration_factor_rand_sneak_in(community_function, migration_type = 'migrate_by_fraction', n_repeats=5):
#     # Number of wells
#     n_wells = len(community_function)
#         # Migration
#     migration_factor_index = np.random.choice(n_wells)
#     migration_factor=np.zeros(n_wells)
#     for i in range(np.random.choice(n_wells)):
#         migration_factor_index = np.random.choice(n_wells)
#         migration_factor[migration_factor_index] = 1
#     
#     return migration_factor
# 
# # Julie ----
# 
# def migration_into_winner(community_function, migration_type = 'migrate_into_winners'):
#     # Number of wells
#     n_wells = len(community_function)
#     
#     # Community function per transfer
#     sorted_community_function = np.sort(community_function)
#     
#     #select top 25%
#     n_select=0.25
#     
#     # 25% cutoff for selecting communities
#     cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
#     
#     # Winner wells
#     winner_index = np.where(community_function >= cut_off)
#     
#     # Migration
#     if migration_type == 'migrate_half': # Migrate half 
#         migration_factor = [1, 0] * int(n_wells/2)
#     elif migration_type == 'random': # Random migration. Each well has probability=0.5 being migrated
#         migration_factor = np.random.binomial(1, 0.5, size = n_wells)
#     elif migration_type == 'random_90': #Random migration with probability=0.9
#         migration_factor = np.random.binomial(1, 0.9, size = n_wells)
#     elif migration_type == 'super_random': #Random migration 
#         migration_factor = np.random.randint(2, size = n_wells)
#     elif migration_type == 'migrate_all': #Migrate only into low performing wells
#         migration_factor = [1] * int(n_wells)
#     elif migration_type == 'migrate_into_winners': #Migrate only into high performing wells
#         migration_factor = np.zeros(n_wells)
#         for i in range(len(winner_index)):
#             for j in range(len(migration_factor)):
#                 if np.any(winner_index[i]==j):
#                     migration_factor[j] = 1
#     elif migration_type == 'migrate_into_losers': #Migrate only into low performing wells
#         migration_factor = np.full(n_wells,1)
#         for i in range(len(winner_index)):
#             for j in range(len(migration_factor)):
#                 if np.any(winner_index[i]==j):
#                     migration_factor[j] = 0
#     else: # No migration 
#         migration_factor = np.zeros(n_wells)
#         
#     return migration_factor
# 
# 
# 
# # Molly Bassette ----
# 
# def migrate_first_half(community_function):
#     # Number of wells
#     n_wells = len(community_function)
#     
#     # Migration
#     migration_factor = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0] * int(n_wells/16)
# 
#     return migration_factor
#     

