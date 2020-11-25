#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 27 2019
@author: changyuchang
"""
import numpy as np
import random
from community_selection.A_experiment_functions import *

def resource_perturb(plate, params_simulation, keep):
    """
    Perturb the communities by shifting the medium composition
    """
    #Remove new fresh media
    plate.R = plate.R - plate.R0
    old_R0 = plate.R0[plate.N.columns[keep]]
    #First construct olist of possible metabolite perturbations (depends on r_type, either list of tuples of index opr simple list of index)
    if params_simulation['r_type'] == 'add': #Remove from top and add to random
        metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y and x ==  old_R0.idxmax()]
    if params_simulation['r_type']  == 'remove': #Remove from random and add to bottom
        metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y and y == old_R0.idxmin() and old_R0[x]>0]
    if params_simulation['r_type'] == 'rescale_add' or params_simulation['r_type'] == 'old':  # add to random
        metabolite_choice = [x for x in old_R0.index]
    if params_simulation['r_type'] == 'rescale_remove':  #remove from random    
        metabolite_choice = [x for x in old_R0.index if old_R0[x] >0]
    else: #default_resource_swap
        metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y]
    
    # If f6_target_resource, avoid target_resource
    if "target_resource" in params_simulation["selected_function"]:
        target_resource = old_R0.index[params_simulation["target_resource"]]
        if params_simulation['r_type'] == 'add': #Remove from top and add to random
            metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y and x ==  old_R0.idxmax() and x != target_resource and y != target_resource]
        if params_simulation['r_type']  == 'remove': #Remove from random and add to bottom
            metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y and y == old_R0.idxmin() and old_R0[x]>0 and x != target_resource and y != target_resource]
        if params_simulation['r_type'] == 'rescale_add' or params_simulation['r_type'] == 'old':  # add to random
            metabolite_choice = [x for x in old_R0.index if x != target_resource and y != target_resource]
        if params_simulation['r_type'] == 'rescale_remove':  #remove from random    
            metabolite_choice = [x for x in old_R0.index if old_R0[x] >0 and x != target_resource and y != target_resource]
        else: #default_resource_swap
            metabolite_choice = [(x,y) for x in old_R0.index for y in old_R0.index if x !=y and x != target_resource and y != target_resource]

    #next randomly pick element in list and apply pertubation 
    for k in plate.R0.columns:
        if k != plate.R0.columns[keep]:
            #So first default to kept media
            plate.R0[k] = old_R0
            if len(metabolite_choice) ==0: #If all possible pertubations have been carried out skip
                continue
            #Pick random pertubation
            r_id = random.choice(metabolite_choice)
            #perform pertubations
            if params_simulation['r_type']  == 'rescale_add': 
                plate.R0[k][r_id] = plate.R0[k][r_id]*(1+params_simulation['r_percent'])
            elif params_simulation['r_type'] == 'rescale_remove':
                plate.R0[k][r_id] = plate.R0[k][r_id]*(1-params_simulation['r_percent']) 
            elif params_simulation['r_type'] == 'old':
                plate.R0[k] = plate.R0[k] * (1-params_simulation['R_percent']) #Dilute old resource
                plate.R0[k][r_id] = plate.R0[k][r_id] + (params_simulation['R0_food']*params_simulation['R_percent']) #Add fixed percent
            else:
                plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*params_simulation['r_percent']) #add new resources
                plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-params_simulation['r_percent']) #remove new resources
            # Remove chosen pertubation as option for subsequent loop
            metabolite_choice = [x for x in metabolite_choice if x != r_id]
    plate.R0 = plate.R0/np.sum(plate.R0)*params_simulation['R0_food'] #Keep this to avoid floating point error and rescale when neeeded.
    #add new fresh environment (so that this round uses R0
    plate.R = plate.R + plate.R0
    return plate
                

def perturb(plate, params_simulation, keep):
    """
    Perturbs all communities except for the one specified by the argument keep. Default is the first well so keep = 0
    Only runs if directed selection is true
    """
    #Bottleneck
    if params_simulation['bottleneck']:
        dilution_matrix = np.eye(params_simulation['n_wells'])*params_simulation['bottleneck_size'] 
        dilution_matrix[keep,keep] = 1
        old_R = plate.R.copy()
        plate.Passage(dilution_matrix)
        plate.R = old_R.copy()  #knock_in isolates absent from all communities
    if params_simulation['knock_in']:
        knock_in_list = np.where(np.logical_and(np.array(np.sum(plate.N,axis=1)==0.0), plate.knock_in_species_function >= np.percentile(plate.knock_in_species_function, q = 100*params_simulation['knock_in_threshold'])))[0]
        # If f5, avoid using invader
        if "invader" in params_simulation["selected_function"]:
            knock_in_list[params_simulation["invader_index"]] = False
        for k in plate.N.columns:
            if k == plate.N.columns[keep] or len(knock_in_list) ==0.0:
                continue
            else:
                s_id = np.random.choice(knock_in_list) 
                plate.N[k][s_id]= 1/params_simulation["dilution"] * 1/params_simulation["scale"] #Knock in enough to survive 1 dilution even with no growth
                knock_in_list = knock_in_list[knock_in_list != s_id] 
    #knock_out isolates present in all communities
    if params_simulation['knock_out']:
        knock_out_list = np.where(np.sum(plate.N>0.0,axis=1) == params_simulation['n_wells'])[0]
        for k in plate.N.columns:
            if k == plate.N.columns[keep] or len(knock_out_list) ==0.0:
                continue
            else:
                s_id = np.random.choice(knock_out_list) 
                plate.N[k][s_id]= 0
                knock_out_list = knock_out_list[knock_out_list != s_id] 
    #Migrate taxa into the best performing community. By default migrations are done using power law model but can tune the diversity of migration using s_migration
    if params_simulation['migration']:
        migration_factor = np.ones(params_simulation['n_wells'])
        migration_factor[keep] = 0
        if np.isfinite(params_simulation['s_migration']):
            plate.N = migrate_from_pool(plate,migration_factor,params_simulation,power_law=False,n=params_simulation['n_migration'])
        else:
            plate.N = migrate_from_pool(plate,migration_factor,params_simulation,power_law = True,n=params_simulation['n_migration'])
        # If f5, avoid using invader
        if "invader" in params_simulation["selected_function"]:
            plate.N.iloc[params_simulation["invader_index"],] = 0
    #Migrate taxa into the best performing community. By default migrations are done using power law model but can tune the diversity of migration using s_migration
    if params_simulation['coalescence']:
        plate.Propagate(params_simulation["n_propagation"])
        plate.N = plate.N*(1-params_simulation['frac_coalescence']) + plate.prior_N*params_simulation['frac_coalescence']
        plate.R = plate.R*(1-params_simulation['frac_coalescence']) + plate.prior_R*params_simulation['frac_coalescence']
        plate.Passage(np.eye(params_simulation['n_wells'])*params_simulation['dilution'] )
    #Shift_R0
    if params_simulation['resource_shift']:
        plate = resource_perturb(plate, params_simulation, keep)
    return plate



# Design migration_factor (a sequence of binary factors)
def no_migration(community_function):
    """
    No migration
    """
    # Number of wells
    n_wells = len(community_function)

    # No migration
    migration_factor = np.zeros(n_wells)

    return migration_factor

def parent_migration(community_function):
    """
    Parent migration, migrate into all wells
    """
    # Number of wells
    n_wells = len(community_function)

    # All migration
    migration_factor = np.ones(n_wells)

    #dont migrate into winner
    winner_index = np.where(community_function >= np.max(community_function))[0][::-1] # Reverse the list so the higher
    migration_factor[winner_index] = 0
    return migration_factor

def directed_selection_migrate(community_function):
    """
    Sample new communities from species pool, coalesce the migrant communities to the species pools
    """
    # Number of wells
    n_wells = len(community_function)

    # Compute the cutoff based on the number of wells
    cut_off_percent = (np.sqrt(n_wells))/n_wells

    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)

    # Community function value cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-cut_off_percent)))]

    # Winner wells
    winner_index = np.where(community_function >= cut_off)[0][::-1]

    # Migration factor. A list of whether to migrate the community or not
    migration_factor = np.ones(n_wells) # Migrate all the wells except for the new wells that contain the winner replicate
    migration_factor[range(len(winner_index))] = 0 # Don't migrate to the winner wells

    return migration_factor

def migrate_half(community_function):
    # Number of wells
    n_wells = len(community_function)

    # Migration
    migration_factor = [1, 0] * int(n_wells/2)

    return migration_factor


def migrate_random(community_function):
    # Number of wells
    n_wells = len(community_function)

    # Migration
    migration_factor = np.random.binomial(1, 0.5, size = n_wells)

    return migration_factor
    
    

