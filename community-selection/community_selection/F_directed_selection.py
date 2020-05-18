#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 14 2020
@author: changyuchang
"""

# Replace all assumptions by params_simulation
# replace all hard coded parameters by params_simulation

  
def directed_selection(params = params, 
    params_simulation = params_simulation, 
    params_algorithm = params_algorithm,
    plate = plate
): 
  """
  Take the community  
  """

# Perturbation section
if 'bottleneck' in params_algorithm["algorithm_name"][0]  and selection_algorithm == 'select_top':
    winning_index = np.where(community_function >= np.max(community_function))[0][0] 
    dilution_matrix = np.eye(assumptions['n_wells'])
    dilution_matrix[winning_index,winning_index] = 1/params_simulation['bottleneck']
    dilution_matrix[winning_index,winning_index] = 1/params_simulation['bottleneck'] #Jean
    plate.Passage(dilution_matrix* params_simulation['bottleneck'])
if  'knock_in' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
    selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice
    winning_index = np.where(community_function >= np.max(community_function))[0][0] 
    for k in plate.N.columns:
        if k != plate.N.columns[winning_index]:
            sel = [x for x in np.where(plate.N[k]==0)[0] if x not in selected] #List to be selected from
            if 'isolates' in params_algorithm["algorithm_name"][0]:
                # Knock in an isolate that is 1) not in the resident community 2) have high monoculture function (above 90% isolates in the pool)  
                temp = np.logical_and(np.array(plate.N[k]==0), plate.isolate_function >= np.percentile(plate.isolate_function, q = 90))
                sel = [x for x in np.where(temp)[0] if x not in selected] #List to be selected from                
            if len(sel) == 0:
                continue
            s_id = np.random.choice(sel)
            selected.append(s_id)
            plate.N[k][s_id]= 1/params_simulation["dilution"] * 1/assumptions["scale"]
if 'knock_out' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
    selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice=
    winning_index = np.where(community_function >= np.max(community_function))[0][0]
    for k in plate.N.columns:
        if k != plate.N.columns[winning_index]:
            sel = [x for x in np.where(plate.N[k]>0)[0] if x not in selected] #List to be selected from
            if len(sel) == 0: #If you've already knocked out every species skip
                continue
            s_id = np.random.choice(sel)
            selected.append(s_id)
            plate.N[k][s_id]=0      
if 'resource' in params_algorithm["algorithm_name"][0] and selection_algorithm == 'select_top':
    selected = [] #Tracks pertubation id's to make sure we are not repeating the same one twice
    winning_index = np.where(community_function >= np.max(community_function))[0][0]
    #Remove fresh environment that was added by passage
    plate.R = plate.R - plate.R0
    #change default fresh renvironment so that all subsequent rounds use R0
    for k in plate.R0.columns:
        if k != plate.R0.columns[winning_index]: 
            #By default pick 2 resources at randomk
            sel = [(x,y) for x in np.where(plate.R0[k]>=0)[0] for y in np.where(plate.R0[k]>=0)[0] if x !=y]
            sel = [x for x in sel if x not in selected]
            if 'add' in params_algorithm["algorithm_name"][0]: #remove from top and add to random
                top = np.where(plate.R0[k]==np.max(plate.R0[k]))[0]
                sel = [x for x in sel if x[0] ==top]
            if 'remove' in params_algorithm["algorithm_name"][0]: #remove from random and add to bottom
                bottom = np.where(plate.R0[k]==np.min(plate.R0[k]))[0]
                sel = [x for x in sel if x[1] ==bottom]
            if len(sel) == 0: #If you've done every possible resource pertubation skip.
                continue
            r_id = sel[np.random.choice(len(sel))]
            selected.append(r_id)   
            if 'rescale_add' in params_algorithm["algorithm_name"][0]:  # Increase Fraction of resource
                plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]]*(1+params_simulation['R_percent']) #increase resource conc by fixed %
            elif 'rescale_remove' in params_algorithm["algorithm_name"][0]: # Decrease Fraction of resource
                plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-params_simulation['R_percent'])  #decrease resource 
            elif 'resource_old' in params_algorithm["algorithm_name"][0]:
                plate.R0[k] = plate.R0[k] * (1-params_simulation['R_percent']) #Dilute old resources
                plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (assumptions['R0_food']*params_simulation['R_percent'])
            else: #Default to resource swap.
                # For resource_add with R_percent 0.1, 0.5, and 1
                if "resource_add_10" in params_algorithm["algorithm_name"][0]:
                    plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*0.1) #add new resources
                    plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-0.1) #remove new resources
                elif "resource_add_2" in params_algorithm["algorithm_name"][0]:
                    plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*0.5) #add new resources
                    plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-0.5) #remove new resources
                elif "resource_add_1" in params_algorithm["algorithm_name"][0]:
                    plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*1) #add new resources
                    plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-1) #remove new resources
                # Else default remain the same
                else:
                    plate.R0[k][r_id[0]] = plate.R0[k][r_id[0]] + (plate.R0[k][r_id[1]]*params_simulation['R_percent']) #add new resources
                    plate.R0[k][r_id[1]] = plate.R0[k][r_id[1]]*(1-params_simulation['R_percent']) #remove new resources

    plate.R0 = plate.R0/np.sum(plate.R0)*assumptions['R0_food'] #Keep this to avoid floating point error and rescale when neeeded.
    #add new fresh environment (so that this round uses R0
    plate.R = plate.R + plate.R0
