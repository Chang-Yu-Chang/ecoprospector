#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 27 2019
@author: changyuchang
"""

"""
Python functions for selection algorithms during plate passage in 96-well plates.
- It takes the community_function
- It takes the communtiy function to growht 

"""

import numpy as np
import scipy as sp

def no_selection(community_function):
    """
    Direct well-to-well transfer without selection
    """
    # Read number of wells 
    n_wells = len(community_function)
    return np.eye(n_wells)
    
def select_top25percent(community_function, p = 0.25):
    """
    Select the top 25% communities 
    """
    # Read number of wells 
    n_wells = len(community_function)
    
    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)
    
    # 25% cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    
    # Winner wells
    winner_index = np.where(community_function >= cut_off)[0][::-1] # Reverse the list so the higher 
    
    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/p)) # Old wells
        
    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
  
    return transfer_matrix


def select_top10percent(community_function, p = 0.1):
    """
    Select the top 10% communities 
    """
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

def select_bottom25percent(community_function):
    """
    Select the bottom 25% communities
    """
    # Read number of wells 
    n_wells = len(community_function)
    
    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)
    
    # 25% cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(0.25)))]
    
    # Winner wells
    winner_index = np.where(community_function <= cut_off)[0][::-1] # Reverse the list so the higher 
    
    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    t_new = range(n_wells) # New wells
    t_old = list(winner_index) * int(np.round(1/0.25)) # Old wells

    # Fill in the transfer matrix
    for i in range(n_wells):
        transfer_matrix[t_new[i], t_old[i]] = 1
  
    return transfer_matrix
  
def pool_top25percent(community_function, p = 0.25):
    """
    Select the top 25% communities, pool them all together, and replicate to all the new wells
    """    
    # Read number of wells 
    n_wells = len(community_function)
    
    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)
    
    # 25% cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    
    # Winner wells
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    
    # Transfer matrix
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
  
    return transfer_matrix 


def pool_top10percent(community_function, p = 0.1):
    """
    Select the top 10% communities, pool them all together, and replicate to all the new wells
    """    
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
    return transfer_matrix 

def pool_top20percent(community_function, p = 0.2):
    """
    Select the top 20% communities, pool them all together, and replicate to all the new wells
    """    
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
    return transfer_matrix 

def pool_top28percent(community_function, p = 0.28):
    """
    Select the top 28% communities, pool them all together, and replicate to all the new wells
    """    
    n_wells = len(community_function)
    sorted_community_function = np.sort(community_function)
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-p))) - 1]
    winner_index = np.where(community_function > cut_off)[0][::-1] # Reverse the list so the higher 
    transfer_matrix = np.zeros((n_wells,n_wells))
    transfer_matrix[:, winner_index] = 1
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






# # Chang-Yu Chang ----
# 
# def select_best_n(community_function, n_select=0.25):
#     # Read number of wells 
#     n = len(community_function)
#     
#     # Community function per transfer
#     sorted_community_function = np.sort(community_function)
#     
#     # 25% cutoff for selecting communities
#     cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
#     
#     # Winner wells
#     winner_index = np.where(community_function >= cut_off)
#     
#     # Empty transfer matrix
#     t = np.zeros((n, n))
#     t_x = range(n)
#     t_y = np.repeat(winner_index, int(np.round(1/n_select)))
#     t_y = t_y[:n]
#     
#     # Fill in the transfer matrix
#     for i in range(len(t_x)):
#         t[t_x[i], t_y[i]] = 1
#   
#     return t
# 
# def no_selection(community_function):
#     # Read number of wells 
#     n = len(community_function)
#     return np.eye(n)

# 
# # Jean Vila ----
# 
# def mixing_matrix(n_wells,n_mix):
#     '''
#     Generates a transfer matrix with the following  rule every community is transfered to itself and combined with 
#     n_mix other randomly chosen communities. Mixing is done in equal proportions
#     '''
#     TM = np.eye(n_wells)
#     for i in range(0,np.shape(TM)[0]):
#         r = np.random.choice(np.arange(0,np.shape(TM)[1]-1),n_mix,replace=False)
#         for j in r:
#             if j < i:
#                 TM[i,j] = 1
#             else:
#                 TM[i,j+1]  = 1
#     return(TM/(n_mix+1))
# 
# 
# 
# 
# # Xinwen Zhu ----
# 
# def pairwise_XZ(community_function, dilution=1/125):
    # import itertools
    # n = len(community_function)
    # # Community function per transfer
    # sorted_community_function = np.sort(community_function)
    # # cutoff top ten
    # top_10 = sorted_community_function[86]
    # winner_index = np.where(community_function >= top_10)
    # pairwise = list(itertools.combinations(winner_index[0], 2))
    # # Empty transfer matrix
    # t = np.zeros((n, n))
    # c=0
    # for i in range (0,90): #populate rows one by one
    #     t[i,pairwise[c][0]]=1
    #     t[i,pairwise[c][1]]=1
    #     c = c+1
    #     if c == 45:
    #         c = 0
    # c=0
    # for i in range (90,96):
    #     t[i,winner_index[0][c]]=1
    #     c = c+1
    # return t
#     
# def mate_top_XZ(community_function,dilution=1/125):
#     n = len(community_function)
#     sorted_community_function = np.sort(community_function)
#     top_index = np.where(community_function >= sorted_community_function[n-1])
#     t = np.eye(n)
#     t[:,top_index] =1
#     return t
# 
# 
# 

# 
# # Rachael Waymack ----
# 
# def select_above_avg_save_losers(community_function, n_select=0.25,dilution=1/125):
#     n=len(community_function)
#     # Community function per transfer
#     sorted_community_function = np.sort(community_function)
#     
#     #find the communities above the mean 
#     mean_community_function = np.mean(community_function)
#     top_communities_index = np.where(sorted_community_function > mean_community_function) #find the communities w fx above mean
#     avg_community = np.argmin(np.abs(sorted_community_function - mean_community_function)) #find index of value closest to the mean
#     #med_community_function = np.median(community_function)
#     #med_community = np.argmin(np.abs(sorted_community_function - med_community_function)) 
#     
#     #index_to_use= np.concatenate(top_communities_index,avg_community) #use those above mean, that closest to mean, and bottom 2 communities
#     worst_index=np.where(community_function == sorted_community_function[0])
#     worst_index=np.append(worst_index, np.where(community_function == sorted_community_function[1]))
#     #create array of indices to use - those above mean, that closest to mean, and bottom 2 communities
#     index_to_use=top_communities_index
#     index_to_use=np.append(top_communities_index, avg_community)
#     #index_to_use=np.append(top_communities_index, med_community)#try with median value
#     index_to_use=np.append(index_to_use,worst_index)
#      # Empty transfer matrix
#     t = np.zeros((n, n))
#     t_x = range(n) #0 to number of wells
#     t_y = np.repeat(index_to_use,2) #np.repeat(winner_index, int(np.round(1/n_select))) #repeat the winners for the fraction your selecting
#     t_y= t_y[:n]
#     #t_y = t_y[:n] #inc all values from beginning to n 
#     
#     # Fill in the transfer matrix
#     for i in range(len(t_x)):
#         #t[t_x[i], t_y[i]] = dilution
#         t[t_x[i],top_communities_index[-1]]=1  #always transfer top community
#         for p in range(len(t_y)): # don't 
#             t[t_x[i],t_y]=1
#    #     print(i)
#             
#     return t
# 
# # Paul Gerald Sanchez ----
# def MixFirstLast(community_function, n_select=0.25):
# 
#     # Read number of wells 
#     n = len(community_function)
#     
#     ascending_list = community_function.argsort() # sorts in ascending order
#     
#     s = int(n*n_select)
#     
#     first = ascending_list[-s:]
#     first = first[::-1] # first element is the community with highest function
#     last = ascending_list[:s] # list of least performing communities
#     
#     transfer_matrix = np.zeros([n,n])
#     c = 0
#     
#     for i in range(0,n):
#         
#         transfer_matrix[i,first[c]] = 1
#         
#         c = c + 1
#         
#         if c == s:
#             c = 0
#             
#     return transfer_matrix
#     
# def ChooseBest_MixFirstLast(community_function, n_select=0.25):
#     # Read number of wells 
#     n = len(community_function)
#     
#     ascending_list = community_function.argsort() # sorts in ascending order
#     
#     s = int(n*n_select)
#     
#     intermediate = ascending_list[-s:]
#     first = intermediate[::-1] # first element is the community with highest function
#     
#     transfer_matrix = np.zeros([n,n])
#     c = 0
#     
#     for i in range(0,n):
#         
#         transfer_matrix[i,first[c]] = 1
#         transfer_matrix[i,intermediate[c]] = 1
#         
#         c = c + 1
#         
#         if c == s:
#             c = 0
#             
#     return transfer_matrix
#     
# def ChooseRandom_NoMix(community_function, n_select=0.25):
#     
#     import random
#     
#     n = len(community_function)
#     s = int(n*n_select)
# 
#     well_index = range(n)
#     random_samples = random.sample(well_index,s)
# 
#     transfer_matrix = np.zeros([n,n])
#     c = 0
#     
#     for i in range(0,n):
#         
#         transfer_matrix[i,random_samples[c]] = 1
#         
#         c = c + 1
#         
#         if c == s:
#             c = 0
#             
#     return transfer_matrix
# 
# def ChooseRandom_Mix(community_function, n_select=0.25):
#     
#     import random
#     
#     n = len(community_function)
#     s = int(n*n_select)
#     well_index = range(n)
#     random_samples1 = random.sample(well_index,s)
#     random_samples2 = random.sample(well_index,s)
# 
#     transfer_matrix = np.zeros([n,n])
#     c = 0
#     
#     for i in range(0,n):
#         
#         transfer_matrix[i,random_samples1[c]] = 1
#         transfer_matrix[i,random_samples2[c]] = 1
#         
#         c = c + 1
#         
#         if c == s:
#             c = 0
#             
#     return transfer_matrix
# 
# def ChooseBest_Shuffle_Mix(community_function, n_select=0.25):
# 
#     import random
#     from sklearn.utils import shuffle
#     
#     n = len(community_function)
#     s = int(n*n_select)
# 
#     ascending_list = community_function.argsort() # sorts in ascending order
#     
#     first = ascending_list[-s:]
#     
#     random_samples1 = shuffle(first,random_state = random.randint(0,s))
#     random_samples2 = shuffle(first,random_state = random.randint(0,s))
# 
#     transfer_matrix = np.zeros([n,n])
#     c = 0
#     
#     for i in range(0,n):
#         
#         transfer_matrix[i,random_samples1[c]] = 1
#         transfer_matrix[i,random_samples2[c]] = 1
#         
#         c = c + 1
#         
#         if c == s:
#             c = 0
#             
#     return transfer_matrix
# 
# 
# # Molly Bassette ----
# def select_mix(community_function):
#   
#     # Read number of wells 
#     n = len(community_function) # n = 96
#     
#     # sort into ascending order 
#     sorted_community_function = np.argsort(community_function)
#     
#     # Top 3
#     top_three = sorted_community_function[-3:]
#     
#     # Second 3
#     second_three = sorted_community_function[90:92]
#     
#     # Third 3
#     third_three = sorted_community_function[87:89]
#     
#     # Fourth 3
#     fourth_three = sorted_community_function[84:86]
#     
#     # Empty transfer matrix
#     t = np.zeros((n, n))
#     
#     t[0:31,top_three] = 1
#     t[0:63,second_three] = 1
#     t[32:63,third_three] = 1
#     t[64:95,top_three] = 1
#     t[64:95,fourth_three] = 1
#         
#     return t
# 
# 
# # Brian von Herzen ----
# def exponealing(community_function, propagation_time = 1, n_propagation=20, n_select=0.06):
#     from scipy.stats import expon
#     
#     # Read number of wells 
#     n = len(community_function)
#     ngens = n_propagation #number of generations
# 
#     # Community function per transfer
#     sorted_community_function = np.sort(community_function)
#     ranked_index = np.argsort(community_function)
#     #print(community_function[ranked_index])
#     #print(sorted_community_function)
#     
#     # cutoff for selecting communities
#     cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
#     
#     # Winner wells
#     winner_index = np.where(community_function >= cut_off)
#     
#     # xyzzy
#     if not hasattr(exponealing, 'xyzzy'): 
#         exponealing.xyzzy = -1 
#         
#     exponealing.xyzzy = propagation_time
#  #   print(exponealing.xyzzy)
#     
#     # Empty transfer matrix
#     t = np.zeros((n, n))
#     t_x = range(n)
#     #t_y = np.repeat(winner_index, int(np.round(1/n_select)))
#     #t_y = t_y[:n]
#     #print(t_y)
#     # Fill in the transfer matrix
#     for i in range(len(t_x)):
#         t[t_x[i], t_x[i]] = 1   
#       
#     # add in the random pairwise exponential tail annealing distribution 
# 
#     annealing_curve = np.zeros(ngens+1)
#     for i in range(ngens+1):
#         annealing_curve[i] = 1.5 + 28.5/((.95*i)+1) 
#     
#     data_expon = np.round(expon.rvs(scale=annealing_curve[exponealing.xyzzy],loc=0,size=90)) % 90
#     # print(data_expon)
#    
#     
#     for i in range(len(t_x)-6):
#         t[ranked_index[t_x[i]], ranked_index[95-int(data_expon[i])]] = 1  
#     if exponealing.xyzzy >= ngens-2:
#         t = np.zeros((n, n))
#         t_x = range(n)
#          # last rond Fill in the transfer matrix with best output
#         for i in range(len(t_x)):
#             t[t_x[i], ranked_index[95]] = 1   
#       
#     # plot_transfer_matrix(t)
#     return t
# 
# #Select for the communities pairwise with an exponential annealing transfer function.
# #keeps the top 6 unaltered but mixes the others with the top performers. -- Brian von Herzen
# def exponealing2(community_function, timestamp=17, n_select=0.06):
#     from scipy.stats import expon
#     
#     # Read number of wells
#     n = len(community_function)
#     ngens = 20 #number of generations
#     # Community function per transfer
#     sorted_community_function = np.sort(community_function)
#     ranked_index = np.argsort(community_function)
#     #print(community_function[ranked_index])
#     #print(sorted_community_function)
#     # cutoff for selecting communities
#     cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
#     # Winner wells
#     winner_index = np.where(community_function >= cut_off)
#     print(timestamp)
#     # Empty transfer matrix
#     t = np.zeros((n, n))
#     t_x = range(n)
#     #t_y = np.repeat(winner_index, int(np.round(1/n_select)))
#     #t_y = t_y[:n]
#     #print(t_y)
#     # Fill in the transfer matrix
#     for i in range(len(t_x)):
#         # t[t_x[i], t_x[i]] = 1
#         t[ranked_index[i], ranked_index[95-(i%48)]] = 1
#     # add in the random pairwise exponential tail annealing distribution
#     annealing_curve = np.zeros(ngens+1)
#     for i in range(ngens+1):
#         annealing_curve[i] = 1.5 + 28.5/((.95*i)+1)
#     data_expon = np.round(expon.rvs(scale=annealing_curve[timestamp],loc=0,size=90)) % 90
#     # print(data_expon)
#     for i in range(len(t_x)-6):
#         t[ranked_index[t_x[i]], ranked_index[95-int(data_expon[i])]] = 1
#     if timestamp >= ngens-2:
#         t = np.zeros((n, n))
#         t_x = range(n)
#         # last rond Fill in the transfer matrix with best output
#         for i in range(len(t_x)):
#             t[t_x[i], ranked_index[95]] = 1
#     #plot_transfer_matrix(t)
#     return t
# 
# 
# # Stefan ----
# def drunkards_walk(community_function, n_select=0.25):
#    np.random.seed(0)
#    community_function = additive_community_function(plate1)
#    n = len(community_function)
#    t_x = []
#    t_y = []
#    for i in range(n):
#        t_x.append(int(np.floor(np.random.uniform()*96)))
#        t_y.append(int(np.floor(np.random.uniform()*96)))
#    t = np.zeros((n, n))
#    # Fill in the transfer matrix
#    for i in range(len(t_x)):
#        t[t_x[i], t_y[i]] = 1
#    return t
# 
# # Julie ----
# def select_one_more(community_function, n_select=0.1, dilution=1/125):
# 
#     # Read number of wells 
#     n = len(community_function)
# 
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
# 
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
# 
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         sample = ascending_list[-i:]
# 
#         #loop through list of performers in sample
#         #take the sample that's the top performer, place it in the ith well of the transfer matrix
#         #set this equal to 1 since we are making the transfer matrix binary
#         for j in range(len(sample)):
# 
#             transfer_matrix[i, sample[j]] = 1
# 
#     return transfer_matrix
#     
# def select_one_less(community_function, n_select=0.1, dilution=1/125):
#     
#     # Read number of wells 
#     n = len(community_function)
#     
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
# 
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
#         
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         sample = ascending_list[-(i-1):]
#         
#         #loop through list of performers in sample
#         #take the j sample that's the top performer, place it in the ith well of the transfer matrix
#         #set this equal to 1 since we are making the transfer matrix binary
#         for j in range(len(sample)):
#             
#             transfer_matrix[i, sample[j]] = 1
#             
#     return transfer_matrix
#     
# def select_best_random(community_function, n_select=0.1, dilution=1/125):
#     
#     # Read number of wells 
#     n = len(community_function)
#     
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
# 
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
#         
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         best_sample = ascending_list[int(len(ascending_list)-1)]
#         random_sample = ascending_list[np.random.choice(ascending_list)]
#         sample = [best_sample, random_sample]
#         
#         for j in range(len(sample)):
#             
#             #loop through list of performers in sample
#             #take the j sample that's the top performer, place it in the ith well of the transfer matrix
#             #set this equal to 1 since we are making the transfer matrix binary
#             transfer_matrix[i, sample[j]] = 1
#             
#     return transfer_matrix
#     
# def select_best(community_function, n_select=0.1, dilution=1/125):
#     
#     # Read number of wells 
#     n = len(community_function)
#     
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
# 
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
#         
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         sample = ascending_list[n-1]
#         
#         transfer_matrix[i, sample] = 1
#             
#     return transfer_matrix
#     
# def select_worst(community_function, n_select=0.1, dilution=1/125):
#     
#     # Read number of wells 
#     n = len(community_function)
#     
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
#     #::-1 inverts list, so now it is descending
#     ascending_list = ascending_list[::-1] 
#     
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
#         
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         sample = ascending_list[n-1]
#         
#         #for j in range(len(sample)):
#             
#         #loop through list of performers in sample
#         #take the j sample that's the top performer, place it in the ith well of the transfer matrix
#         #set this equal to 1 since we are making the transfer matrix binary
#         transfer_matrix[i, sample] = 1
#             
#     return transfer_matrix
#     
# def select_one_more_worst(community_function, n_select=0.1, dilution=1/125):
# 
#     # Read number of wells 
#     n = len(community_function)
# 
#     #sort community_function in ascending (worst to best) order, return list of wells in this order
#     ascending_list = community_function.argsort() 
#     #::-1 inverts list, so now it is descending
#     ascending_list = ascending_list[::-1] 
# 
#     #make transfer matrix filled with zeros
#     transfer_matrix = np.zeros((n,n))
# 
#     for i in range(len(ascending_list)):
# 
#         #select i top performers (negative because it's ascending)
#         #returns a list of these top performers
#         sample = ascending_list[-i:]
# 
#         #loop through list of performers in sample
#         #take the sample that's the top performer, place it in the ith well of the transfer matrix
#         #set this equal to 1 since we are making the transfer matrix binary
#         for j in range(len(sample)):
# 
#             transfer_matrix[i, sample[j]] = 1
# 
#     return transfer_matrix
