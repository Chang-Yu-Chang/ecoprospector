#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
# Example selection function algorithms; select for the top 25% of measured community functions
def select_best_n(community_function, n_select=0.25):
    # Read number of wells 
    n = len(community_function)
    
    # Community function per transfer
    sorted_community_function = np.sort(community_function)
    
    # 25% cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-n_select))) - 1]
    
    # Winner wells
    winner_index = np.where(community_function >= cut_off)
    
    # Empty transfer matrix
    t = np.zeros((n, n))
    t_x = range(n)
    t_y = np.repeat(winner_index, int(np.round(1/n_select)))
    t_y = t_y[:n]
    
    # Fill in the transfer matrix
    for i in range(len(t_x)):
        t[t_x[i], t_y[i]] = 1
  
    return t


# Plot the transfer matrix
def plot_transfer_matrix(t):
    t    
    fig,ax=plt.subplots()
    sns.heatmap(t,ax=ax)
    ax.set_xlabel('Old well',fontsize=14)
    ax.set_ylabel('New well',fontsize=14)
    ax.set_title(r'Transfer Matrix $f$',fontsize=14)
    plt.show()



# Simulate with different community function measurement error
def simulate_community_function_error(measurement_error_sigma = 0.01):

    def additive_community_function(plate,traits,sigma = measurement_error_sigma): #Sigma is the measurement error
        N_tot = plate.N.shape[1]
        return(np.sum(plate.N.values*traits[:,None],axis=0)*(1+ np.random.normal(0,sigma,N_tot)))
    
    # Set up
    np.random.seed(0) #you are all going to get the same species pool
    n = 96 #No of wells
    assumptions=a_default.copy() #Start with default parameters
    assumptions.update({'n_wells':n, # To start only 1 community
                        'c1' :.01,
                        'muc':0.1, #This is just a rescaling of the uptake rate
                        'm':0}) #switch off mortality to make things simpler
    init_state = MakeInitialState(assumptions)
    params = MakeParams(assumptions)
    species_pool = RegionalSpeciesPool(assumptions) # Generate a species pool
    trait_map = species_traits(assumptions)

    # Make plate
    np.random.seed(0) # Change to a unique number so that you each start the experiment with a slightly different species pool
    plate1 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community
    plate1.N = SampleFromPool(plate1.N,species_pool)# Populate the well by sampling from the species pool


    #Propagate for 20 passage each lasting '48hr' and record community function every passage
    function_df = list()
    for i in range(0,10):
        plate1.Propagate(48)
        community_function = additive_community_function(plate1,trait_map)
        t = select_best_n(community_function, n_select=0.25)
        plate1.Passage(t*1/125)
        function_df.append(community_function)
        
    return plate1, function_df


# Test selection function algorithms
def test_selection_function(selection_algorithm, assumptions):
    n=96
    np.random.seed(0) #you are all going to get the same species pool
 #   assumptions=a_default.copy() #Start with default parameters
    assumptions.update({'n_wells':n, 'c1' :.01, 'muc':0.1, 'm':0}) 
    init_state = MakeInitialState(assumptions)
    params = MakeParams(assumptions)
    species_pool = RegionalSpeciesPool(assumptions) # Generate a species pool

    # Make plate
    np.random.seed(0) # Change to a unique number so that you each start the experiment with a slightly different species pool
    plate1 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community
    plate1.N = SampleFromPool(plate1.N,species_pool)# Populate the well by sampling from the species pool
    
    plate1.Propagate(24)
    community_function = additive_community_function(plate1)
    t = globals()[selection_algorithm](community_function) # Call function name from selection_function_list
    
    return any(pd.DataFrame(t) != 1)








