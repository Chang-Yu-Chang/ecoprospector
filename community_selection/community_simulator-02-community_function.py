"""
# Assign additive neutral traits to each of species
def species_traits(assumptions, mu=0, sigma=0.1): # mean and standard deviation
    # Number of species in pools
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    # Assign additive traits to each species
    additive_traits = np.random.normal(mu, sigma, size=S_tot)
    
    return additive_traits
        

# Sum up the community function for a well
def additive_community_function(plate, sigma = 0.001): #Sigma is the measurement error
    N_tot = plate.N.shape[1]
    return(np.sum(plate.N.values*traits[:,None],axis=0)*(1+ np.random.normal(0,sigma,N_tot)))
"""
    
# Plot community function as a function of time
def plot_community_function(function_df):
    import matplotlib.pyplot as plt
    #function_df
    time = range(0, len(function_df))
    plt.plot(time,function_df, 'ko', markersize=2)
    ax = plt.gca()
    ax.set_xlabel("transfer")
    ax.set_ylabel("Community function")
    plt.show()


# Update the plate with only one single consumer and multiple resource types
def plate_single_consumer(plate, single_consumer_index = 1):
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen'])
    F = len(assumptions['SA'])
    
    #Construct lists of names of resources, consumers, resource types, consumer families and wells:
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S_tot)]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    well_names = ['W'+str(k) for k in range(assumptions['n_wells'])]
        
    N0 = np.zeros((S_tot,assumptions['n_wells']))

    for k in range(assumptions['n_wells']):
        N0[single_consumer_index,k] = 1. 

    N0 = pd.DataFrame(N0,index=consumer_index,columns=well_names)
    
    return N0
    

# Community function as farmer
def farmer_community_function(plate): #Sigma is the measurement error
    # Number of wells in the 
    n_wells = plate.N.shape[1]
    np.random.seed(1) 
    assumptions2=a_default.copy() #Start with default parameters
    assumptions2.update({'n_wells':n_wells, 'c1' :.01, 'muc':0.1, 'm':0, #switch off mortality to make things simpler
                        'Sgen': 0, 'S': 1}) 
                    
    S_tot = int(np.sum(assumptions2['SA']) + assumptions2['Sgen'])

    init_state = MakeInitialState(assumptions2)
    params = MakeParams(assumptions2)
    params['c'].iloc[:,0] = 0 # Force the probe not able to consumer R0

     
    # Make plate
    plate2 = Community(init_state,dynamics,params,scale = 10**6,parallel=True) # Reset the community

    ## Inoculate probe consumer 
    plate2.N = plate_single_consumer(plate2, single_consumer_index = 0) # ALways use the first consumer (F1S1) as the probe

    ## Supernatant
    plate2.R = plate1.R

    # Propagation for one passage
    for i in range(1):
        plate2.Propagate(48)

    return np.array(plate2.N.iloc[0,:])



# Test commnuity phenotype algorithms
def test_community_phenotype(phenotype_algorithm, assumptions):
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
    plate1.N = SampleFromPool(plate1.N,species_pool) # Populate the well by sampling from the species pool
    
    # Propagate
    plate1.Propagate(24)
    community_function = globals()[phenotype_algorithm](plate1) # Call function name from selection_function_list
    
    return len(community_function) == n







