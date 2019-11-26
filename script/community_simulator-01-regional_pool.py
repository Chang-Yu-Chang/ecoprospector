#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
# Create dynamcis
def dNdt(N,R,params):
    return MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return MakeResourceDynamics(assumptions)(N,R,params)
dynamics = [dNdt,dRdt]

# Create a regional species pool
def RegionalSpeciesPool(assumptions):
    # Parameters
    S_tot = int(np.sum(assumptions['SA'])+assumptions['Sgen']) # Total number of species (specialist + generalist)
    F = len(assumptions['SA']) # Number of consumer families
    consumer_names = ['S'+str(k) for k in range(S_tot)]
    family_names = ['F'+str(k) for k in range(F)]
    well_names = ['W'+str(k) for k in range(assumptions['n_wells'])]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
    +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    
    # Assign drawn values based on power-law distribution
    pool = np.random.power(1, size  = S_tot) 
    return pool/np.sum(pool)  # Relative species abundance in regional pool

# Sample communities from regional species pool
def SampleFromPool(plate_N, pool, scale=10**6, inocula=10**6):
    N0 = np.zeros((plate_N.shape))
    consumer_index = plate_N.index
    well_names = plate_N.columns
    
    for k in range(plate_N.shape[1]):
        consumer_list = np.random.choice(len(pool), size=len(pool), replace=True, p=pool)
        my_tab = pd.crosstab(index=consumer_list, columns="count") # Calculate the biomass count
        N0[my_tab.index.values,k] = np.ravel(my_tab.values / np.sum(my_tab.values) * inocula / scale) # Scale to sum

    # Make data.frame
    N0 = pd.DataFrame(N0,index=consumer_index,columns=well_names)
    return N0
