#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp
import pandas as pd
import random

from community_simulator import *
from community_simulator.usertools import *
from community_selection.B_community_phenotypes import *
from community_selection.C_selection_algorithms import *
from community_selection.D_perturbation_algorithms import *

def MakeMatrices(assumptions):
    """
    Inherited function from community-simulator package
    
    Changes:
    
    - Add BINARY_GAMMA SAMPLING
    """
    #PREPARE VARIABLES
    #Force number of species to be an array:
    if isinstance(assumptions['MA'],numbers.Number):
        assumptions['MA'] = [assumptions['MA']]
    if isinstance(assumptions['SA'],numbers.Number):
        assumptions['SA'] = [assumptions['SA']]
    #Force numbers of species to be integers:
    assumptions['MA'] = np.asarray(assumptions['MA'],dtype=int)
    assumptions['SA'] = np.asarray(assumptions['SA'],dtype=int)
    assumptions['Sgen'] = int(assumptions['Sgen'])
    #Default waste type is last type in list:
    if 'waste_type' not in assumptions.keys():
        assumptions['waste_type']=len(assumptions['MA'])-1

    #Extract total numbers of resources, consumers, resource types, and consumer families:
    M = np.sum(assumptions['MA'])
    T = len(assumptions['MA'])
    S = np.sum(assumptions['SA'])+assumptions['Sgen']
    F = len(assumptions['SA'])
    M_waste = assumptions['MA'][assumptions['waste_type']]
    #Construct lists of names of resources, consumers, resource types, and consumer families:
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    waste_name = type_names[assumptions['waste_type']]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                      resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    
    #PERFORM GAUSSIAN SAMPLING
    if assumptions['sampling'] == 'Gaussian':
        #Initialize dataframe:
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add Gaussian-sampled values, biasing consumption of each family towards its preferred resource:
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                c.loc['F'+str(k)]['T'+str(j)] = c_mean + np.random.randn(assumptions['SA'][k],assumptions['MA'][j])*np.sqrt(c_var)
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c_var = assumptions['sigc']**2/M
            c.loc['GEN'] = c_mean + np.random.randn(assumptions['Sgen'],M)*np.sqrt(c_var)
                    
    #PERFORM BINARY SAMPLING
    elif assumptions['sampling'] == 'Binary':
        assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
        #Construct uniform matrix at total background consumption rate c0:
        c = pd.DataFrame(np.ones((S,M))*assumptions['c0']/M,columns=resource_index,index=consumer_index)
        #Sample binary random matrix blocks for each pair of family/resource type:
        for k in range(F):
            for j in range(T):
                if k==j:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1-assumptions['q'])
                    
                c.loc['F'+str(k)]['T'+str(j)] = (c.loc['F'+str(k)]['T'+str(j)].values 
                                                + assumptions['c1']*BinaryRandomMatrix(assumptions['SA'][k],assumptions['MA'][j],p))
        #Sample uniform binary random matrix for generalists:
        if 'GEN' in c.index:
            p = assumptions['muc']/(M*assumptions['c1'])
            c.loc['GEN'] = c.loc['GEN'].values + assumptions['c1']*BinaryRandomMatrix(assumptions['Sgen'],M,p)

    elif assumptions['sampling'] == 'Gamma':
        #Initialize dataframe
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add Gamma-sampled values, biasing consumption of each family towards its preferred resource
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                    thetac = c_var/c_mean
                    kc = c_mean**2/c_var
                    c.loc['F'+str(k)]['T'+str(j)] = np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c_var = assumptions['sigc']**2/M
            thetac = c_var/c_mean
            kc = c_mean**2/c_var
            c.loc['GEN'] = np.random.gamma(kc,scale=thetac,size=(assumptions['Sgen'],M))
    
    #PERFORM UNIFORM SAMPLING
    elif assumptions['sampling'] == 'Uniform':
        #Initialize dataframe:
        c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
        #Add uniformly sampled values, biasing consumption of each family towards its preferred resource:
        for k in range(F):
            for j in range(T):
                if k==j:
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                c.loc['F'+str(k)]['T'+str(j)] = c_mean + (np.random.rand(assumptions['SA'][k],assumptions['MA'][j])-0.5)*assumptions['b']
        if 'GEN' in c.index:
            c_mean = assumptions['muc']/M
            c.loc['GEN'] = c_mean + (np.random.rand(assumptions['Sgen'],M)-0.5)*assumptions['b']
    
    #PERFORM BINARY_GAMMA SAMPLING
    elif assumptions['sampling'] == 'Binary_Gamma':
        assert assumptions['muc'] < M*assumptions['c1'], 'muc not attainable with given M and c1.'
        #Construct uniform matrix at total background consumption rate c0:
        c = pd.DataFrame(np.ones((S,M))*assumptions['c0']/M,columns=resource_index,index=consumer_index)
        #Sample binary random matrix blocks for each pair of family/resource type:
        for k in range(F):
            for j in range(T):
                if k==j:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_mean = (assumptions['muc']/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                    c_var = (assumptions['sigc']**2/M)*(1+assumptions['q']*(M-assumptions['MA'][j])/assumptions['MA'][j])
                else:
                    p = (assumptions['muc']/(M*assumptions['c1']))*(1-assumptions['q'])
                    c_mean = (assumptions['muc']/M)*(1-assumptions['q'])
                    c_var = (assumptions['sigc']**2/M)*(1-assumptions['q'])
                c_mean_binary = assumptions['c0']+ assumptions['c1']*p
                c_var_binary = assumptions['c1']**2 *p*(1-p)
                c_mean_gamma = c_mean/c_mean_binary
                c_var_gamma = (c_var - c_var_binary*(c_mean_gamma**2))/(c_var_binary + c_mean_binary**2)
                thetac = c_var_gamma/c_mean_gamma
                kc = c_mean_gamma**2/c_var_gamma
                c.loc['F'+str(k)]['T'+str(j)] = (c.loc['F'+str(k)]['T'+str(j)].values + assumptions['c1']*BinaryRandomMatrix(assumptions['SA'][k],assumptions['MA'][j],p))*np.random.gamma(kc,scale=thetac,size=(assumptions['SA'][k],assumptions['MA'][j]))
        #Sample uniform binary random matrix for generalists:
        if 'GEN' in c.index:
            p = assumptions['muc']/(M*assumptions['c1'])
            c_mean = assumptions['muc']/M
            c_var = assumptions['sigc']**2/M
            c_mean_binary = assumptions['c0']+ assumptions['c1']*p
            c_var_binary = assumptions['c1']**2 *p*(1-p)
            c_mean_gamma = c_mean/c_mean_binary
            c_var_gamma = (c_var - c_var_binary*(c_mean_gamma**2))/(c_var_binary + c_mean_binary**2)
            thetac = c_var_gamma/c_mean_gamma
            kc = c_mean_gamma**2/c_var_gamma
            c.loc['GEN'] = (c.loc['GEN'].values + assumptions['c1']*BinaryRandomMatrix(assumptions['Sgen'],M,p))*np.random.gamma(kc,scale=thetac,size=(assumptions['Sgen'],M))
    else:
        print('Invalid distribution choice. Valid choices are kind=Gaussian and kind=Binary.')
        return 'Error'

    #SAMPLE METABOLIC MATRIX FROM DIRICHLET DISTRIBUTION
    DT = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
    for type_name in type_names:
        MA = len(DT.loc[type_name])
        if type_name is not waste_name:
            #Set background secretion levels
            p = pd.Series(np.ones(M)*(1-assumptions['fs']-assumptions['fw'])/(M-MA-M_waste),index = DT.keys())
            #Set self-secretion level
            p.loc[type_name] = assumptions['fs']/MA
            #Set waste secretion level
            p.loc[waste_name] = assumptions['fw']/M_waste
            #Sample from dirichlet
            DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=MA)
        else:
            if M > MA:
                #Set background secretion levels
                p = pd.Series(np.ones(M)*(1-assumptions['fw']-assumptions['fs'])/(M-MA),index = DT.keys())
                #Set self-secretion level
                p.loc[type_name] = (assumptions['fw']+assumptions['fs'])/MA
            else:
                p = pd.Series(np.ones(M)/M,index = DT.keys())
            #Sample from dirichlet
            DT.loc[type_name] = dirichlet(p/assumptions['sparsity'],size=MA)
        
    return c, DT.T

def reshape_plate_data(plate, params_simulation,transfer_loop_index):
    """
    Reshape the plate resource and consumer matrices (wider form) into a melted data.frame (longer form)
    """
    # Temporary function for adding variables to and melting df
    def melt_df(plate_df, data_type = "consumer"):
        # Consumers
        temp_df = pd.DataFrame(plate_df)
        total_number = temp_df.shape[0]
        
        ## Add variables
        temp_df["Type"] = np.repeat(data_type, total_number)
        temp_df["ID"] = range(total_number)
        temp_df["Transfer"] = np.repeat(str(transfer_loop_index), total_number)
        temp_df["exp_id"] = np.repeat(params_simulation['exp_id'] , total_number)

        ## Melt the df
        temp_df = pd.melt(temp_df, id_vars = ["exp_id","Transfer", "Type", "ID"], var_name = "Well", value_name = "Abundance")
        temp_df = temp_df[temp_df.Abundance != 0] # Remove zero abundances
        return temp_df
        
    # Melt the df
    temp_plate = plate.copy() # Copy the original plate 
    df_N = melt_df(temp_plate.N, data_type = "consumer")
    df_R = melt_df(temp_plate.R, data_type = "resource")
    df_R0 = melt_df(temp_plate.R0,data_type = "R0")
    
    # Concatenate dataframes
    merged_df = pd.concat([df_N, df_R,df_R0]) 
    merged_df["Index"] = list(range(0, merged_df.shape[0]))
    merged_df.set_index("Index", inplace = True)

    return merged_df # Return concatenated dataframe


def reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index):
    """
    Reshape the community function, richness, biomass into a melted data.frame
    """
    temp_vector1 = community_function.copy()
    temp_vector2 = richness.copy()
    temp_vector3 = biomass.copy()
    
    # Number of wells
    number_well = len(richness)

    # Make data.frame
    temp_df = pd.DataFrame({
        "exp_id": np.repeat(params_simulation['exp_id'], number_well),
        "Well": ["W" + str(i) for i in range(number_well)], 
        "Transfer": np.repeat(str(transfer_loop_index), number_well), 
        "CommunityPhenotype": temp_vector1,
        "Richness": temp_vector2,
        "Biomass": temp_vector3})
    
    # Turn the transfer columns as numeric
    temp_df[["Transfer"]] = temp_df[["Transfer"]].apply(pd.to_numeric)
    
    return temp_df 
    

def migrate_from_pool(plate,migration_factor,params_simulation, power_law = True, n = None):
    """
    Migrate from species pool to the plate mainly for directed selection)
    If power_law pool is true than sample n cells from species pool following power law distribution (default is same as inoculum)
    If power_law is false sample s_migration species from isolates with each total number of cells equivalent to n
    """
    from community_selection.usertools import sample_from_pool
    if n is None:
        n = params_simulation['n_inoc']
    if power_law:
        if np.sum(migration_factor) != 0:
            temp_params_simulation = params_simulation.copy() 
            migration_plate = sample_from_pool(plate.N, params_simulation,n=n) * migration_factor # Migration factor is a list determined by migration algorithms and community function
            plate_migrated = plate.N + migration_plate 
        else:
            plate_migrated = plate.N
    else: 
        if np.sum(migration_factor) != 0:
            migration_plate = plate.N.copy()
            migration_plate[:]  = 0
            for k in plate.N.columns:
                if migration_factor[np.where(plate.N.columns == k)[0]]>0: 
                    for j in range(0,params_simulation['s_migration']):
                        s_id = np.random.choice(np.where(plate.N[k]==0)[0])
                        migration_plate[k][s_id]= n * 1/params_simulation["scale"] * 1/params_simulation['s_migration']
            plate_migrated = plate.N + migration_plate
        else:
            plate_migrated = plate.N
    return plate_migrated


def passage_monoculture(plate_mono, f, scale = None, refresh_resource=True):
    """
    Reduced version of Passage(), for passaging a large set of wells without multinomial sampling
    Most code adapted from community-simulator
    """
    self = plate_mono.copy()
    #HOUSEKEEPING
    if scale == None:
        scale = self.scale #Use scale from initialization by default
    self.N[self.N<0] = 0 #Remove any negative values that may have crept in
    self.R[self.R<0] = 0
    
    #DEFINE NEW VARIABLES
    N_tot = np.sum(self.N)
    R_tot = np.sum(self.R)
    N = np.zeros(np.shape(self.N))
    
    #Poisson sample cells
    self.N = self.N * f *scale
    self.N.applymap(np.random.poisson)   
    self.N = self.N/scale

    if refresh_resource:
        self.R = self.R * f
        self.R = self.R+self.R0
        
    #In continuous culture, it is useful to eliminate the resources that are
    #going extinct, to avoid numerical instability
    else:
        R_tot = np.sum(self.R)
        R = np.zeros(np.shape(self.R))
        for k in range(self.n_wells):
            if f[k,k] > 0 and R_tot[k] > 0:
                R[:,k] += np.random.multinomial(int(scale*R_tot[k]*f[k,k]),(self.R/R_tot).values[:,k])*1./scale
        self.R = pd.DataFrame(R, index = self.R.index, columns = self.R.keys())

    return self


def simulate_community(params, params_simulation, params_algorithm, plate):
    """
    Simulate community dynamics by given experimental regimes
    
    params = parameter passed from community-simulator
    params_simulation = dictionary of parameters for running experiment
    params_algorithm = dictionary of algorithms that determine the selection regime, migration regime, and community pheotypes
    plate = Plate object specified by community-simulator
    
    Return:
    community_composition = concatenated, melted panda dataframe of community and resource composition in each transfer
    community_function = melted panda dataframe of community function
    """

    # Test the community function
    try:
        community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
    except:
        print('\n Community phenotype test failed')
        raise SystemExit

    # Save the inocula composition
    if params_simulation['save_composition']:
        plate_data_list = list() # Plate composition
        plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=0)  # Initial state
        plate_data_list.append(plate_data)
        composition_filename = params_simulation['output_dir'] + params_simulation['exp_id'] + '_composition.txt'   
        
    # Save the initial community function + richness + biomass
    if params_simulation['save_function']:
        community_function_list = list() # Plate composition
        richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
        biomass = list(np.sum(plate.N, axis = 0)) # Biomass
        function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index =0)
        community_function_list.append(function_data)
        function_filename = params_simulation['output_dir'] + params_simulation['exp_id'] + '_function.txt'   


    print("\nStart propogation")
    # Run simulation
    for i in range(0, params_simulation["n_transfer"]):
        # Algorithms used in this transfer
        phenotype_algorithm = params_algorithm["community_phenotype"][i]
        selection_algorithm = params_algorithm["selection_algorithm"][i]
#        migration_algorithm = params_algorithm["migration_algorithm"][i]
        print("Transfer " + str(i+1))

        # Propagation
        plate.Propagate(params_simulation["n_propagation"])
    
        # Measure Community phenotype
        community_function = globals()[params_algorithm["community_phenotype"][0]](plate, params_simulation = params_simulation) # Community phenotype
        
        # Append the composition to a list
        if params_simulation['save_composition'] and ((i+1) % params_simulation['composition_lograte'] == 0):
            plate_data = reshape_plate_data(plate, params_simulation,transfer_loop_index=i+1)  # Initial state
            plate_data_list.append(plate_data)

        if params_simulation['save_function'] and ((i+1) % params_simulation['function_lograte'] == 0):
            richness = np.sum(plate.N >= 1/params_simulation["scale"], axis = 0) # Richness
            biomass = list(np.sum(plate.N, axis = 0)) # Biomass
            function_data = reshape_function_data(params_simulation,community_function, richness, biomass, transfer_loop_index =i+1)
            community_function_list.append(function_data)

		#Store prior state before passaging (For coalescence)
        setattr(plate, "prior_N", plate.N)
        setattr(plate, "prior_R", plate.R)
        setattr(plate, "prior_R0", plate.R0)

        # Passage and tranfer matrix
        transfer_matrix = globals()[selection_algorithm](community_function)
        if params_simulation['monoculture']:
            plate = passage_monoculture(plate, params_simulation["dilution"])
        else:
            plate.Passage(transfer_matrix * params_simulation["dilution"])
        
        # Migration
        # m = globals()[migration_algorithm](community_function) 
        # plate.N = migrate_from_pool(plate, migration_factor = m, params_simulation = params_simulation) # By default, n_migration is the same as n_inoc
        
        # Perturbation
        if selection_algorithm == 'select_top' and params_simulation['directed_selection']:
        #if (i+1) % params_simulation['n_transfer_selection'] == 0 and params_simulation['directed_selection'] and params_simulation['n_transfer'] != (i+1):
            #keep = 0 # Always keep the first communtiy, by the design of selection transfer matrix, this is always the top of selected communities
            #plate = perturb(plate, params_simulation, keep = keep)
            plate  = perturb(plate,params_simulation,keep =  np.where(community_function >= np.max(community_function))[0][0])

    if params_simulation['save_composition']:
        pd.concat(plate_data_list).to_csv(composition_filename ,index=False)
    if params_simulation['save_function']:
        pd.concat(community_function_list).to_csv(function_filename ,index=False)
    print("\n"+ params_simulation["exp_id"]+ " finished")



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
				

def perturb(plate,params_simulation, keep):
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
		knock_in_list = np.where(np.logical_and(np.array(np.sum(plate.N,axis=1) ==0.0) , plate.isolate_function >= np.percentile(plate.isolate_function, q = 100*params_simulation['knock_in_threshold'])))[0]
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
	#Migrate taxa into the best performing community. By default migrations are done using power law model but can tune the diversity of migration using s_migration
	if params_simulation['coalescence']:
		plate.Propagate(params_simulation["n_propagation"])
		plate.N = plate.N*(1-params_simulation['frac_coalescence']) + plate.prior_N*params_simulation['frac_coalescence']
		plate.R = plate.R*(1-params_simulation['frac_coalescence']) + plate.prior_R*params_simulation['frac_coalescence']
		plate.Passage(np.eye(params_simulation['n_wells'])*params_simulation['dilution'] )
	#Shift_R0
	if params_simulation['resource_shift']:
		plate = resource_perturb(plate,params_simulation,keep)
	return plate


