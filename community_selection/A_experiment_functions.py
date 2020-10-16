#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""
import numpy as np
from community_simulator import *
from community_simulator.usertools import *
from community_selection.__init__ import *

# Species features

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

def draw_species_function(assumptions):
	"""
	Draw species-specific functions
	
	assumptions = dictionary of metaparameters from community-simulator
	
	Return:
	function_species, function_interaction
	"""
	# Total number of species in the pool
	S_tot = int(np.sum(assumptions['SA']) + assumptions['Sgen']) 
	
	# Species-specific function, 1-D array
	f1_species_smooth = np.random.normal(0, assumptions["sigma_func"], size = S_tot)
	f1_species_rugged = np.random.normal(0, assumptions["sigma_func"], size = S_tot) * np.random.binomial(1, 0.2, size = S_tot)
	
	# Interaction-specific function, 2-D n by n array
	f2_species_smooth = np.random.normal(0, assumptions["sigma_func"] * assumptions["alpha_func"], size = S_tot * S_tot).reshape(S_tot, S_tot)
	f2_species_rugged = np.random.binomial(1, 0.2, S_tot**2).reshape(S_tot, S_tot) * np.array(np.random.normal(0, assumptions["sigma_func"] * assumptions["alpha_func"], size = S_tot * S_tot)).reshape(S_tot, S_tot)
	
	# Remove diagonals in the interation matrix
	np.fill_diagonal(f2_species_smooth, 0)
	np.fill_diagonal(f2_species_rugged, 0)

	return f1_species_smooth, f1_species_rugged, f2_species_smooth, f2_species_rugged

def draw_species_cost(per_capita_function, assumptions):
	"""
	Draw species-specific function cost
	k_i is a conversion factor that specifies cost per function 
	"""
	
	if assumptions["cost_mean"] !=0:
		cost_var = assumptions["cost_sd"]**2
		cost_k = assumptions["cost_mean"]**2/cost_var
		cost_theta = cost_var/assumptions["cost_mean"]
		cost = np.random.gamma(shape = cost_k, scale = cost_theta, size = len(per_capita_function))
		g0 = assumptions["g0"]
		gi = g0/(1-per_capita_function*cost)
	else: 
		gi = np.repeat(assumptions["g0"], len(per_capita_function))
	
	return gi

def add_community_function(plate, assumptions, params):
	"""
	Add the function attribute to the community
	
	For f1 and f3, add species_function 
	For f2 and f4, add interaction_function
	For f5, add invasion_plate_t0 and invasion_plate_t1
	For f6, f7, and f8, add resident_plate_t0_N, resident_plate_t1_N, resident_plate_t0_R, and resident_plate_t1_R
	
	if isolates calculate function for every isolate in monoculture.
	"""
	
	#Generate per capita species function
	np.random.seed(assumptions['seed']) 
	f1_species_smooth, f1_species_rugged, f2_species_smooth, f2_species_rugged = draw_species_function(assumptions)
	
	# Species function for f1 additive community function
	setattr(plate, "f1_species_smooth", f1_species_smooth)
	setattr(plate, "f1_species_rugged", f1_species_rugged)

	# Species interaction function for f2 Interactive function
	setattr(plate, "f2_species_smooth", f2_species_smooth)
	setattr(plate, "f2_species_rugged", f2_species_rugged)


	# Invasion function f5 or knock_in with a threshold requires us to grow isolates in monoculture to obtain their abundance.
	if (assumptions["selected_function"] == 'f5_invader_growth') | (assumptions['knock_in']):
		print("\nStabilizing monoculture plate")
		# Keep the initial plate R0 for function f7 
		setattr(plate, "R0_initial", plate.R0)
		
		assumptions_invasion = assumptions.copy()
		params_invasion = params.copy()
		
		#Update assumptions
		assumptions_invasion.update({"n_wells": np.sum(assumptions["SA"])  + assumptions["Sgen"]})
		assumptions_invasion.update({"monoculture":True})

		# Make plates
		plate_invasion = make_plate(assumptions_invasion,params_invasion)
		
		# Species function, f1 and f3 (to calculate function at end)
		setattr(plate_invasion, "species_function", function_species) # Species function for additive community function

		# Interactive functions, f2 , f2b and f4
		setattr(plate_invasion, "interaction_function",function_interaction) # Interactive function for interactive community function
		setattr(plate_invasion, "interaction_function_p25", function_interaction_p25)   
		
		
		# Grow the invader plate  to equilibrium
		for i in range(assumptions_invasion["n_transfer"] - assumptions_invasion["n_transfer_selection"]):
			plate_invasion.Propagate(assumptions_invasion["n_propagation"])
			plate_invasion = passage_monoculture(plate_invasion, assumptions_invasion["dilution"])
		
		#  1 final growth cycle before storing data
		plate_invasion.Propagate(assumptions_invasion["n_propagation"])
		
		# find well with highest biomass
		dominant_index = np.where(np.sum(plate_invasion.N, axis = 0) == np.max(np.sum(plate_invasion.N, axis = 0)))[0][0] # Find the well with the highest biomass

		# Duplicate the chosen community  to the entire plate and save this in a data.frame to be add to as an attribute of the plate
		invader_N = pd.DataFrame()
		invader_R = pd.DataFrame()
		invader_R0 = pd.DataFrame()

		for i in range(assumptions["n_wells"]):
			invader_N["W" + str(i)] = plate_invasion.N["W" + str(dominant_index)]
			invader_R["W" + str(i)] = plate_invasion.R["W" + str(dominant_index)]
			invader_R0["W" + str(i)] = plate_invasion.R0["W" + str(dominant_index)]

		#Add the invasion plate to the attr of community
		setattr(plate, "invader_N", invader_N)
		setattr(plate, "invader_R", invader_R)
		setattr(plate, "invader_R0", invader_R0)
		setattr(plate, "isolate_abundance", np.sum(plate_invasion.N,axis=1)) 
		setattr(plate, "isolate_function", globals()[assumptions["selected_function"]](plate_invasion, params_simulation = assumptions))     
	
		print("\nFinished Stabilizing monoculture plate")
	return plate

# Plate

def sample_from_pool(plate_N, assumptions, n = None):
	"""
	Sample communities from regional species pool.

	plate_N = consumer data.frame
	"""
	S_tot = plate_N.shape[0] # Total number of species in the pool
	N0 = np.zeros((plate_N.shape)) # Make empty plate
	consumer_index = plate_N.index
	well_names = plate_N.columns
	if n is None:
		n = assumptions['n_inoc'] #if not specified n is n_inoc
	# Draw community
	if assumptions['monoculture'] == False:
		# Sample initial community for each well
		for k in range(plate_N.shape[1]):
			pool = np.random.power(0.01, size = S_tot) # Power-law distribution
			pool = pool/np.sum(pool) # Normalize the pool
			consumer_list = np.random.choice(S_tot, size = n , replace = True, p = pool) # Draw from the pool
			my_tab = pd.crosstab(index = consumer_list, columns = "count") # Calculate the cell count
			N0[my_tab.index.values,k] = np.ravel(my_tab.values / assumptions['scale']) # Scale to biomass

		# Make data.frame
		N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

	# Monoculture plate
	elif assumptions['monoculture'] == True:
		N0 = np.eye(plate_N.shape[0]) *assumptions['n_inoc']/assumptions['scale']
		N0 = pd.DataFrame(N0, index = consumer_index, columns = ["W" + str(i) for i in range(plate_N.shape[0])])
	
	return N0

def sample_from_pool2(plate_N, assumptions, synthetic_community_size = 2, n = None):
	"""
	Make synthetic communities with given initial richness
	"""
	S_tot = plate_N.shape[0] 
	N0 = np.zeros((plate_N.shape))
	consumer_index = plate_N.index
	well_names = plate_N.columns
	
	if n is None:
		n = assumptions['n_inoc']
		
	for k in range(plate_N.shape[1]):
		consumer_list = np.random.choice(S_tot, size = synthetic_community_size, replace = False) 
		
		for v in range(synthetic_community_size):
				N0[consumer_list[v], k] = n / synthetic_community_size / assumptions["scale"]

	N0 = pd.DataFrame(N0, index = consumer_index, columns = well_names)

	return N0

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

def passage_monoculture(plate_mono, f, scale = None, refresh_resource = True):
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

def make_medium(plate_R, assumptions):
	"""
	Design medium for the plate
	if assumptions['rich_medium'] == True, make rich medium
	"""
	if assumptions['rich_medium'] == True:
		np.random.seed(1)
	
		# Total number of resource in this universe
		R_tot = plate_R.shape[0] 
	
		# Make empty plate
		R0 = np.zeros((plate_R.shape)) # Make empty plate
	
		# Resource index
		resource_index = plate_R.index 
	
		# Well index
		well_names = plate_R.columns
	
		resource_pool = np.random.uniform(0, 1, size = R_tot) # Uniform distribution
		resource_pool = resource_pool/np.sum(resource_pool)
		resource_list = np.random.choice(R_tot, size = assumptions["R0_food"], replace = True, p = resource_pool) # Draw from the pool
		my_tab = pd.crosstab(index = resource_list, columns = "count")
		food_compostion = np.ravel(my_tab.values)
		for i in range(plate_R.shape[1]):
			R0[my_tab.index.values,i] = food_compostion
		R0 = pd.DataFrame(R0, index = resource_index, columns = well_names)
	else:
		R0 = plate_R
	return R0

def make_plate(assumptions, params):
	"""
	prepares the plate
	"""
	
	# Make dynamical equations
	def dNdt(N,R,params):
		return MakeConsumerDynamics(assumptions)(N,R,params)
	def dRdt(N,R,params):
		return MakeResourceDynamics(assumptions)(N,R,params)
	dynamics = [dNdt,dRdt]
	
	# Make initial state
	init_state = MakeInitialState(assumptions)
	
	plate = Metacommunity(init_state, dynamics, params, scale = assumptions["scale"], parallel = False) 
	
	# Add media to plate (overrides community simulator)
	plate.R = make_medium(plate.R, assumptions)
	plate.R0 = make_medium(plate.R0, assumptions)  
	  
	# If plate is to be replaced by overwritting plate, skip the sampling
	if pd.isnull(assumptions["overwrite_plate"]):
		plate.N = sample_from_pool(plate.N, assumptions)
	
	return plate

# Data operation

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

def overwrite_plate(plate, assumptions):
	""" 
	Overwrite the plate N, R, and R0 dataframe by the input composition file
	"""
	import os
	assert(os.path.isfile(assumptions['overwrite_plate'])), "The overwrite_plate does not exist"
	# Read the input data file
	df = pd.read_csv(assumptions["overwrite_plate"])
	
	# By default, use the latest transfer to avoid well name conflict
	df = df[df.Transfer == np.max(df.Transfer)]
	
	# If only one community, repeat filling this community into n_wells wells
	if len(df["Well"].unique()) == 1:
		temp_df = df.copy()
		for i in range(assumptions["n_wells"]):
			temp_df["Well"] = "W" + str(i)
			temp_df.assign(Well = "W" + str(i))
			df = pd.concat([df, temp_df])
	# If the input overwrite file has multiple communities, check if it has the same number as n_wells
	assert len(df["Well"].unique()) == assumptions["n_wells"], "overwrite_plate does not have the same number of wells as n_wells"
	# Check if the input file type has consumer, resurce and R0
	assert all(pd.Series(df["Type"].unique()).isin(["consumer", "resource", "R0"])), "overwrite_plate must have three types of rows: consumer, resource, R0"
	# Make empty dataframes
	N = plate.N.copy()
	R = plate.R.copy()
	R0 = plate.R.copy()
	# N0
	for w in range(assumptions["n_wells"]):
		temp_comm = df[(df["Well"] == ("W" + str(w))) & (df["Type"] == "consumer")][["ID", "Abundance"]]
		temp = np.zeros(N.shape[0])
		for i in range(temp_comm.shape[0]):
			temp[int(temp_comm.iloc[i]["ID"])] = temp_comm.iloc[i]["Abundance"]
			N["W" + str(w)] = temp

	# R
	for w in range(assumptions["n_wells"]):
		temp_res = df[(df["Well"] == ("W" + str(w))) & (df["Type"] == "resource")][["ID", "Abundance"]]
		temp = np.zeros(R.shape[0])
		for i in range(temp_res.shape[0]):
			temp[int(temp_res.iloc[i]["ID"])] = temp_res.iloc[i]["Abundance"]
			R["W" + str(w)] = temp
	# R0
	for w in range(assumptions["n_wells"]):
		temp_R0 = df[(df["Well"] == ("W" + str(w))) & (df["Type"] == "R0")][["ID", "Abundance"]]
		temp = np.zeros(R0.shape[0])
		for i in range(temp_R0.shape[0]):
			temp[int(temp_R0.iloc[i]["ID"])] = temp_R0.iloc[i]["Abundance"]
			R0["W" + str(w)] = temp
	plate.N = N
	plate.N0 = N
	plate.R = R
	plate.R0 = R0
	
	# Passaage the overwrite plate
	if assumptions["passage_overwrite_plate"]:
		plate.Passage(np.eye(assumptions["n_wells"]) * assumptions["dilution"])
	
	return(plate)
