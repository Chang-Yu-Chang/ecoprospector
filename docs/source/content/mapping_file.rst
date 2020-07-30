Input Mapping File
==================

The input mapping file centralizes parameters in columns and lists (indepdendent) selection experiments in rows, allowing trackable changes among different experiments.

Each line of the csv file contains 66 entries that fully determine how the simulation is to be conducted. Despite the seemingly large number of input parameters, many of these do not need to be specified in a typical run (and can be simply set as NA). In addition many are boolean flags determining the type of simulation to be conducted. In Table S2-S5 we provided a full list of parameters as described in the main text. Table S6 provides a comprehensive list of all 66 inputs in the csv file which includes the string form of most parameters in Table S2-S5.

The mapping file has four main types of parameters:

* File operation
* Protocol-specific parameters
* Species contribution to function
* Directed selection
* Parameters inherited from community-simulator


File operation
---------------

* **selected_trait**: Function under selection
* **protocol**:  Protocol to implement
* **seed**: Integer. Random seed to initiate pseudorandom number generator
* **exp_id**: experiment-specific ID, which will also determine output filenames
* **overwrite_plate**: If not NA, there must be a text file. Overwrite the initial plate with this composition saved in this text file"
* **passage_overwrite_plate**: If overwrite_plate != NA, set TRUE if the overwrite_plate is at equilibrium and need an addititonal transfer"
* **output_dir**: Output directory. Default is :code:`data` subfolder
* **save_function**: Logical. Set True to save function data
* **save_composition**: Logical. Set True to save composition Data
* **save_plate**: Logical. Set True to save initial plate
* **function_lograte**: How often you save the function in transfers
* **composition_lograte**: How often do you save the composition in transfers
    
| 

Protocol-specific parameters
----------------------------

* **a**: Exponent parameter in power-law distribution that determines the species abundance in regional pool
* **scale**: Number of cells equivalent to N_i = 1
* **n_inoc**: Number of cells in the initial inoculum
* **rich_medium**: Logical. Whether to generate a rich medium sampled from a a random distribution or a minimal media with only a single resource
* **monoculture**: Logical. Whether to run simple screening with all monocultures from pool
* **d**: Dilution factor in the batch culture
* **t_incubation**: Incubation time
* **n_wells**: Number of wells; number of communities
* **n_transfer_total**: Number of total transfers (generations)
* **n_transfer_selection**: Number of selection transfers (generations)

|

Species contribution to function    
--------------------------------

* **sigma_func**: Standard deviation for drawing specifc speices/interaction function", 1,
* **alpha_func**: Relative functional contribution of species interaction to the additive case", 1,
* **binary_threshold**: Threshold for binary functions", 1,
* **g0**: The baseline conversion factor of biomass per energy", 1,
* **cost_mean**: Mean fraction of cost feeded into a gamma distribution. Suggested up to 0.05", 0,
* **cost_sd**: Sd of fraction of cost feeded into a gamma distribution. cost_sd = 0 if cost_mean = 0, cost_sd= 0.01 if cost_mean >0", 0,

|

Directed evolution
------------------

* **directed_selection**: If True, directed selection", F,
* **knock_out**: If True performs knock out pertubations", F,
* **knock_in**: If True performs knock in pertubation", F,
* **knock_in_threshold**: The percentile determining the high-performing species in the species pool used to knock in", 0.95,
* **bottleneck**: If True perform bottleneck pertubations", F, 
* **d_bottleneck**: Bottleneck size", 10^5, 
* **migration**: If true perform migration pertubations", F,
* **n_migration**: Number of cells in the migrant community", 10^6, 
* **coalescence**: If true perform coalescence pertubation", F,
* **f_coalescence**: Mixing ratio of coalescence; biomass of immigrant community relative to that of a perturbed community copy", 0.5,
* **resource_shift**: If true performs resource pertubations", F,
* **type_resource**: Type of resource pertubation. rescale_add, rescale_remove, add, remove, old. if NA defaults to resource swap", "add",
* **p_resource**: Tunes the magnitude of resource perturbation. The fraction from depleting a resource and move the same amount to another", 1,
    
Community-simulator parameters
-------------------------------

* **sampling**: {'Gaussian','Binary','Gamma', 'Binary_Gamma'} specifies choice of sampling algorithm", "Binary_Gamma",
* **sn**: Number of microbial species in global pool", 2100,
* **sf**: Number of specialist family
* **s_gen**: Number of generalists", 0,
* **rn**: Number of resources", 90,
* **rf**: Number of resource classes", 1,
* **R0_food**: Total resource abundance", 1000,
* **muc**: Mean sum over a row of the preference matrix ciα", 10,
* **sigc**: Standard deviation of sum over a row of the preference matrix ciα", 3,
* **c0**: Low consumption level for binary ciα", 0,
* **c1**: High consumption level for binary ciα", 1,
* **q**: Fraction of consumption capacity allocated to preferred resource class", 0,
* **s**: Sparsity of metabolic matrix", 0.2,
* **fw**: Fraction of secreted byproducts allocated to waste resource class",	0.45,
* **fs**: Fraction of secreted byproducts allocated to the same resource class", 0.45,
* **gi**: Conversion factor from energy uptake to growth rate (1/energy)", 1,
* **w**: Energy content of resource α (energy/mass)", 1,
* **l**: Leakage fraction", 0,
* **mi**: Minimal energy uptake for maintenance of species i (energy/time)", 0,
* **n**: Hill coefficient for functional response (unitless)", 2,
* **m**: Mortality", 0,
* **response**: Functional response", "type III",
* **sigma_max**: Maximum input flux (mass/time)", 1
