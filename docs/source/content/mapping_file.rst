Input Mapping File
==================

The input mapping ``.csv`` lists 65 essential parameters in columns and (indepdendent) selection experiments in rows.

An example is here [Attach a example csv file]. 


.. contents:: The mapping file has four main types of parameters

File operation
---------------

.. confval:: selected_trait

    :type: string
    :default: ``f1_additive``

    Function under selection.


.. confval:: protocol

    :type: string
    :default: ``simple_screening``

    Protocol to implement.
    

.. confval:: seed

    :type: integer
    :default: ``1``

    Random seed to initiate pseudorandom number generator.


.. confval:: exp_id

    :type: string
    :default: ``f1_additive-simple_screening-1``

    Experiment-specific ID, which will also determine the output filenames


.. confval:: overwrite_plate

    :type: string
    :default: ``NA``

    To replace the initial plate composition with an arbitrary plate. It must be a text file that specifies the community composition for initial community, containin four columns: Type, ID, Well, and Abundance.


.. confval:: passage_overwrite_plate

    :type: boolean
    :default: ``False``

    If overwrite_plate != NA, set TRUE if the community from overwrite_plate is at equilibrium and need an addititonal transfer.


.. confval:: output_dir

    :type: string
    :default: ``data/``
    
    Output directory. 


.. confval:: save_function

    :type: boolean
    :default: ``True``
    
    Set True to save function data. 


.. confval:: save_composition

    :type: boolean
    :default: ``True``
    
    Set True to save composition data.


.. confval:: save_plate

    :type: boolean
    :default: ``False``
    
    Set True to save initial plate in a ``pickle`` file.


.. confval:: function_lograte

    :type: integer
    :default: ``1``
    
    How often you save the function in transfers.

.. confval:: composition_lograte

    :type: integer
    :default: ``20``
    
    How often do you save the composition in transfers.
    
| 

Protocol-specific parameters
----------------------------

.. confval:: a

    :type: float
    :default: ``0.01``
    
    Exponent parameter in power-law distribution that determines the species abundance in regional pool.


.. confval:: scale

    :type: integer
    :default: ``1000000``
    
    Number of cells equivalent to :math:`N_i = 1`.


.. confval:: n_inoc

    :type: integer
    :default: ``1000000``
    
    Number of cells in the initial inoculum.


.. confval:: rich_medium

    :type: boolean
    :default: ``True``
    
    Whether to generate a rich medium sampled from a a random distribution or a minimal media with only a single resource. Set True to allow rich medium.


.. confval:: monoculture

    :type: boolean
    :default: ``False``
    
    Whether to run simple screening with all monocultures from pool. Set True to run monocultures with the number of wells equal to the number of species in the regional pool.


.. confval:: d

    :type: float
    :default: ``0.001``
    
    Dilution factor in the batch culture.


.. confval:: t_incubation

    :type: float
    :default: ``1``
    
    Incubation time in a generation/transfer.


.. confval:: n_wells

    :type: integer
    :default: ``96``
    
    Number of wells (communities).
    

.. confval:: n_transfer_total

    :type: integer
    :default: ``40``
    
    Number of total transfers (generations).


.. confval:: n_transfer_selection

    :type: interger
    :default: ``20``
    
    Number of selection transfers (generations).

|

Species contribution to function    
--------------------------------

.. confval:: sigma_func

    :type: float
    :default: ``1``
    
    Standard deviation for drawing specifc speices/interaction function.


.. confval:: alpha_func

    :type: float
    :default: ``1``
    
    Relative functional contribution of species interaction to the additive case.


.. confval:: binary_threshold

    :type: float
    :default: ``1``
    
    Threshold for binary functions.


.. confval:: g0

    :type: float
    :default: ``1``
    
    The baseline conversion factor of biomass per energy.


.. confval:: cost_mean

    :type: float
    :default: ``0``
    
    Mean fraction of cost feeded into a gamma distribution. Suggested maximum to 0.05.


.. confval:: cost_sd

    :type: float
    :default: ``0``
    
    Standard deviation of fraction of cost feeded into a gamma distribution. cost_sd = 0 if cost_mean = 0, cost_sd= 0.01 if cost_mean > 0.


|

Directed evolution
------------------

.. confval:: directed_selection

    :type: boolean
    :default: ``False``
    
    If set True, run directed selection. One of the other flag in directed evolution selection has to be set True.


.. confval:: knock_out

    :type: boolean
    :default: ``False``
    
    If True performs knock out pertubation.


.. confval:: knock_in

    :type: boolean
    :default: ``F``
    
    If True performs knock in pertubation. 


.. confval:: knock_in_threshold

    :type: float 
    :default: ``0.95``
    
    The percentile determining the high-performing species in the species pool used to knock in. Default means top 5% species in the pool is prepared to be knocked in a community, whereas the rest 95% of are not used.


.. confval:: bottleneck

    :type: boolean
    :default: ``False``
    
    If True perform bottleneck pertubations.


.. confval:: d_bottleneck

    :type: float
    :default: ``0.00001``
    
    Bottleneck size.


.. confval:: migration

    :type: boolean
    :default: ``False``
    
    If True perform migration pertubations.


.. confval:: n_migration

    :type: integer
    :default: ``1000000``
    
    Number of cells in the migrant community.


.. confval:: s_migration

    :type: integer
    :default: ``NA``
    
    Number of species in the migrant community. If NA, the migrant community is sample from a regional pool where the species abundance follows power-law distribution. If set into an integer, ``n_migration`` cells will be equally allocated to ``s_migrations`` species.


.. confval:: coalescence

    :type: boolean
    :default: ``False``
    
    If True perform coalescence pertubation.


.. confval:: f_coalescence

    :type: float
    :default: ``0.5``
    
    Mixing ratio of coalescence; The fraction of immigrant community relative to that of a perturbed community. copy.


.. confval:: resource_shift

    :type: boolean
    :default: ``False``
    
    If True performs resource pertubations.


.. confval:: type_resource

    :type: string
    :default: ``add``
    
    Type of resource pertubation. rescale_add, rescale_remove, add, remove, old. If NA defaults to resource swap, 


.. confval:: p_resource

    :type: float
    :default: ``1``
    
    Tunes the magnitude of resource perturbation. The fraction from depleting a resource and move the same amount to another", 1,

|

Community-simulator parameters
-------------------------------

The parameters in this section are inherited and some with differnt values from community-simulator.

.. confval:: sampling

    :type: string
    :default: ``Binary_Gamma``
    
    Specify choice of sampling algorithm to generate the consumer uptake rate vector. Options are 'Gaussian','Binary','Gamma', 'Binary_Gamma'.


.. confval:: sn

    :type: integer
    :default: ``2100``
    
    Number of microbial species in global pool.


.. confval:: sf

    :type: integer
    :default: ``1``
    
    Number of specialist family.


.. confval:: s_gen

    :type: integer
    :default: ``0``
    
    Number of generalists.


.. confval:: rn

    :type: integer 
    :default: ``90``
    
    Number of resources.


.. confval:: rf     

    :type: integer
    :default: ``1``
    
    Number of resource classes.


.. confval:: R0_food

    :type: float
    :default: ``1000``
    
    Total resource abundance.
    

.. confval:: muc

    :type: float
    :default: ``10``
    
    Mean sum over a row of the preference matrix ciα.
    
    
.. confval:: sigc

    :type: float
    :default: ``3``
    
    Standard deviation of sum over a row of the preference matrix ciα.
    
    
.. confval:: c0

    :type: float
    :default: ``0``
    
    Low consumption level for binary ciα.
    
    
.. confval:: c1

    :type: integer
    :default: ``1``: 
    
    High consumption level for binary ciα.


.. confval:: q

    :type: float
    :default: ``0``
    
    Fraction of consumption capacity allocated to preferred resource class.


.. confval:: s
    
    :type: float
    :default: ``0.2``
    
    Sparsity of metabolic matrix.


.. confval:: fw

    :type: float
    :default: ``0.45``
    
    Fraction of secreted byproducts allocated to waste resource class.


.. confval:: fs
    
    :type: float
    :default: ``0.45``
    
    Fraction of secreted byproducts allocated to the same resource class.


.. confval:: gi
    
    :type: float
    :default: ``1``
    
    Conversion factor from energy uptake to growth rate (1/energy).


.. confval:: w
    
    :type: float 
    :default: ``1``
    
    Energy content of resource α (energy/mass).


.. confval:: l

    :type: float
    :default: ``0``
    
    Leakage fraction.


.. confval:: mi

    :type: float
    :default: ``0``
    
    Minimal energy uptake for maintenance of species i (energy/time).


.. confval:: n

    :type: integer 
    :default: ``2``
    
    Hill coefficient for functional response (unitless).


.. confval:: m
    
    :type: float
    :default: ``0``
    
    Mortality rate.


.. confval:: response
    
    :type: string
    :default: ``type III``
    
    Functional response.


.. confval:: sigma_max
    
    :type: float 
    :default: ``1``
    
    Maximum input flux (mass/time)
