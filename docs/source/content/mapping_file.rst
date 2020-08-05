Input Mapping File
==================

The input mapping ``.csv`` lists 69 essential parameters in columns and (indepdendent) selection experiments in rows.

An example is here [Attach a example csv file]. 


.. contents:: The mapping file has four main types of parameters

File operation
---------------

.. confval:: selected_function

    :type: string
    :default: ``f1_additive``

    Function under selection. Available options are ``f1_additive`` and ``f2_interaction``, ``f2a_interaction``, ``f3_additive_binary``, ``f4_interaction_binary``, ``f5_invader_growth``, and ``resource_distance_community_function``.


.. confval:: protocol

    :type: string
    :default: ``simple_screening``

    Protocol to implement. Only the protocols listed in ``E_protocols.py`` can be used.
    

.. confval:: seed

    :type: integer
    :default: ``1``

    Random seed to initiate pseudorandom number generator.


.. confval:: exp_id

    :type: string
    :default: ``f1_additive-simple_screening-1``

    Experiment-specific ID, which will also determine the naming convention of output files. For example, the community function is saved in ``f1_additive-simple_screening-1_function.txt`` if ``save_function=True``, whereas community compostition is saved in ``f1_additive-simple_screening-1_compostition.txt`` if ``save_composition=True``.


.. confval:: overwrite_plate

    :type: string
    :default: ``NA``

    To replace the initial plate composition with an arbitrary plate, specify a text file of the community composition that containes four columns: Type, ID, Well, and Abundance. If an output text file (e.g., ``f1_additive-simple_screening-1_compostition.txt``) is specified and it contains composition for more than two transfers, by default only the metacommunity compostition of the latter tranfer is read.


.. confval:: passage_overwrite_plate

    :type: boolean
    :default: ``False``

    If overwrite_plate != NA, set TRUE if the community from overwrite_plate is at equilibrium and need an addititonal transfer.


.. confval:: output_dir

    :type: string
    :default: ``data/``
    
    Directory where the output files will be stored. 


.. confval:: save_function

    :type: boolean
    :default: ``True``
    
    Set True to save function data. 


.. confval:: save_composition

    :type: boolean
    :default: ``True``
    
    Set ``True`` to save composition data.


.. confval:: save_plate

    :type: boolean
    :default: ``False``
    
    Set ``True`` to save initial Metacommunity in a ``pickle`` file.


.. confval:: function_lograte

    :type: integer
    :default: ``1``
    
    How often you save the function in transfers. Default is saving functional data from every transfer.

.. confval:: composition_lograte

    :type: integer
    :default: ``20``
    
    How often do you save the composition in transfers. 
    
| 

Protocol-specific parameters
----------------------------

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
    
    Set ``True`` to generate a rich medium sampled from an uniform distribution. Set ``False`` to generate a minimal medium with only the first resource is supplied. 


.. confval:: monoculture

    :type: boolean
    :default: ``False``
    
    Set ``True`` to run simple screening with all monocultures from the regional species pool. The number of wells is equal to the number of species in the regional pool.


.. confval:: dilution

    :type: float
    :default: ``0.001``
    
    Dilution factor in the batch culture.


.. confval:: n_wells

    :type: integer
    :default: ``96``
    
    Number of wells (communities) in a plate (metacommunity).


.. confval:: n_propagation

    :type: float
    :default: ``1``
    
    Incubation time of a transfer. 
    

.. confval:: n_transfer

    :type: integer
    :default: ``40``
    
    Number of total transfers (generations) to be run in the protocol.


.. confval:: n_transfer_selection

    :type: interger
    :default: ``20``
    
    Number of transfers (generations) that consecutively executes selection matrices from the start of an experiment. The number of stabilizaiton transfer equals to the difference between ``n_transfer_total`` and ``n_transfer_selection``.

|

Species contribution to function    
--------------------------------

.. confval:: sigma_func

    :type: float
    :default: ``1``
    
    Standard deviation for drawing speices-specific per-capita contribution to community function.


.. confval:: alpha_func

    :type: float
    :default: ``1``
    
    Contribution of species interaction to community function relative to the additive case.


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
    
    Standard deviation of fraction of cost feeded into a gamma distribution. ``cost_sd = 0`` if ``cost_mean = 0``, ``cost_sd = 0.01`` if ``cost_mean > 0``.


|

Directed evolution
------------------

.. confval:: directed_selection

    :type: boolean
    :default: ``False``
    
    Set ``True`` to run directed selection, one of flags below in directed evolution has to be also set ``True``.


.. confval:: knock_out

    :type: boolean
    :default: ``False``
    
    Set ``True`` to perform knock out pertubation.


.. confval:: knock_in

    :type: boolean
    :default: ``F``
    
    Set ``True`` performs knock in pertubation. 


.. confval:: knock_in_threshold

    :type: float 
    :default: ``0.95``
    
    If ``knock_in = True``, use the default ``knock_in_threshold=0.95``, which means that top 5% species in the pool is prepared to be knocked in a community, whereas the rest 95% of are not used.


.. confval:: bottleneck

    :type: boolean
    :default: ``False``
    
    Set ``True`` to perform bottleneck pertubations.


.. confval:: bottleneck_size

    :type: float
    :default: ``0.00001``
    
    If ``bottleneck=T``, perform an bottleneck shock to the specified communities by a dilution factor default to ``bottleneck_size=0.00001``. This bottleneck dilutoon is in addition to the regular dilution factor in the batch culture ``dilution=0.001``.


.. confval:: migration

    :type: boolean
    :default: ``False``
    
    Set ``True`` to perform migration pertubations.


.. confval:: n_migration

    :type: integer
    :default: ``1000000``
    
    Number of cells in the migrant community.


.. confval:: s_migration

    :type: integer
    :default: ``NA``
    
    Number of species in the migrant community. If ``NA`` (as default), the migrant community is sampled from a regional pool where the species abundance follows power-law distribution. If set into an integer, ``n_migration`` cells will be equally allocated to ``s_migrations`` species from the pool to build the migrant community.


.. confval:: coalescence

    :type: boolean
    :default: ``False``
    
    Set ``True`` to perform coalescence pertubation.


.. confval:: f_coalescence

    :type: float
    :default: ``0.5``
    
    Between 0 and 1. Fraction of migrant community during coalescence. The fraction of a perturbed community is ``1-f_coalescence``. 


.. confval:: resource_shift

    :type: boolean
    :default: ``False``
    
    Set ``True`` performs resource pertubations.


.. confval:: r_type

    :type: string
    :default: ``add``
    
    Type of resource pertubation. Available options are ``rescale_add``, ``rescale_remove``, ``add``, ``remove``, ``old``. A fraction ``r_percent`` of resource A is removed, and that amount of resource is added to another resource B.


.. confval:: r_percent

    :type: float
    :default: ``1``
    
    Fraction of specified resource that is removed. ``r_percent=1`` means all resource A is removed. 

|

Community-simulator parameters
-------------------------------

The parameters in this section are inherited and some with differnt values from community-simulator.

.. confval:: sampling

    :type: string
    :default: ``Binary_Gamma``
    
    Specify choice of sampling algorithm to generate the consumer uptake rate vector. Options are ``Gaussian``,``Binary``,``Gamma``, ``Binary_Gamma``.


.. confval:: sn

    :type: integer
    :default: ``2100``
    
    Number of microbial species in the global pool.


.. confval:: sf

    :type: integer
    :default: ``1``
    
    Number of specialist family.


.. confval:: s_gen

    :type: integer
    :default: ``0``
    
    Number/Richness of generalist taxa.


.. confval:: rn

    :type: integer 
    :default: ``90``
    
    Number of resource types. 


.. confval:: rf     

    :type: integer
    :default: ``1``
    
    Number of resource classes.


.. confval:: R0_food

    :type: float
    :default: ``1000``
    
    Total resource abundance.
    
    
.. confval:: food

    :type: float
    :default: ``1000``
    
    Index of food source being supplied in the minimal medium. Only works when ``rich_medium=False``.


.. confval:: supply

    :type: string
    :default: ``off``
    
    Choice of intrinsic resoruce dynamics. Set ``off`` for batch culture where resource is not renewing within a transfer. 
    

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


.. confval:: sparsity
    
    :type: float
    :default: ``0.2``
    
    Sparsity of metabolic matrix.


.. confval:: fs
    
    :type: float
    :default: ``0.45``
    
    Fraction of secreted byproducts allocated to the same resource class.


.. confval:: fw

    :type: float
    :default: ``0.45``
    
    Fraction of secreted byproducts allocated to waste resource class.


.. confval:: g
    
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


.. confval:: m

    :type: float
    :default: ``0``
    
    Minimal energy uptake for maintenance of species i (energy/time). Mortality.


.. confval:: n

    :type: integer 
    :default: ``2``
    
    Hill coefficient for functional response (unitless).


.. confval:: response
    
    :type: string
    :default: ``type III``
    
    Functional response of uptaking rates.


.. confval:: sigma_max
    
    :type: float 
    :default: ``1``
    
    Maximum input flux (mass/time) for type III functional response.


.. confval:: regulation
    
    :type: string
    :default: ``independent``
    
    Metabolic regulation.


.. confval:: nreg
    
    :type: integer
    :default: ``10``
    
    Hill coefficient that tunes steepness of metabolic regulation.


.. confval:: tau
    
    :type: float
    :default: ``1``
    
    External resource supply rate when ``supply="external"`` for chemostat setting.


.. confval:: r
    
    :type: string
    :default: ``independent``
    
    Renewal rate for self renewing resources when ``supply="self-renewing"`` for chemostat setting.


    
    
    
    
