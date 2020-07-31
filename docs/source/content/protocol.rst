Selection Protocol
==================

.. contents:: 


Protocol overview
------------------

Similar to conducting a selection experiment in batch culture, *selection protocol* specifies how the experiment should be followed in each generation (transfer). Selection protocol is a DataFrame, where each row is a generation where communities are ranked by their *Selected function* and are passed to the next generation according to *Selection matrix*.

For example, a Swenson2000a protocol that repeatedly conducted for 20 will be like. 

[Placehold for the img]

A selection protocol is auto-generated from the imported parameters from mapping file ``.csv``. Currently Ecoprospector does not allow imported protocol. However, an advanced user should be able to make own protocol by adding a user-defined DataFrame in the source code.


Generation/Transfer
---------------------

A series of positive integers from 1 to the number of total generations, which is specified as ``n_transfer`` in the mapping :code:`.csv`. 

|

Selected function
------------------

The community funtion that is under selection at each generation. Inputs in this column should map to python functions in B_community_phenotype.py. 

[show list of options]

f1_additive
f2_interaction
distance from an arbitrary resource target (similar to Williams et al 2007)
resistance to a fast growing invader (similar to Wilson 1992)

An advanced user should be able to add functions that depend on any model and code it into the package as a function in B_community_phenotype.py

|

Selection matrix
----------------

Each cell has a string name that refers to a Python function that computes the selection matrix. A library of possible selection matrix options is saved in C_selection_matrix.py. The selection matrix represents what selection regime (i.e. which and what fraction of communities are selected, how much of a parent community is transferred into an offspring community) is used at the end of each generation. See the selection matrix section for how to make the selection matrix.




