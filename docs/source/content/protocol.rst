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







