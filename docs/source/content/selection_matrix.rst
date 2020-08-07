Selection Matrix
===================

What is a selection matrix?
------------------------------------
A selection matrix is a map specifying how the parental communities are selected according to their function ranks, and how the selected communities are pooled or distributed to seed the offspring communities. It is a square matrix of size ``n_wells``. The columns are ranked parental communities and the rows are offspring communities. Each element in the selection matrix specifies the dilution factor used for the batch culture. 


[insert the image of an example of a selection matrix]

The selection matrices allow us to standardize most strategies of artificial community selection, for example, propagule and migrant pool approaches, into a regular form.

[Example of selection matrices of propagule and migrant pool strategies]


How does selection matrix work in ecoprospector?
----------------------------------------------------------------------

A selection matrix is written in the form of a Python function to accommodate a varied number of communities in different independent experiments. These functions take a vector of values (the default output of :ref:`Community Function` functions) as input. The selection matrix function will read the length of the input vector, and construct a selection matrix of that length. The selection matrix is then used to guide the passaging of the metacommunity.

A selection matrix must be defined during the simulation setup, i.e. stored in the ``C_selection_matrices.py``. During simulation, any particular selection matrix will be called according to :ref:`Selection Protocol`.

A library of selection matrices
----------------------------------------------------------------------

We saved all the predefined selection matrices in ``C_selection_matrices.py``. These selection matrices were adapted from the selection protocols in the prior empirical and theoretical studies. 


