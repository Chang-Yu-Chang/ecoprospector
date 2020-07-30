Selection Protocol
==================

Protocols overview
------------------

A selection protocol is basically a data frame with three columns (Generation, Selected trait, and Selection matrix) organized in a generation-wise manner.

A selection protocol is auto-generated from the imported parameters, which specifies genearal protocol number of generations. In the mapping :code:`.csv`, one needs to specify the name of the protocol to call it in Ecoprospector. For example, Swenson2000a is a protocol that repeats a "select-top-20percent-and-pool" selection regime for a couple of genearations. 

An advanced user can also make own protocol, which explicitly imported from the mapping :code:`.csv` but rather

To implement a selection protocol (hereafter protocol) on a set of microbial communities, an experimenter has to determine the community function under selection, which and what fraction of communities meet the criteria of selection, and whether and how to mix those selected communities to seed the next generation. These elements of experimental protocol are put together into a single DataFrame in Ecoprospector to capture the selection regimes in a generation-wise manner (FigSXX). The protocol DataFrame includes three columns: Generation, Selected trait, and Selection matrix.

All predefined protocols used in the main text are saved in the E_protocol.py. During the simulation of a protocol, Ecoprospector will calculate the community functions of target and return a vector of community function (length of metacommunity size). The resulting vector of community function is then fed into the Python function that specifies the selection matrix, which determines the transfer regime from parent metacommunity to the offspring metacommunity.

This protocol format allows different community traits to be selected for at different generations. For example, one may iteratively impose different selection targets to ask whether phenotypic variation maintained by balancing selective forces can potentially increase the effectiveness of a selection protocol. Likewise, the communities can be subject to different selection matrices over generations in order to explore an enormous space of possible selection regimes. The modular design of the protocol allows users to flexibly tune the protocol for their needs.


Generation/Transfer
---------------------

A series of positive integers from 1 to the number of total generations, which is specified as **n_transfer** in the mapping :code:`.csv`. 

|

Selected function
------------------

Each cell is a string name with references to a Python function that computes trait value from the community composition. A library of possible trait options is saved in B_selected_traits.py. At the end of each generation the scores of targeted functions under selection will be used to rank the communities. Community ranks are required as input of a selection matrix. See the :ref:`Community Function` section for details.

|

Selection matrix
----------------

Each cell has a string name that refers to a Python function that computes the selection matrix. A library of possible selection matrix options is saved in C_selection_matrix.py. The selection matrix represents what selection regime (i.e. which and what fraction of communities are selected, how much of a parent community is transferred into an offspring community) is used at the end of each generation. See the selection matrix section for how to make the selection matrix.






Community Function
--------------------------------

The Selected_function column in the csv file specifies the community trait that is under selection at each generation. Inputs in this column should map to the name of a python function in B_community_traits. The python functions in this script take as an input the plate object and parameters, and output a list of the function for each community in the metacommunity.  Any function of this form is compatible with the package. As well as the community function mentioned in the main text. Ecoprospector  currently includes the ability to simulate selection for biological inspired  functions such as i) distance from an arbitrary resource target (similar to Williams et al 2007) an	d ii) resistance to a fast growing invader (similar to Wilson 1992) . Additional functions that depend on any model parameter can also be coded  into the package as a function in B_community_phenotype.py
