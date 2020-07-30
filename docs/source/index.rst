Ecoprospector's Tutorial
=========================================

This site covers Ecoprospector 

What is Ecoprospector?
======================

Ecoprospector is a Python package designed to simulate protocols of artificial selection on microbial communities, using a single ::code::`mapping_file.csv` as input (with specified parameters) and specified experiment (the first row in this case)::

    $ Ecoprospector mapping_file.csv 0

| 

Main features
=============

It aims to flexibly adapt major componenets of an experiemntal protocol, including features like batch culture, multi-well plate, selection strategy, and allows an user to design arbitrarily protocols with own needs.
 
* **Single input file**: Ecoprospector takes a ::code::`.csv` file where each row of the csv file specifies the parameters for one experiment.
* **Consumer-resource dynamics**: virtual microbial species with idiosyncratic metabolic properties interact with others in a community through secretion and uptakes. 
* **Batch-culture**: the community generation is divided by serial batch culture.
* **Community function**: various community functions can be applied, inclduing additive, intereative function.
* **Selection matrix**: the selection regimes (i.e., which parental communitues to select and how to seed the offspring communities) are standardized by selection matrix.
* **Modular protocol design**: the feature mentioned above can be assembled into a user-manual designed protocol. 

Our package is designed with three types of user in mind. 

* **Beginners** who have no python experience would be able to re-run all the simulations presented in this paper using the csv file alone and should be able to repeat pre-coded experimental protocols under varying parameter choices. 
* **Intermediate users** who have basic knowledge of python, should be able to code up their own protocols and may also be able to perform simple extensions to the package (such as introducing new types of community function, or selection matrices). 
* **Advanced users** who are familiar with python coding should be able to add additional functionality to the package, including carrying over several features of community-simulator that are not currently in use. This includes but is not limited to  a) introducing intrinsic resource dynamics (for chemostat simulations) b) alternative dynamical models such as Lotka-Volterra models.

|

Key contributors
================

Jean Vila and Chang-Yu Chang (both at Yale working with `Alvaro Sanchez <http://www.sanchezlaboratory.com/>`_) started to build Ecoprospector in collaboration with students from Physical Biology of Cells Course at Marine Biology Laboratoy in Woods Hole, who provided feedback on the modular design of Ecoprospector.

|

Documentation
=====================
.. toctree::
   :maxdepth: 2
   
   content/installation
   content/quickstart
   content/usertools
   content/mapping_file
   content/protocol

