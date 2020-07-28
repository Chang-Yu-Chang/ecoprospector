# Ecoprospector

## Introduction

This package is designed for executing protocols on artificially selecting microbial communities of desired selected functions. The simulation is based on batch culture in which the microbial consumer-resource dynamics undergoes.

See our [preprint](https://www.biorxiv.org/content/10.1101/2020.07.24.214775v2) that uses Ecoprospector to simulate directed evolution of microbial communities.

## Dependencies

Ecoprospector package depends on [community-simulator](https://github.com/Emergent-Behaviors-in-Biology/community-simulator) (developed by the Mehta group and described in their [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230430)), which depends on Numpy, Pandas, Matplotlib, SciPy that are all included in Anaconda distribution. 

A python development setup by [Anaconda](https://docs.anaconda.com/anaconda/install/) will be sufficient to implement Ecoprospector.

## Installation

### Mac and Linux

Download the code or clone this github repository to a local directory 

`$ git clone https://github.com/Chang-Yu-Chang/Ecoprospector .`

Browse to the Ecoprospector directory and install package 

`$ pip install -e .`


## Package structure

![Ecoprospector](outline.png)

Ecoprospector depends on the input experimental protocols to artificially select a set of microbial communties (metacommunity). Each part of the package can be analogous to protocol steps in a bench experiment. 


- `input.csv` is a mapping file that contains all specified parameters used by Ecoproepector.

- Experimental setup are deployed before the simulation.
    
    - Selection protocol
    
    - Species features
    
    - Initial plate (metacommunity) composition

- Simulation deployes generation-wise protocol. A community generation includes three cyclic steps:

    - Batch culture: consumer-resource dynamics
    
    - Community function and rank: communities are ranked according to their traits
    
    - Selection matrix: selection regimes to seed the offspring communities are standardized by selection matricex

- Outputs: simple text files with specified naming conventions

    - Function: data frame with columns exp_id, Well, Transfer, CommunityTrait, Richness, Biomass
    
    - Composition: data frame with columns exp_id, Well, Transfer, Type, Abundance










