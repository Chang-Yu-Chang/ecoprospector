# Ecoprospector

## Introduction

This package is designed for executing protocols on artificially selecting microbial communities of desired selected functions. The simulation is based on batch culture in which the microbial consumer-resource dynamics undergoes.

See the [preprint](https://www.biorxiv.org/content/10.1101/2020.07.24.214775v2) that uses `Ecoprospector` to simulate directed evolution of microbial communities.

## Dependencies

`Ecoprospector` package depends on [community-simulator](https://github.com/Emergent-Behaviors-in-Biology/community-simulator) (developed by the Mehta group and this package was described in their [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230430)), which depends on Numpy, Pandas, Matplotlib, SciPy that are all included in Anaconda distribution. 

A python development setup by [Anaconda](https://docs.anaconda.com/anaconda/install/) will be sufficient to implement `Ecoprospector`.

## Installation

### Mac and Linux

Clone this github repository to local using 
`$ git clone https://github.com/Chang-Yu-Chang/Ecoprospector`

Browse to the `Ecoprospector` directory in Terminal, and enter 

`$ pip install -e .`

