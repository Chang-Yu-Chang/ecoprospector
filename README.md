# Ecoprospector

> Design protocol for artificially selecting microbial communities

This package is designed for executing community-level selection protocols to engineer microbial communities toward desired selected functions. The simulation is based on batch culture in which the microbial consumer-resource dynamics undergoes.

See our [preprint](https://www.biorxiv.org/content/10.1101/2020.07.24.214775v2) that uses Ecoprospector to study a range of selection strategies of directed evolution of microbial communities.

![](outline.png)

## Installation

### Mac and Linux

Download the code or clone this github repository to a local directory 

```sh
$ git clone https://github.com/Chang-Yu-Chang/Ecoprospector .
```

Browse to the Ecoprospector directory and install package 
```sh
$ pip install -e .
```

### Windows


## Usage example

With the mapping file (csv), executing one experiment is simple as 

```sh
$ Ecoprospector mapping_file.csv 0
```

For more examples and usage, please refer to the [documentation][https://ecoprospector.readthedocs.io/en/latest/].

## Development setup

Ecoprospector package depends on [community-simulator](https://github.com/Emergent-Behaviors-in-Biology/community-simulator) (developed by the Mehta group and described in their [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230430)), which depends on Numpy, Pandas, Matplotlib, SciPy that are all included in Anaconda distribution. 

A python development setup by [Anaconda](https://docs.anaconda.com/anaconda/install/) will be sufficient to implement Ecoprospector.


## Release History

* 0.0.1
    * Work in progress

## Documentation

Ecoprospector's documentation lives at [ecoprospector.readthedocs.io](https://ecoprospector.readthedocs.io/en/latest/)

## Meta

Chang-Yu Chang – [@changyu_chang](https://twitter.com/changyu_chang) – chang-yu.chang@yale.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/Chang-Yu-Chang/Ecoprospector](https://github.com/Chang-Yu-Chang/Ecoprospector)

## Contributing

1. Fork it (<https://github.com/Chang-Yu-Chang/Ecoprospector/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

