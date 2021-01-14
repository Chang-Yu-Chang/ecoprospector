# ecoprospector

> Simulating protocols for the artificial selection of microbial communities

This package is designed to simulate arbitrary community-level selection protocols on microbial metacommunitities in order to determine whether they can effectively engineer communitities with desired functions. The simulation are based on batch culture and the microbes in each community interact via consumer-resource dynamics.

See our [preprint](https://www.biorxiv.org/content/10.1101/2020.07.24.214775v2) that uses ecoprospector to study a range of selection strategies that direct the evolution of microbial communities.

![](outline.png)



# Installation

## Anaconda

A python development setup by [Anaconda](https://docs.anaconda.com/anaconda/install/) will be sufficient to implement community-simulator and ecoprospector.

### Mac and Linux

Download the code or clone the github repository of community simulator to a local directory and browse to the community-simulator directory and install the package

```sh
$ cd <your_local_directory>
$ git clone https://github.com/Emergent-Behaviors-in-Biology/community-simulator .
$ pip install -e .
```

Download the code or clone this github repository to a local directory and browse to the ecoprospector directory and install package 

```sh
$ cd <your_local_directory>
$ git clone https://github.com/Chang-Yu-Chang/ecoprospector .
$ pip install -e .
```

### Windows

The parallelization features in community-simulator are not currently supported on Windows and as such we cannot guarantee that the current version of ecoprospector will work in a windows environment. We would recommend using a linux emulator for windows such as Cygwin instead.

## pipenv

Alternative to Anaconda, one can choose to use `pipenv`.

```sh
$ pip install pipenv
$ brew install pipenv # Alternative for Mac users
```

Browse to the local directory and install dependency from Pipfile. Make sure that the Pipfile is in the current directory.

```sh
$ cd <your_local_directory>
$ pipenv install
```

Check dependency. ecoprospector should depend on community-selection and all other dependent packages.

```sh
$ pipenv graph
```

To activate the Pipenv shell:

```sh
$ pipenv shell
$ exit
```



## Usage example

With the mapping file (csv), executing one experiment is simple as 

```sh
$ ecoprospector mapping_file.csv 0
```

For more examples and usage, please refer to the [documentation](https://ecoprospector.readthedocs.io/en/latest/).


## Release History

* 0.0.2
    * Include other non-additive functions
* 0.0.1
    * Work in progress

## Documentation

Ecoprospector's documentation lives at [ecoprospector.readthedocs.io](https://ecoprospector.readthedocs.io/en/latest/)

## Meta

Chang-Yu Chang – [@changyu_chang](https://twitter.com/changyu_chang) – chang-yu.chang@yale.edu

Jean Vila – [@jccvila](https://twitter.com/jccvila) – jean.vila@yale.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/Chang-Yu-Chang/ecoprospector](https://github.com/Chang-Yu-Chang/ecoprospector)


