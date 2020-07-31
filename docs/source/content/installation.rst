Installation
============


System requirement
------------------

* Python 3.7.3
* `community-simulator <https://github.com/Emergent-Behaviors-in-Biology/community-simulator>`_
* Scipy, Numpy, Pandas, Matplotlib, functools, itertools, random. CVXPY is not required for our simulations as batch-culture simulations does not use the steadystate method in community-simulator.

Ecoprospector package depends on [community-simulator](https://github.com/Emergent-Behaviors-in-Biology/community-simulator) (developed by the Mehta group and described in their [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230430)), which depends on Numpy, Pandas, Matplotlib, SciPy that are all included in Anaconda distribution. 


| 

Install the development version
-------------------------------

Clone the github repository  to a local directory ::

    $ git clone https://github.com/Chang-Yu-Chang/ecoprospector .

Then browse to the Ecoprospector directory and install package ::

    $ pip install -e .
