Quick Start Guide 
=================

This page provides a quick start guide for using Ecoprospector with an single input :code:`.csv` file. If you do not know how to configurate :code:`.csv` file, please read :ref:`Input Mapping File` for more details on each parameter.

We created a commandline tool so you can run Ecoprospector , For example in Terminal on Mac, enter

.. code-block:: bash

    $ Ecoprospector mapping_file.csv 0

Where mapping_file.csv is the input :code:`csv` file and i is the row (0-indexed) specifying the experiment to be run. 

You can also run the above code in python. The line above is equivalent as:

.. code-block:: python

    from community_selection import *
    from community_selection.usertools import *
    assumptions = make_assumptions("mapping_file.csv", 0)
    params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)
    simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm, plate = plate)

The functons are described in :ref:`User Tools`