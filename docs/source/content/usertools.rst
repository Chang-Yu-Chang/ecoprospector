User Tools
==========

Main functions in ecoprospector



Make parameters
---------------

.. code-block:: python

    make_assumptions(input_csv, row_number)


.. confval:: input_csv

    :type: DataFrame
    :default: ``input_csv``

    mapping csv file

.. confval:: row_number

    :type: integer
    :default: ``0``

    The row number that specifies the experiment to run (0-indexed)

|

Prepare and set up expeirments
------------------------------

.. code-block:: python

    params, params_simulation , params_algorithm, plate = prepare_experiment(assumptions)


.. confval:: assumptions

    :type: List
    :default: ``assumptions``

    A comprehensive list read from the input csv file


|

Simulate the protocol
----------------------


.. code-block:: python

    simulate_community(params = params, params_simulation = params_simulation, params_algorithm = params_algorithm, plate = plate)


.. confval:: params

    :type: List
    :default: ``assumptions``

    A comprehensive list read from the input csv file


.. confval:: params_simulation

    :type: List
    :default: ``params_simulation``

    Parameters related to simulating batch culture 

.. confval:: params_algorithm

    :type: List
    :default: ``params_algorithm``

    Parameters related to protocol, community function,  selection matrices, and


.. confval:: plate

    :type: Metacommunity object
    :default: ``plate``

    Object defined in this project 
    


