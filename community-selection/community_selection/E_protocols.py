#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""

"""
This python script contains the different protocols (that combine selection and migration regimes
"""

import numpy as np
import scipy as sp
import pandas as pd

# Make library of algorithms
def make_algorithm_library():
    """
    Show the table of algorithms in this package
    """
    import re
    import pandas as pd
    
    # Find directory of community_selection modultes
    import community_selection
    module_dir = community_selection.__file__
    module_dir = re.sub("__init__.py", "", module_dir) 
    
    # 
    algorithm_types = ["community_phenotypes", "selection_algorithms", "migration_algorithms"]
    algorithms = list()
    
    for i in range(len(algorithm_types)):
    
        # Open files
        file_algorithm_phenotype = open(module_dir + ["B", "C", "D"][i] + "_" + algorithm_types[i] + ".py", "r")
        
        # Read lines
        line_list = list()
        line = file_algorithm_phenotype.readline()
        cnt = 1
        
        while line:
            line = file_algorithm_phenotype.readline()
            line_list.append(line.strip())
            cnt += 1
        
        # Regular expression
        algorithm_names = re.findall("def \w+", " ".join(line_list))
        list_algorithm = [re.sub("^def ", "", x) for x in algorithm_names]
        
        # Write the files
        algorithms.append(pd.DataFrame({"AlgorithmType": re.sub("s$", "", algorithm_types[i]), "AlgorithmName": list_algorithm}))
     
    return pd.concat(algorithms)
    
def make_algorithms(params_simulation):
    # Algorithms
    ## Simple screening
    simple_screening = pd.DataFrame({
        "algorithm_name": "simple_screening",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": "no_selection",
        "migration_algorithm": "no_migration"
    })

    ## Select top 25%
    select_top25 = pd.DataFrame({
        "algorithm_name": "select_top25",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top25percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  

    ## Select top 10%
    select_top10 = pd.DataFrame({
        "algorithm_name": "select_top10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top10percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    ## Pool top 25%
    pool_top25 = pd.DataFrame({
        "algorithm_name": "pool_top25",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top25percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  

    ## Pool top 10%
    pool_top10 = pd.DataFrame({
        "algorithm_name": "pool_top10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top10percent"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  

    # Select top 25% and bottleneck or pool top 25% and bottleneck
    select_top25_bottleneck_10 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_10"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    select_top25_bottleneck_100 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_100",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_100"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    select_top25_bottleneck_1000 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_1000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_1000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    select_top25_bottleneck_10000 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_10000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_10000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    select_top25_bottleneck_100000 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_100000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_100000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    select_top25_bottleneck_1000000 = pd.DataFrame({
        "algorithm_name": "select_top25_bottleneck_1000000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top_25_bottleneck_1000000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_10 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_10",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_10"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_100 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_100",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_100"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_1000 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_1000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_1000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_10000 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_10000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_10000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_100000 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_100000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_100000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    pool_top25_bottleneck_1000000 = pd.DataFrame({
        "algorithm_name": "pool_top25_bottleneck_1000000",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["pool_top_25_bottleneck_1000000"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })  
    
    # Arora2019
    Arora2019 = pd.DataFrame({
        "algorithm_name": "Arora2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top33percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Blouin2015
    Blouin2015 = pd.DataFrame({
        "algorithm_name": "Blouin2015",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Blouin2015 control
    Blouin2015_control = pd.DataFrame({
        "algorithm_name": "Blouin2015_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Chang2020a
    Chang2020a = pd.DataFrame({
        "algorithm_name": "Chang2020a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top16percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    # Chang2020a_control
    Chang2020a_control = pd.DataFrame({
        "algorithm_name": "Chang2020a_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top16percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })    
    
    # Chang2020b
    Chang2020b = pd.DataFrame({
        "algorithm_name": "Chang2020b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    # Chang2020b_control
    Chang2020b_control = pd.DataFrame({
        "algorithm_name": "Chang2020b_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Jochum2019
    Jochum2019 = pd.DataFrame({
        "algorithm_name": "Jochum2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Mueller2019
    Mueller2019 = pd.DataFrame({
        "algorithm_name": "Mueller2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Panke-Buisse2015
    Panke_Buisse2015 = pd.DataFrame({
        "algorithm_name": "Panke_Buisse2015",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top28percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Penn2004
    Penn2004 = pd.DataFrame({
        "algorithm_name": "Penn2004",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Raynaud2019a
    Raynaud2019a = pd.DataFrame({
        "algorithm_name": "Raynaud2019a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Raynaud2019b
    Raynaud2019b = pd.DataFrame({
        "algorithm_name": "Raynaud2019b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000a
    Swenson2000a = pd.DataFrame({
        "algorithm_name": "Swenson2000a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000a control
    Swenson2000a_control = pd.DataFrame({
        "algorithm_name": "Swenson2000a_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })

    # Swenson2000b
    Swenson2000b = pd.DataFrame({
        "algorithm_name": "Swenson2000b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000b_control
    Swenson2000b_control = pd.DataFrame({
        "algorithm_name": "Swenson2000b_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top25percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Swenson2000c
    Swenson2000c = pd.DataFrame({
        "algorithm_name": "Swenson2000c",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top20percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Williams2007a
    Williams2007a = pd.DataFrame({
        "algorithm_name": "Williams2007a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Williams2007b
    Williams2007b = pd.DataFrame({
        "algorithm_name": "Williams2007b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Williams2007b" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
        
    # Wright2019
    Wright2019 = pd.DataFrame({
        "algorithm_name": "Wright2019",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Wright2019_control
    Wright2019_control = pd.DataFrame({
        "algorithm_name": "Wright2019_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["pool_top10percent_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Xie2019a
    Xie2019a = pd.DataFrame({
        "algorithm_name": "Xie2019a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top_dog" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Xie2019b
    Xie2019b = pd.DataFrame({
        "algorithm_name": "Xie2019b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["select_top10percent" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    # Arora2019
    Arora2019 = pd.DataFrame({
        "algorithm_name": "Arora2019_V2",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Arora2019" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })    
    Arora2019_control = pd.DataFrame({
        "algorithm_name": "Arora2019_V2_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Arora2019_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })       
    # Raynaud2019a
    Raynaud2019a	= pd.DataFrame({
        "algorithm_name": "Raynaud2019a",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019a" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    Raynaud2019a_control	= pd.DataFrame({
        "algorithm_name": "Raynaud2019a_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019a_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    }) 
    # Raynaud2019b
    Raynaud2019b = pd.DataFrame({
        "algorithm_name": "Raynaud2019b",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019b" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    Raynaud2019b_control = pd.DataFrame({
        "algorithm_name": "Raynaud2019b_control",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["Raynaud2019b_control" for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    
    #directed_selection
    directed_selection = pd.DataFrame({
        "algorithm_name": "directed_selection",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])], 
        "migration_algorithm": "no_migration"
    })
    
    #long experiments
    directed_selection_long = pd.DataFrame({
        "algorithm_name": "directed_selection_long",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": np.tile(["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"],int(params_simulation["n_transfer"]/params_simulation["n_transfer_selection"])).tolist(), 
        "migration_algorithm": "no_migration"
    })  
    
    
    # Save the algorithms
    algorithms = pd.concat([
        # Controls
        simple_screening, 
        select_top25, select_top10, 
        pool_top25, pool_top10,
        select_top25_bottleneck_10, select_top25_bottleneck_100, select_top25_bottleneck_1000, 
        select_top25_bottleneck_10000, select_top25_bottleneck_100000, select_top25_bottleneck_1000000,
        pool_top25_bottleneck_10, pool_top25_bottleneck_100, pool_top25_bottleneck_1000, 
        pool_top25_bottleneck_10000, pool_top25_bottleneck_100000, pool_top25_bottleneck_1000000,
        # Literature
        Arora2019, Blouin2015, Blouin2015_control, Chang2020a, Chang2020a_control, Chang2020b, Chang2020b_control, 
        Jochum2019, Mueller2019, Panke_Buisse2015, Penn2004,
        Raynaud2019a, Raynaud2019b, Swenson2000a, Swenson2000a_control, Swenson2000b, Swenson2000b_control, Swenson2000c,
        Williams2007a, Williams2007b, Wright2019, Wright2019_control, Xie2019a, Xie2019b,
        Arora2019, Arora2019_control, Raynaud2019a, Raynaud2019a_control, Raynaud2019b, Raynaud2019b_control,
        # directed selection
        directed_selection,directed_selection_long])
    
    return algorithms
