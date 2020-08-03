#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 26 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp
import pandas as pd


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
    algorithm_types = ["community_phenotypes", "selection_algorithms", "perturbation_algorithms"]
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
    
    
def make_protocol(params_simulation, protocol_name, selection_algorithm = None, repeated_selection = False):
    """
    Make protocol for one experimental protocol 
    """
    temp_df = pd.DataFrame({
        "algorithm_name": protocol_name,
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": "no_selection"
        })
    if protocol_name != "simple_screening":
        if repeated_selection: 
            temp_df["selection_algorithm"] = [selection_algorithm for i in range(params_simulation["n_transfer_selection"])] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])]
        elif repeated_selection == False:
            temp_df["selection_algorithm"] = ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + [selection_algorithm] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])]
    
    return temp_df
    

def make_algorithms(params_simulation):
    """
    Make a comprehensive dataframe of all protocols 
    """
    
    
    # Control
    simple_screening = make_protocol(params_simulation, "simple_screening")
    select_top25 = make_protocol(params_simulation, protocol_name = "select_top25", selection_algorithm = "select_top25percent", repeated_selection = False)
    select_top10 = make_protocol(params_simulation, protocol_name = "select_top10", selection_algorithm = "select_top10percent", repeated_selection = False)
    pool_top25 = make_protocol(params_simulation, protocol_name = "pool_top25", selection_algorithm = "pool_top25percent", repeated_selection = False)
    pool_top10 = make_protocol(params_simulation, protocol_name = "pool_top10", selection_algorithm = "pool_top10percent", repeated_selection = False)
    
    # Experimental protocols
    Blouin2015 = make_protocol(params_simulation, protocol_name = "Blouin2015", selection_algorithm = "pool_top10percent", repeated_selection = True)
    Blouin2015_control = make_protocol(params_simulation, protocol_name = "Blouin2015_control", selection_algorithm = "pool_top10percent_control", repeated_selection = True)
    Chang2020a = make_protocol(params_simulation, protocol_name = "Chang2020a", selection_algorithm = "select_top16percent", repeated_selection = True)
    Chang2020a_control = make_protocol(params_simulation, protocol_name = "Chang2020a_control", selection_algorithm = "select_top16percent_control", repeated_selection = True)
    Chang2020b = make_protocol(params_simulation, protocol_name = "Chang2020b", selection_algorithm = "select_top25percent", repeated_selection = True)
    Chang2020b_control = make_protocol(params_simulation, protocol_name = "Chang2020b_control", selection_algorithm = "select_top25percent_control", repeated_selection = True)
    Jochum2019 = make_protocol(params_simulation, protocol_name = "Jochum2019", selection_algorithm = "pool_top10percent", repeated_selection = True)
    Mueller2019 = make_protocol(params_simulation, protocol_name = "Mueller2019", selection_algorithm = "pool_top25percent", repeated_selection = True)
    Panke_Buisse2015 = make_protocol(params_simulation, protocol_name = "Panke_Buisse2015", selection_algorithm = "pool_top28percent", repeated_selection = True)
    Swenson2000a = make_protocol(params_simulation, protocol_name = "Swenson2000a", selection_algorithm = "pool_top20percent", repeated_selection = True)
    Swenson2000a_control = make_protocol(params_simulation, protocol_name = "Swenson2000a_control", selection_algorithm = "pool_top20percent_control", repeated_selection = True)
    Swenson2000b = make_protocol(params_simulation, protocol_name = "Swenson2000b", selection_algorithm = "select_top25percent", repeated_selection = True)
    Swenson2000b_control = make_protocol(params_simulation, protocol_name = "Swenson2000b_control", selection_algorithm = "select_top25percent_control", repeated_selection = True)
    Swenson2000c = make_protocol(params_simulation, protocol_name = "Swenson2000c", selection_algorithm = "pool_top20percent", repeated_selection = True)
    Wright2019 = make_protocol(params_simulation, protocol_name = "Wright2019", selection_algorithm = "pool_top10percent", repeated_selection = True)
    Wright2019_control = make_protocol(params_simulation, protocol_name = "Wright2019_control", selection_algorithm = "pool_top10percent_control", repeated_selection = True)
    
    # Sub-lineage protocols
    Arora2019 = make_protocol(params_simulation, protocol_name = "Arora2019", selection_algorithm = "Arora2019", repeated_selection = True)
    Arora2019_control = make_protocol(params_simulation, protocol_name = "Arora2019_control", selection_algorithm = "Arora2019_control", repeated_selection = True)
    Raynaud2019a = make_protocol(params_simulation, protocol_name = "Raynaud2019a", selection_algorithm = "Raynaud2019a", repeated_selection = True)
    Raynaud2019a_control = make_protocol(params_simulation, protocol_name = "Raynaud2019a_control", selection_algorithm = "Raynaud2019a_control", repeated_selection = True)
    Raynaud2019b = make_protocol(params_simulation, protocol_name = "Raynaud2019b", selection_algorithm = "Raynaud2019b", repeated_selection = True)
    Raynaud2019b_control = make_protocol(params_simulation, protocol_name = "Raynaud2019b_control", selection_algorithm = "Raynaud2019b_control", repeated_selection = True)
    
    # Theory
    Penn2004 = make_protocol(params_simulation, protocol_name = "Penn2004", selection_algorithm = "Williams2007a", repeated_selection = True)
    Williams2007a = make_protocol(params_simulation, protocol_name = "Williams2007a", selection_algorithm = "Williams2007a", repeated_selection = True)
    Williams2007b = make_protocol(params_simulation, protocol_name = "Williams2007b", selection_algorithm = "Williams2007b", repeated_selection = True)
    Xie2019a = make_protocol(params_simulation, protocol_name = "Xie2019a", selection_algorithm = "select_top_dog", repeated_selection = True)
    Xie2019b = make_protocol(params_simulation, protocol_name = "Xie2019b", selection_algorithm = "select_top10percent", repeated_selection = True)
    
    
    #directed_selection
    directed_selection = pd.DataFrame({
        "algorithm_name": "directed_selection",
        "transfer": range(1, params_simulation["n_transfer"] + 1),
        "community_phenotype": params_simulation["selected_function"],
        "selection_algorithm": ["no_selection" for i in range(params_simulation["n_transfer_selection"]-1)] + ["select_top"] + ["no_selection" for i in range(params_simulation["n_transfer"] - params_simulation["n_transfer_selection"])]
    })
    
    algorithms = pd.concat([
        # Control
        simple_screening, select_top25, select_top10, pool_top25, pool_top10,
        # Experimental protocols
        Blouin2015, Blouin2015_control, Chang2020a, Chang2020a_control, Chang2020b, Chang2020b_control, 
        Jochum2019, Mueller2019, Panke_Buisse2015, 
        Swenson2000a, Swenson2000a_control, Swenson2000b, Swenson2000b_control, Swenson2000c,
        Wright2019, Wright2019_control,
        # Sub-lineage protocols
        Arora2019, Arora2019_control, Raynaud2019a, Raynaud2019a_control, Raynaud2019b, Raynaud2019b_control, 
        # Theory
        Penn2004, Williams2007a, Williams2007b, Xie2019a, Xie2019b,
        directed_selection
        ])

    
    return algorithms
    
    
    
    
    
    
    

