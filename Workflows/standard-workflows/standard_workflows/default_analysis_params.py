""" Analysis Parameters 
Definition of default parameters and dataset specific parameters in a python dictionary that is saved as yaml file. 
This is a non redundant representation of all analysis paths. 

Please stay with the standard folder structure. If you do this, you don't need to adjust the parameters:  
- proj_id
- version

"""
import yaml
from os import path
# Caution: Do not split paths and analysis_params on the level of Analysis Obj. 
# Splitting leads to problems with merging and replacing. 
# 1. Execute additions and replacements
# 2. merge_dicts: proj_params are merged with the params of a specific dataset (dataset_params/organism/seq_type/datasetname)
# 3. Extend dataset specific paths

# ToDo: Implement use_pickle_data
analysis_params = {
    'proj_params': { 
        'proj_id': path.basename(path.normpath(path.abspath("./../../../"))),    # always one folder above, projID should be 2 to 8 capital letters or underscores
        'version': path.basename(path.normpath(path.abspath("./../../"))), 
        'paths': {
            'analysis_path': path.abspath('./../../../../'),   # path to 'projects' folder or the folder where the proj results shall be saved
            'data_root_path': path.abspath('./../../../../') + '/<default>'  # for example path to SDS mounted location: .../mounted/projects/
        },
        'use_pickle_data': True,    # h5ad files are read and then saved as pickle, the pickle files are used from there on
        'priorKnowledge':{
            'dorothea':{   # levels
                'levels': [('A', 'B', 'C')]
            },  
            'progeny':{
                'top': [300]
            }
        },
        'decoupler':{
            'methods': [('mlm', ),],
            #minsize
            #numberOfPermutations
            'meanacts': {
                'groupby': ['vars'],
                'minstd': [0.0]
            }
        },
        'liana': {
            "methods": [["natmi", "connectome", "logfc", "sca", "cellphonedb"]],
            "base": "exp(1)",
            "lig_rec": [["all"]]
        }
    },
    'dataset_params': {
        'mouse': { # organism (human, mouse)
            'scRNA': {
                '01':{},
                'priorKnowledge':{
                    'dorothea': 
                        {   
                            'levels': [('D', 'ADD')]
                        },  
                    'progeny':
                        {
                            'top':['ADD']
                        }
                }
            }
        }
    }  
}

def init():
    with open('./../../analysis/analysis_params.yaml', 'w+') as file: 
        yaml.dump(analysis_params, file, sort_keys=False)

if __name__ == '__main__':
    init()

