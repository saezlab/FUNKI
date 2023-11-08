# FUNKI <img style="float: right; margin-left: 3em" width="300" alt="funki logo" src="./images/funki_humanAndMice_whitebackground.jpg">

FUNKI is a web application that connects prior knowledge with your datasets. It is written in python and uses the package *streamlit* to build the graphical user interface.  
To provide feedback or ask for help, you can open an *issue* in this repository or write to *hanna.schumacher@uni-heidelberg.de* . 
  
</br>
</br>
</br>
</br>
<p align="center">
  <span>The following diagram shows the options that hide behind the dropdown menus in the user interface. Dotted lines show options that are visible but not implemented, yet.</span></br>   
  <span><i>Update: The phospho option is implemented. The single cell option is implemented in the development version only. The website will probably crash with single cell data. If this is the case, please install FUNKI locally.</i></span></br></br>
  <img width="680"src="./images/funki_ui.jpg">
</p>

</br>
  
<p align="center">
  <span> The following diagram explains how the input from the user interface is processed by FUNKI.   </span> </br></br>
  <img width="700"src="./images/funki_processing.jpg">
</p>

## Local Installation
If you want to install FUNKI locally, for example to use it without internet connection, you can do this in four simple steps: 
1. Install poetry for managing virtual python environments
2. Clone the repository 
3. Initiate the virtual environment by using poetry
4. Run FUNKI

The following code should do this. If you encounter any problems with the poetry installation, you can look through the [official instructions](https://python-poetry.org/docs/#installation). 
```bash
# install poetry 
curl -sSL https://install.python-poetry.org | python3 -
# if the poetry command is not found on your Mac 
# (for other systems, see the offical instructions)
echo 'export PATH="$HOME/.local/bin/:$PATH"' >> ~/.zshrc  
source ~/.zshrc   
# clone repo and enter it
git clone https://github.com/saezlab/FUNKI.git
cd FUNKI
# initiate virtual environment
poetry update
# run FUNKI
poetry run streamlit run Analysis.py
```  
  
## Usage of the *standard-workflows* Python package

While the Graphical User Interface of FUNKI is rather limited in its functionality so far, the code base behind FUNKI is much more developed. Inlcuded in the FUNKI repository is the *standard-workflows* Python package. This offers the following functionality for bulk-RNA and single cell RNA data.   

- Analysis project management
  - Config file (= analysis_params.py)
  - Automatic path definitions for input and output data based on the config file
  - An *analysis object* that holds *dataset objects*. Those hold the data, params and paths based on the config file.
  - A loop decorator and implemented loops for the utility functions (see below) for easy processing of multiple datasets at once
- Utility functions for running nf-core pipelines
  - Implemented for the rnaseq pipeline, other pipelines can be added in the same way
  - Pipeline configuration can be done with the config file (see above)
  - Needed files for the nf-core pipeline are created: sample sheet, config file, params file and execution script
  - Files are saved in the standardised folder structure. Different parameter combinations are saved separately.
  - Paths are provided (for example to the resulting counts file)
- Utility functions for running Decoupler on single cell RNA data
  - progeny/collectri/dorothea
  - mean activity calculation
  - output: csv files and plots as pdf (heatmaps, umaps)
- Utility functions for running DeSeq2 on bulk-RNA data (in progress)
  
### Example
*Caution: The following explanations assume some familiarity with programming. If you're not familiar with some of the used technologies (github, bash, poetry, python) and/or you'd prefer some assistance in using the functionalities of the *standard-workflows* package for your project, feel free to write to hanna.schumacher@uni-heidelberg.de. Working through this together via emails or video chat is no problem. On the contrary, it helps to improve the package and the instructions.*  
  
#### Environment setup
1. Clone FUNKI    
2. Create a new project folder next to your other projects (in this example we call it 'PROJ') and set up your poetry environment:   
  ```bash
  cd PROJ
  poetry init
  ```
3. Add the standard-workflows package to your pyproject.toml file by using the relative or full path to where it is placed:  
  ```python
  standard-workflows = {path = "/Users/name/Documents/projects/FUNKI/Workflows/standard-workflows", develop = true}
  ```

#### Main Script
Create a new jupyter notebook under *PROJ/v00/scripts/Python/main.ipynb*.   
*v00* stands for *version 0*. This allows us to be very explorative in the first analysis round and then move on to *v01* for the second round which could be the final version of the code and results for a paper. The input data (.h5ad) doesn't need to be copied from the first to the second version. The initialisation of the *analysis object* automatically searches for the data in the previous version when there is no data in the current one.    

Let's start with importing all functionality of the *standard_workflows* package:  
```python
# import
from standard_workflows import sc_analysis_baseclass as sc_classes
from standard_workflows import sc_analysis_loops as scl
from standard_workflows import sc_decoupler_utility as dcu
from standard_workflows import nfcore_utility as nfu
from standard_workflows import diffexpr_utility as diffexpr
from standard_workflows import sc_funcs
```
To start with a new project, the first step is to create a *dataset* that is used as template for your actual datasets. The *Analysis* class has a function that creates a *dataset object*. This object gets its functionality from the Analysis classes that it is inhering from. Every *dataset* must inherit from the class *Baseanalysis*. The other classes are optional based on the functionality that you want to use. Each additional class adds new functions, parameters and paths to your *datasets*.  

```python
dataset_template = sc_classes.Analysis.new_dataset(sc_classes.Baseanalysis, dcu.Decoupler, nfu.NfCore, diffexpr.DiffExpr) 
```
The second step is to initialise the *analysis object* by providing some details about your datasets. Each dataset gets described by a tuple in a list. Each tuple contains the name of the dataset, the sequencing type, the organism and the dataset_template that you want to use with this dataset. Additionally, you have to provide the path to your *analysis_params.py* file. This file is usually placed in *PROJ/v00/analysis/*. We'll come to this later.

```python
analysis = sc_classes.Analysis(datasets=[
            ('01', 'bulkRNA', 'mouse', dataset_template),
            ('02', 'bulkRNA', 'mouse', dataset_template)
            ], params_path = path.abspath("./../../analysis/"))
```  
The sc_analysis_loops module provides a couple of functions that wrap analysis functions into loops. Using these functions is very convenient as they automatically run over all *datasets* that the *analysis object* holds. You can enable this functionality by providing the sc_analysis_loops (scl) module with a reference to the *analysis object*.  

```python
scl.analysis = analysis
```
#### Config File (analysis_params.py)

Create a new python script here: *PROJ/v00/analysis/analysis_params.py*. While the name of the main script doesn't matter, the  name of this file does. You can use the following content for the file or just copy the analysis_params file from within the FUNKI repository. 

```python

from os import path
import yaml
analysis_params = {
#####################
## Project Params ##
#####################
    'proj_params': { 
        # name of the project
        'proj_id': path.basename(path.normpath(path.abspath("./../../../"))),    # always one folder above, projID should be 2 to 8 capital letters or underscores
        # Don't change the version. A new version must have a new config file. 
        'version': path.basename(path.normpath(path.abspath("./../../"))), 
        'paths': {
            # path to 'projects' folder or the folder where the proj results shall be saved
            'analysis_path': path.abspath('./../../../../'), #"/Volumes/sd22b002/guest/
            # Your data is already placed correctly inside the folder structure and this file is placed in the same folder structure
            # -> same as analysis path: path.abspath('./../../../../') + '/<default>'
            # Your data is placed in the standard folder structure but on a mounted volume, while this file is placed locally 
            # -> 'pathTo/FolderAbove/ProjectFolder' + '/<default>'
            # Your data is not placed in the standard folder structure (be aware that fastq files will not be copied into the folder structure in contrast to other input data)
            # -> 'path/to/data'
            'data_root_path': '/Volumes/sd22b002/guest/' + '/<default>', #path.abspath('./../../../../') + '/<default>',  # for example path to SDS mounted location: .../mounted/projects/
            # name of the metadata file
            'metadata': 'metadata.tsv',
            ###################
            ## Nf-Core Paths ##
            ###################
            # Path to where the reference genome is placed. If you choose different references from the provided ones, make sure to overwrite the nf-core params accordingly.
            'references_path': '/mnt/sds-hd/sd22b002/projects/references',
            # Path to the project at the timepoint of nfcore pipeline execution. (When you work with the project locally but want to run the nf-core pipeline on the cluster, provide here the project path on the cluster)
            # The nf-core sample sheet will be based on this path so that the run script can find the fastq files.
            'exec_env_path': '/mnt/sds-hd/sd22b002/guest/', 
            # Your raw data is expected to be in the data folder. Provide the path to the folder that contains the samples starting from within the data foolder. (path to folders named by sample and containing paired end fastq files)
            'rawpath': 'raw/F22FTSEUBEB0037_MUSrkqhR/soapnuke/clean/', 
            'nfcore':{
                # Before you can use the generated run script for the nf-core pipeline, you have to install nextflow. Provide here the path to the executable.  
                # This path is used by the run script
                'nextflow_executable': '/home/hd/hd_hd/hd_ac294/nextflow', # binac: /beegfs/work/hd_ac294/nextflow 
        #'samplesheet_name': '', # without filetype! It is saved as .csv
            }
        },
        'use_pickle_data': False,    # if True h5ad files are read and then saved as pickle, the pickle files are used from there on
        'cluster':{
            # Your email address for getting emails from the jobs
            'email': 'hanna.schumacher@uni-heidelberg.de'
        },
        ####################
        ## Nf-Core Params ##
        ####################
        'nfcore':{
            'pipeline':{                       
                # name of the nf-core pipeline
                'rnaseq':{
                    # used in the nf-core config file and must match the version given by nf-core
                    'version': '3.12.0',      
                    # used as foldername 
                    'version_name': 'v031200',  
                    'executor': 'slurm',
                    # check the modules of your cluster and write here the name of the java module
                    'module_java': 'devel/java_jdk/1.18', # helix; binac would be devel/java_jdk/11.0.4
                    # if more modules are needed, write them here. For example singularity.  
                    'load_modules': "module load system/singularity/3.11.3", 
                    # nf-core profile parameter
                    'profile': 'singularity', # BinAC cluster would have: binac,singularity as they have their own profile for nf-core jobs.
                    'process_config': "queue = 'single'", # queue = 'single' is needed for slurm on helix (otherwise you get the error that the time limit is not set correctly)
                    'max_cpus': 25,
                    'max_memory': '210.GB',
                    'max_time': '6.h'
                },
                'atacseq':{

                },
                'hic':{

                }
            }
        },
        #####################
        ## Prior Knowledge ##
        #####################
        'priorKnowledge':{
            'transcription_factors':{
                'collectri':{
                    'split_complexes': ['False']
                },
                #'dorothea':{   # levels
                #    'levels': [('A', 'B', 'C')]
                #}, 
            },
            'pathways': { 
                'progeny':{
                    'top': [500]
                }
            }
        },
        ###############
        ## Decoupler ##
        ###############
        'decoupler':{
            'methods': [('ulm', ),],
            #minsize
            #numberOfPermutations
            'meanacts': {
                'groupby': ['vars'],
                'minstd': [0.0]
            }
        },
        ###########
        ## Liana ##
        ###########
        'liana': {
            "methods": [["natmi", "connectome", "logfc", "sca", "cellphonedb"]],
            "base": "exp(1)",
            "lig_rec": [["all"]]
        }, 
        ###################
        ## Preprocessing ##
        ###################
        'preprocess': {
            'basicFilt': {
                'group': None, 
                'min_count': 10, 
                'min_total_count': 15, 
                'large_n': 1, 
                'min_prop': 1
            }
        },
        #############################
        ## Differential Expression ##
        #############################
        'diffExpr':{
            'deseq2': {

            }
        }
    },
#####################
## Dataset Params ##
#####################
    'dataset_params': {
        'mouse': { # organism (human, mouse)
            'scRNA': {
                'priorKnowledge':{
                    'collectri':{ # can be left out?

                    },
                    'progeny':
                        {
                            'top':['ADD']
                        }
                }
            },
            'bulkRNA':{
                '02': {
                    ####################
                    ## Nf-Core Params ##
                    ####################
                    'hasUmi': False,
                    'nfcore':{
                        'pipeline':{   # names of the nf-core pipelines
                            'rnaseq':{
                                'params':{
                                    'skip_trimming': True,
                                    'skip_deseq2_qc': True # because data is already cleaned by soapnuke
                                }
                            }
                        }
                    },
                    'paths':{
                        'rawpath': 'raw/F22FTSEUBEB0037_MUSrkqhR/soapnuke/clean/',
                        'nfcore':{
                            'samplesheet_name': 'samples_soapnuke_clean',
                        }
                    }
                },
                '01': {
                    'hasUmi': False,
                    'paths':{
                        'rawpath': 'raw/F22FTSEUBEB0037_MUSrkqhR/soapnuke/raw/',
                        'nfcore':{
                            'samplesheet_name': 'samples_soapnuke_raw', 
                        }
                    }
                },
            }
        }
    }  
}

```

#### Nf-Core Utilities Module
If you have filled in the needed parameters in the config file, you can use the following code to prepare the nf-core run(s):  
```python
  scl.prepare_nfcore('rnaseq')

  # Choose the dataset that you want to do the preparations for
  ds = analysis.datasets[0]

  # We call our basic parameter combination option 0 and run the init function
  ds.init_nfcore_run('opt0')
  # The default parameters are adjusted by the custom params from the config file.
  # The resulting parameter combination is saved in the nfc_custom_params attribute under 'opt0'.
  # Therefore, we can not just manipulate the parameters via the config file but also via code as shown below. 
  # In this case the parameters from opt0 are taken, but the aligner is changed to star_salmon.
  # We call this parameter combination opt1_salmon.
  ds.nfc_custom_params['opt1_salmon']= ds.nfc_custom_params['opt0'].copy()
  ds.nfc_custom_params['opt1_salmon'].update({"aligner": "star_salmon"})
  # Now we run the init function with the new parameter combination
  ds.init_nfcore_run('opt1_salmon')
  # The whole process generates new paths that we can save as follows:
  analysis.save_paths()
  # Double check the parameter combinations:
  ds.nfc_custom_params
  del ds
```

