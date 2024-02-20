""" Analysis Parameters / Config File
Definition of default parameters and dataset specific parameters in a python dictionary that is saved as yaml file. 
This is a non redundant representation of all analysis paths. 

Please stay with the standard folder structure. If you do this, you don't need to adjust the parameters:  
- proj_id
- version

"""

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
            'analysis_path': path.abspath('./../../../../'), 
            # Cases:
            # Your data is already placed correctly inside the folder structure and this file is placed in the same folder structure
            # -> same as analysis path: path.abspath('./../../../../') + '/<default>'
            # Your data is placed in the standard folder structure but on a mounted volume, while this file is placed locally 
            # -> 'pathTo/FolderAbove/ProjectFolder' + '/<default>'
            # Your data is not placed in the standard folder structure (be aware that fastq files will not be copied into the folder structure in contrast to other input data)
            # -> 'path/to/data'
            'data_root_path': '/Volumes/sdxxxxxx/guest/' + '/<default>', #path.abspath('./../../../../') + '/<default>',  # for example path to SDS mounted location
            # name of the metadata file
            'metadata': 'metadata.tsv',
            ###################
            ## Nf-Core Paths ##
            ###################
            # Path to where the reference genome is placed. If you choose different references from the provided ones, make sure to overwrite the nf-core params accordingly.
            'references_path': '/mnt/sds-hd/sdxxxxxx/projects/references',
            # Path to the project at the timepoint of nfcore pipeline execution. (When you work with the project locally but want to run the nf-core pipeline on the cluster, provide here the project path on the cluster)
            # The nf-core sample sheet will be based on this path so that the run script can find the fastq files.
            'exec_env_path': '/mnt/sds-hd/sdxxxxxx/guest/', 
            'nfcore':{
                # Before you can use the generated run script for the nf-core pipeline, you have to install nextflow. Provide here the path to the executable.  
                # This path is used by the run script
                'nextflow_executable': '/home/user_id/nextflow', 
        #'samplesheet_name': '', # name without filetype, it is saved as .csv
            }
        },
        'use_pickle_data': True,    # if True h5ad files are read and then saved as pickle, the pickle files are used from there on
        'cluster':{
            # Your email address for getting emails from the jobs
            'email': ''
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
            },
            'gene_sets':{
                'MSigDB':{ # list here the files that you downloaded, they should be placed in the prior knowledge folder
                    'files': ['m2.cp.reactome.v2023.2.Mm.symbols.gmt', 'm5.go.v2023.2.Mm.symbols.gmt', 'mh.all.v2023.2.Mm.symbols.gmt'],
                    'substr':['immune', 'inflammatory'] # only datasets that include one of these substrings are used
                }
            },
            'ligand_receptor': {
                'liana':{
                    'resources': [('mouseconsensus',)]
                }
            }
        },
        ###############
        ## Decoupler ##
        ###############
        'decoupler':{
            'methods': [('ulm', ),],
            'min_n': [5],
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
        'preprocessing': {
            'reference_genome': 'GRCm39', # not implemented
            'reference_genome_release': 109, # only implemented for getting gene names, not used for nfcore
            'gene_id': 'gene_id',
            'qc_cols_obs': ['n_genes_by_counts', 'total_counts'],
            'qc_cols_var': ['n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'],
            'basicFilt': {
                #'group': None, 
                # filter used for X attribute (filtered)
                # shortcut in foldername is 'prev' for 'prevalence' as it shows in how many samples a gene must be expressed
                # preprocessing is executed per value, the last value defines the results that stay in the adata object
                'large_n': [5, 10, 20], # used by dc filter expr
                # plot_filter_expr params and dc filter expr params
                'min_count': 10, # filter for raw attribute (counts)
                'min_total_count': 15, 
                'min_sample_size': 5, # not used by dc filter expr
                'min_cells': 4, # sc.pp.filter_genes -> a gene must be expressed in at least min_cells
                'min_genes': 1000, # sc.pp.filter_cells -> a cell must express at least min_genes
                'downsample_to': '', # for example 10000 reads per cell in sc data; if left empty, the size of the smallest sample/cell is used
                'pca_dot_size': 1500 # good for bulk/rather big
            }
        },
        #############################
        ## Differential Expression ##
        #############################
        # In this example config file we use the following diffExpr block for the bulk dataset. We don't use this block in any other dataset so we could move this block directly to the dataset specific params part of bulk dataset 01. 
        # The single cell dataset and the pseudobulk dataset overwrite the diffExpr block defined here.
        'diffExpr':{
            'conditions': ["sex_isCond_sample", "sex_subj", "isCond_sample", "sampleID"],
            'deseq2': {
                'design_factors': [
                    ['sex_isCond_sample', 'f_0'],   # leads to m1_f0, m0_f0, f1_f0
                        ['sex_isCond_sample', 'f_1'], 
                        ['sex_isCond_sample', 'm_0'],    # leads to m1_f1 as the other levels are sorted out because the zero is not at the end
                        ['sex_subj', 'f'],
                        ['isCond_sample', '0']
                ]
            }
        }
    },
    #####################
    ## Dataset Params ##
    #####################
    'dataset_params': {
        'mouse': { # organism (human, mouse)
            'scRNA': {
                '01':{
                    'rawType': 'log',
                    'xType': 'filtered',
                    'subs': {
                        'sampleID': ['sample', ('C1',)]
                    },
                    'nfcore':{'take_only': []},
                    'paths':{
                        'celltype_markers': 'markers.yaml',
                        'rawpath': 'raw_preprocessed/'
                    },
                    'preprocessing': {
                        'reference_genome': 'GRCm39', # not implemented
                        'reference_genome_release': 109, # only implemented for getting gene names, not used for nfcore
                        #'gene_id': 'gene_id',
                        'qc_cols_obs': ['n_genes_by_counts', 'total_counts'],
                        'qc_cols_var': ['n_cells_by_counts', 'mean_counts', 'pct_counts_mt', 'pct_dropout_by_counts', 'total_counts'],
                        'basicFilt': {
                            #'group': None, 
                            # filter used for X attribute (filtered)
                            # shortcut in foldername is 'prev' for 'prevalence' as it shows in how many samples a gene must be expressed
                            # preprocessing is executed per value, the last value defines the results that stay in the adata object
                            'large_n': [2], # used by dc filter expr
                            'min_sample_size': 2, # not used by dc filter expr
                            'min_cells': 3, # sc.pp.filter_genes -> a gene must be expressed in at least min_cells
                            'min_genes': 200, # sc.pp.filter_cells -> a cell must express at least min_genes
                            'max_genes': 6000,
                            'downsample_to': 75000, # for example 10000 reads per cell in sc data; if left empty, the size of the smallest sample/cell is used
                            'pca_dot_size': 500
                        }
                    },
                    'diffExpr':{
                        'conditions': ['isCond_sample', 'sampleID'],
                    }                   
                }
            },
            'bulkRNA':{
                '01_sc':{ # pseudobulk dataset
                    'rawType': 'counts',
                    'xType': 'filtered',
                    'paths':{
  
                    },
                    'priorKnowledge':{
                        'transcription_factors': { 
                            'collectri':{ # can be left out?

                            }
                        },
                        'pathways': { 
                            'progeny':
                                {
                                   
                                }
                        },
                        'gene_sets':{
                            'MSigDB':{

                            }
                         }
                    },
                    'preprocessing': {
                        'reference_genome': 'GRCm39', # not implemented
                        'reference_genome_release': 109, # only implemented for getting gene names, not used for nfcore
                        'gene_id': 'gene_id',
                        'qc_cols_obs': ['n_genes_by_counts', 'total_counts','psbulk_n_cells', 'psbulk_counts'], 
                        'qc_cols_var': ['n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'],
                        'basicFilt': {
                            #'group': None, 
                            # filter used for X attribute (filtered)
                            # shortcut in foldername is 'prev' for 'prevalence' as it shows in how many samples a gene must be expressed
                            # preprocessing is executed per value, the last value defines the results that stay in the adata object
                            'large_n': [5], # used by dc filter expr
                            # plot_filter_expr params and dc filter expr params
                            'min_count': 10, # filter for raw attribute (counts)
                            'min_total_count': 15, 
                            'min_sample_size': 5, # not used by dc filter expr
                            'min_cells': 4, # sc.pp.filter_genes -> a gene must be expressed in at least min_cells
                            'min_genes': 1000, # sc.pp.filter_cells -> a cell must express at least min_genes
                            'max_genes': '',
                            'downsample_to': '', # for example 10000 reads per cell in sc data; if left empty, the size of the smallest sample/cell is used
                            'pca_dot_size': 800 # good for bulk/rather big
                        }
                    },
                    'diffExpr':{
                        'conditions': ["cluster_isCond_sample"],
                        'deseq2': {
                            'design_factors': [
                                    ['cluster_isCond_sample', '0_0'],
                                    ['cluster_isCond_sample', '1_0'],
                                    ['cluster_isCond_sample', '2_0'],
                                    ['cluster_isCond_sample', '3_0']
                            ]
                        }
                    }
                },
                '01': {
                    ####################
                    ## Nf-Core Params ##
                    ####################
                    'hasUmi': False,
                    'rawType': 'counts',
                    'xType': 'filtered',
                    'nfcore':{
                        'pipeline':{
                            'take_only': ['rnaseq'],   # names of the nf-core pipelines
                            'rnaseq':{
                                'params':{
                                    'skip_trimming': True,
                                    'skip_deseq2_qc': True # because data is already cleaned by soapnuke
                                }
                            }
                        }
                    },
                    'paths':{
                        # Your raw data is expected to be in the data folder. Provide the path to the folder that contains the samples starting from within the data foolder. (path to folders named by sample and containing paired end fastq files)
                        'rawpath': 'raw/name_of_seq_run/',
                        # When you have different h5ad files for the same dataset, you can choose here which to take
                        #'data_root_filename': '01_opt0_merged_gene_counts.h5ad',
                        'nfcore':{
                            'samplesheet_name': 'samples',
                        }
                    },
                    'priorKnowledge':{
                        'transcription_factors': { 
                            'collectri':{ # can be left out?

                            }
                        },
                        'pathways': { 
                            'progeny':
                                {
                                    'top':['300', 'ADD']
                                }
                        },
                        'gene_sets':{
                            'MSigDB':{

                            }
                         }
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

