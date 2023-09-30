"""Contains the classes for initializing an analysis. 
Creating a new *Analysis* leads to reading the parameters, initializing the *Datasets*. 
A dataset can inherit from various ...analysis classes but inherits always from the *Baseanalysis* class contained in this file."""
from os import path, makedirs
from .sc_analysis_loops import loop
from .scfunctions import replace_dictvalues, merge_dicts, add_nested_key
from . import memoize
from copy import deepcopy
import scanpy as sc
import dill, yaml

############################
#### Baseanalysis Class ####
############################

class Baseanalysis:
    """The class that every *Dataset* inherits from."""
    def __init__(self, name, seq_type, organism, analysis_params, paths):
        """ Adds basic analysis information to a dataset (name, ..., params, paths)

        Pseudocode:   
            1. Set dataset properties (organism, seq_type, name)
            2. Add dataset_params:organism:seq_tpye to ap
            3. If ds name %in% ap:dataset_params:organism:seq_type: 
                merge ap:proj_params into it
                merge ap:proj_params into dataset paths
            Else 
                take ap and paths from ap:proj_params
            4. Add dataset properties to ap
            5. Add paths
            6. Add datapath 
                If data_root_path ends with '<default>' replace end with default path
                Else take path as is
            7. If data_root_filename in ap, use it
            Else use default datasetname
            8. Add paths
            9. Call super()      
        """

        # Dataset properties
        self.organism = organism
        self.seq_type = seq_type
        self.name     = name

        add_nested_key(analysis_params, ['dataset_params', self.organism, self.seq_type]) # important when analysis_params doesn't contain the dataset, yet
        
        # Combine the "project_params" with the "dataset_params".
        self.analysis_params = deepcopy(analysis_params["proj_params"])
        ap_ds_org_seqtype = analysis_params["dataset_params"][self.organism][self.seq_type]
        if(self.name in ap_ds_org_seqtype):
            self.analysis_params = merge_dicts(self.analysis_params,           ap_ds_org_seqtype[self.name])
            self._paths          = merge_dicts(analysis_params["proj_params"], ap_ds_org_seqtype[self.name])["paths"]
        else:
            #self.analysis_params = merge_dicts(self.analysis_params, {})
            self._paths = deepcopy(paths)
        
        # Add dataset properties to analysis_params
        self.analysis_params["organism"] = self.organism
        self.analysis_params["seq_type"] = self.seq_type
        self.analysis_params["name"]     = self.name
        
        # middlepath is the path between the path to the proj location and the file. 
        middlepath         = path.join(self.analysis_params["proj_id"], self.analysis_params["version"], "analysis", self.organism, self.seq_type) # MBEN_T/v00/analysis/human/sn
        middlepath_dataset = path.join(middlepath, self.name)             # MBEN_T/v00/analysis/human/sn/all

        datafoldername = 'data'

        # Add more paths
        for dictentryname, foldername in zip(["datapath_tmp", "figpath", "resultpath", "loggingpath", "subsetspath"],  [datafoldername, "figures", "results", "logging", "subsets"]):
            self._paths.update({dictentryname: path.join(self._paths["analysis_path"], middlepath_dataset, foldername)})
        
        # Set datapath (depends on data_root_path)
        if not path.basename(path.normpath(self._paths["data_root_path"])) == "<default>": #if path not ends with <default>
            self._paths.update({"datapath": self._paths["data_root_path"]}) # -> datapath is same as data_root_path, no default folder structure
        else: 
            self._paths["data_root_path"] = path.dirname(path.normpath(self._paths["data_root_path"])) # remove "<default>" 
            self._paths.update({ 
                "datapath": path.join(self._paths["data_root_path"], middlepath_dataset, datafoldername) # add MBEN_T/v00/analysis/human/sn/all/data
                })  
            self._paths["data_root_path"] = self._paths["datapath"]

        # Set data_root_filename (either given one or default)
        if "data_root_filename" in self._paths.keys(): 
            self._paths["data_root_filename"] = self._paths["data_root_filename"] 
        else: 
            self._paths["data_root_filename"] = self.name + ".h5ad"    

        # Set datafilepath_tmp based on use_pickle_data
        fileextension=""
        if bool(self.analysis_params['use_pickle_data']) == True:
            fileextension = ".pickle"
        else: 
            fileextension = ".h5ad"
        datafilepath_tmp = path.join(self._paths["datapath_tmp"], f"{self.name}{fileextension}")

        self._paths.update({
            "datasetpath":        middlepath_dataset, # path from projectname to datasetname: MBEN_T/v00/analysis/human/sn/all
            "datafilepath":       path.join(self._paths["datapath"], self._paths["data_root_filename"]),
            "datafilepath_tmp":   datafilepath_tmp,
            "priorknowledge":     path.join(self._paths["analysis_path"], middlepath, "priorKnowledge"), # TODO: change this to storage path 
            "priorknowledge_tmp": path.join(self._paths["analysis_path"], middlepath, "priorKnowledge") 
        })
        
        self.data = "" # is set in __init__ of analysis obj
        super().__init__()
    
    def get_paths(self) -> dict:
        return self._paths

########################
#### Analysis Class ####
########################  
class Analysis: 
    """The first step in an analysis is to create an *Analysis* object with this class."""
    @memoize.Memoize
    def new_dataset(*bases):
        """Initialises new dataset by inheriting from the given classes."""
        class Dataset(*bases):
            def __init__(self, name, seq_type, organism, analysis_params, paths):
                """Set up of new dataset. 

                Args:
                    name (str): name/id of the dataset, like "all", "tumor", "pat1To4"
                    seq_type (str): seq_type of sequencing data, sn, sc, atac, 16s etc.
                    organism (str): "human" or "mouse"
                    analysis_params (dict): from "Analysis" obj
                    paths (dict): from "Analysis" obj
                """                
                super().__init__(name, seq_type, organism, analysis_params, paths)
        return Dataset

    @staticmethod
    def read_analysis(proj, proj_path, analysisversion):
        dill.load(path.join(proj_path, "projects", proj, "results", "analysis_{analysisversion}.pickle"))


    ### Init functions ###

    def __init__(self, datasets, params_path):       
        """Initializes analysis.
        Paths are generated and dataset objects are created. 

        Args:
            datasets (tuple): It contains datasetname (str), sequencingseq_type (str), organism (str), datasetClass (class)
            params_path (str): path to analysis params
        """
        # Read analysis_params.yaml
        with open(path.join(params_path, "analysis_params.yaml")) as stream:
            self.analysis_params = yaml.load(stream, Loader=yaml.BaseLoader)

        self.analysis_params = replace_dictvalues(self.analysis_params)

        # Init paths from the information in analysis_params
        self.__paths = deepcopy(self.analysis_params["proj_params"]["paths"])
        # Create datasets
        self.datasets = [constructor(name, seq_type, organism, self.analysis_params, deepcopy(self.__paths)) for name, seq_type, organism, constructor in datasets]        
        self.init_datasets()

    def get_paths(self) -> dict:
        return self.__paths

    def save_paths(self):
        """ Saves analysis params and paths from all datasets to yaml file. """
        self.analysis_params["default"]["paths"] = self.get_paths()
        for data in self.datasets:
            data.analysis_params["dataset_params"][data.organism][data.seq_type][data.name]["paths"] <- data.get_paths()
            with open(path.join(data.get_paths()["datasetpath"], "analysis_params_extended.yaml"), "w+") as file:
                yaml.dump(data.analysis_params, file)
    
    def init_datasets (self) -> None :
        """ Reads and cleans datasets """
        @loop(self.datasets, True)
        def init(self:Analysis) -> None:
            """ Read hf5ad data if no pickle exists. Reads into 'data' property of dataset. 

            ### Pseudocode
            ----------
            read datafilepath_tmp (either pickle or h5ad file)
                else read datafilepath
                    else assume new version copy from mounted path of old version to new version, then read
                save as pickle file
            """
            # Create tmp variables
            datapath = self._paths["datapath"] 
            datapath_tmp = self._paths["datapath_tmp"]
            datafilepath_tmp = self._paths["datafilepath_tmp"]
            datafilepath = self._paths["datafilepath"]
            print(self._paths)
            # Read data
            fileextension = path.splitext(datafilepath_tmp)[1]
            if(path.exists(datafilepath_tmp)):
                if  fileextension == '.pickle':
                    with open(datafilepath_tmp, "rb") as f:
                        self.data = dill.load(f)
                    print("Data was read in from pickle file.")
                elif fileextension == '.h5ad':
                    self.data = sc.read(datafilepath_tmp, cache = True)
                else: 
                    print(f"Please make sure that datafilepath_tmp ({datafilepath_tmp}) either ends with '.pickle' or with '.h5ad'.")
            else:
                if not path.exists(datapath_tmp):
                    makedirs(datapath_tmp)
                if datapath != '':
                    try: 
                        self.data = sc.read(datafilepath, cache = True)
                    except OSError as e:
                        # probably a new version number
                        if not path.exists(self._paths["datapath"]):
                            # get previous version number
                            import re
                            datapath = str(self._paths["datapath"])
                            r = re.compile(".*v(..).*")
                            numb = int((r.match(datapath)).group(1)) - 1
                            numb = "v" + str(numb).zfill(2) 
                            datapath_old = re.sub(r"v..", numb, datapath)                        
                            # take data from old version
                            import distutils.dir_util
                            distutils.dir_util.copy_tree(datapath_old, datapath)
                            # retry
                            self.data = sc.read(datafilepath, cache = True)
                        else:
                            print(f"Datafilepath_tmp ({datafilepath_tmp}) does not exist. Datapath ({datapath}) exists but datafilepath ({datafilepath}) does not. If you wanted to create a new version, delete datapath.")
                    if  fileextension == '.pickle':
                        with open(datafilepath_tmp, "wb") as dill_file:
                            dill.dump(self.data, dill_file)
                        print("Data was saved as pickle file.")
                    elif fileextension == '.h5ad':
                        self.data = sc.write(datafilepath_tmp)
                        print("Data was saved as h5ad file.")
                    else: 
                        print(f"Please make sure that datafilepath_tmp ({datafilepath_tmp}) either ends with '.pickle' or with '.h5ad'.")
         
                    print("Data was read in from datapath and is now saved in datapath_tmp.")
                else: 
                    "Please be aware that no data was read in as no data_root_path was provided."
        init()
        self.clean_datasets()

    ### Processing functions ###       
    def clean_datasets(self) :
        """ correct datatypes """
        @loop(self.datasets, True)
        def clean (dataset) : 
            if "seurat_clusters" in dataset.data.obs.columns : dataset.data.obs["seurat_clusters"] = dataset.data.obs.seurat_clusters.astype("category")
                        # to be sure that the sources and targets from the model match with the names used in the data, convert to lowercase
            dataset.data.var_names = [x.lower() for x in list(dataset.data.var_names)]
            # dataset.data.raw.var.index = dataset.data.raw.var._index
            dataset.data.raw.var.index = [x.lower() for x in list(dataset.data.raw.var.index)]
        clean()


    ### Save & Load ###
    # These functions are for saving an intermediate status as pickle file. 
    def save_datasets(self):
        with open(path.join(self.paths["datapath"], "datasets.pickle"), "wb") as dill_file:
            dill.dump(self.datasets, dill_file)

    def load_datasets(self):
        return dill.load(path.join(self.paths["datapath"], "datasets.pickle"))


                     



