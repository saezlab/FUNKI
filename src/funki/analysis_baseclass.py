"""Contains the classes for initializing an analysis. 
Creating a new *Analysis* leads to reading the parameters, initializing the *Datasets*. 
A dataset can inherit from various analysis classes but inherits always from the *Baseanalysis* class contained in this file."""
from os import path, makedirs
from .analysis_loops import loop
from .utility_functions import replace_dictvalues, merge_dicts, add_nested_key, print_paths
from . import memoize
from copy import deepcopy
import scanpy as sc, pandas as pd
from abc import ABC
from IPython.display import display, Markdown 
import dill, yaml, json, yaml , pathlib
    
class AnalysisI(ABC):
        data:type
        analysis_params:dict
        paths:dict


############################
#### Baseanalysis Class ####
############################

class Baseanalysis(AnalysisI):
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
            self.paths          = merge_dicts(analysis_params["proj_params"], ap_ds_org_seqtype[self.name])["paths"]
        else:
            #self.analysis_params = merge_dicts(self.analysis_params, {})
            self.paths = deepcopy(paths)
        
        # Add dataset properties to analysis_params
        self.analysis_params["organism"] = self.organism
        self.analysis_params["seq_type"] = self.seq_type
        self.analysis_params["name"]     = self.name
        
        # middlepath is the path between the path to the proj location and the file. 
        middlepath         = path.join(self.analysis_params["proj_id"], self.analysis_params["version"], "analysis", self.organism, self.seq_type) # MBEN_T/v00/analysis/human/sn

        middlepath_dataset = path.join(middlepath, self.name)             # MBEN_T/v00/analysis/human/sn/all
        datafoldername = 'data'

        # Add more paths
        for dictentryname, foldername in zip(["datapath_tmp", "figpath", "resultpath", "loggingpath"],  [datafoldername, "figures", "results", "logging"]):
            self.paths.update({dictentryname: path.join(self.paths["analysis_path"], middlepath_dataset, foldername)})
        self.paths['subsetspath'] = path.join(self.paths["analysis_path"], middlepath_dataset)
        self.paths["exec_env_data_path"] = path.join(self.paths["exec_env_path"], middlepath_dataset, datafoldername)

        # Set datapath (depends on data_root_path)
        if not path.basename(path.normpath(self.paths["data_root_path"])) == "<default>": #if path not ends with <default>
            self.paths.update({"datapath": self.paths["data_root_path"]}) # -> datapath is same as data_root_path, no default folder structure
        else: 
            self.paths["data_root_path"] = path.dirname(path.normpath(self.paths["data_root_path"])) # remove "<default>" 
            self.paths["datapath"] = path.join(self.paths["data_root_path"], middlepath_dataset, datafoldername) # add MBEN_T/v00/analysis/human/sn/all/data
        
        self.paths["metadata_orig_path"] = path.join(self.paths["analysis_path"], middlepath, "metadata")
        self.paths.update({ "metadata_orig_filepath": path.join(self.paths["metadata_orig_path"],"metadata.tsv")})  
        
        # read, merge and write meta  
        self.paths["metadata_subj_path"] = path.join(self.paths["metadata_orig_path"], "metadata_subj.xlsx")
        self.paths["metadata_sample_path"] = path.join(self.paths["metadata_orig_path"], "metadata_sample.xlsx")
        meta_subj, meta_sample = '', ''
        if path.exists(self.paths["metadata_subj_path"]):
            meta_subj = pd.read_excel(self.paths["metadata_subj_path"])
        if path.exists(self.paths["metadata_sample_path"]):
            meta_sample = pd.read_excel(self.paths["metadata_sample_path"])
        if (len(meta_subj) > 0) & (len(meta_sample) > 0):
            meta = pd.merge(meta_sample, meta_subj, how="left", on = "subjID")
            meta.to_csv(self.paths["metadata_orig_filepath"], sep = "\t", index = False)
    
        self.paths["data_root_path"] = self.paths["datapath"]

        # Set data_root_filename (either given one or default)
        if "data_root_filename" not in self.paths.keys(): 
            self.paths["data_root_filename"] = self.name + ".h5ad"    

        # Set datafilepath_tmp based on use_pickle_data
        fileextension=""
        if bool(self.analysis_params['use_pickle_data']) == True:
            fileextension = ".pickle"
        else: 
            fileextension = ".h5ad"

        self.paths.update({
            "datasetpath":        middlepath_dataset, # path from projectname to datasetname: MBEN_T/v00/analysis/human/sn/all
            "datafilepath":       path.join(self.paths["datapath"], self.paths["data_root_filename"]),
            "metadatapath":       path.join(self.paths["datapath_tmp"], "metadata"),
            "priorknowledge":     path.join(self.paths["analysis_path"], middlepath, "priorKnowledge"), # TODO: change this to storage path 
            "priorknowledge_tmp": path.join(self.paths["analysis_path"], middlepath, "priorKnowledge") 
        })
        self.paths["datafilepath_tmp"] = path.join(self.paths["datapath_tmp"], f"{pathlib.Path(self.paths['datafilepath']).stem}{fileextension}")
        self.data = "" # is set in __init__ of analysis obj
        super().__init__()

    def print_info(self):
        display(Markdown(f'### Dataset {self.name}  '))
        display(Markdown('**Analysis Parameters**  '))
        print(json.dumps(self.analysis_params, indent=4, sort_keys=True, default=str))
        display(Markdown('**Paths**  '))
        print(print_paths(self.paths))

    def get_paths(self) -> dict:
        return self.paths

    def read_data(self, filepath):
        fileextension = path.splitext(filepath)[1]
        if  fileextension == '.pickle':
            with open(filepath, "rb") as f:
                self.data = dill.load(f)
            print("Data was read in from pickle file.")
        elif fileextension == '.h5ad':
            self.data = sc.read(filepath, cache = True)
        elif fileextension == '.tsv':
            self.data = pd.read_csv(filepath, sep='\t')
            # set row names
            gene_id = self.analysis_params['preprocessing']['gene_id']
            self.data.index = self.data[gene_id]
            # create table for gene name translations
            if self.data.dtypes[0] == object and self.data.dtypes[1] == object:
                self.gene_translations = self.data.iloc[:, 0:2]
                print('The attribute "gene_translations" was added.')
                # remove non count columns (first two)
                #self.data = self.data.iloc[: , 2:]
            self.data = self.data.select_dtypes(exclude=['object']) 
            # if cols are samples
            if self.data.shape[0] > self.data.shape[1]: 
                self.data = self.data.transpose()     
            # to int 
            self.data = self.data.astype(int)         # done after the formatting so that there are no string columns anymore. 

            # to anndata
            from anndata import AnnData
            import numpy as np
            self.data = AnnData(self.data, dtype=np.float32)
            self.data.var_names_make_unique()
            # test
            #(data.var.index != self.gene_translations.index).sum()==0
        else: 
            return False
        return True

    def save_data(self, filepath):
        fileextension = path.splitext(filepath)[1]
        if  fileextension == '.pickle':
            with open(filepath, "wb") as dill_file:
                dill.dump(self.data, dill_file)
            print("Data was saved as pickle file.")
        elif fileextension == '.h5ad':
            sc.write(filepath, self.data)
            print("Data was saved as h5ad file.")
    
    def set_gene_names(self, ensembldata, gene_symbol = 'gene_name'):
        """Replaces var_names with gene_names from Ensembl
        
        Pseudocode: 
            Save current index as 'gene_id' column
            Translate index and save result in gene_symbol column
                If not in release, put gene_id on exclude list
            Remove exclude genes based on index
            Replace empty gene_symbols with gene_id
            Set gene_symbols as index, convert to str and make unique

        Args:
            ensembldata (EnsemblRelease): from pyensembl import EnsemblRelease. Install the needed release like this !pyensembl install --release 109 --species mouse
            gene_symbol (str, optional):  Defaults to 'gene_name'.
        """
        exclude = []
        self.data.var['gene_id'] = self.data.var.index
        for id in self.data.var.index:
            try:
                self.data.var.loc[id, gene_symbol] = ensembldata.gene_by_id(id).gene_name
            except(ValueError): 
                exclude += [id]
        self.data = self.data[:, ~self.data.var_names.isin(exclude)]
        print(f'Dropped the following genes because they are not part of the chosen release: {exclude}')
        
        # replace empty gene_names with gene_id. If this is not done and the gene_names are used for the index, the empty names are replaced with negative integers.
        for i in range(len(self.data.var[gene_symbol])):
            if self.data.var[gene_symbol][i] =='':
                self.data.var[gene_symbol][i] = self.data.var['gene_id'][i]
        self.data.var.index = self.data.var[gene_symbol]
        self.data.var_names = self.data.var_names.astype(str)
        self.data.var_names_make_unique()


########################
#### Analysis Class ####
########################  
class Analysis(AnalysisI): 
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
            params_path (str): (relative) path to analysis params
        """
        # Read analysis_params.yaml
        with open(path.join(params_path, "analysis_params.yaml")) as stream:
            self.analysis_params = yaml.load(stream, Loader=yaml.FullLoader)

        self.analysis_params = replace_dictvalues(self.analysis_params)

        # Init paths from the information in analysis_params
        self.__paths = deepcopy(self.analysis_params["proj_params"]["paths"])
        # Create datasets
        self.datasets = [constructor(name, seq_type, organism, self.analysis_params, deepcopy(self.__paths)) for name, seq_type, organism, constructor in datasets]        
        self.init_datasets()

    def get_paths(self) -> dict:
        return self.__paths

    def save_paths(self):
        """ Saves analysis params together with paths to yaml file. """
        self.analysis_params['proj_params']['paths'] = self.get_paths()
        with open(path.join(self.get_paths()['analysis_path'], self.analysis_params["proj_params"]["proj_id"], self.analysis_params["proj_params"]["version"], 'analysis_params_paths.yaml'), 'w+') as file:
            yaml.dump(self.analysis_params, file)
        for dataset in self.datasets:
            dataset.analysis_params['paths'] = dataset.get_paths()
            with open(path.join(dataset.paths['analysis_path'], dataset.paths['datasetpath'], 'analysis_params.yaml'), 'w+') as file:
                yaml.dump(dataset.analysis_params, file)

    def print_info(self):
        """ Prints basic information about the analysis. """
        display(Markdown('### Analysis Object '))
        display(Markdown('**Analysis Parameters**  '))
        print(json.dumps(self.analysis_params['proj_params'], indent=4, sort_keys=True, default=str))
        display(Markdown('**Paths**  '))
        print(json.dumps(self.get_paths(), indent=4, sort_keys=True, default=str))
        for dataset in self.datasets:
            dataset.print_info()

        
    def init_datasets (self) -> None :
        """ Reads and cleans datasets """
        @loop(self.datasets, True)
        def init(self:Analysis) -> None:
            """ Read hf5ad data if no pickle exists. Reads into 'data' property of dataset. 

            ### Pseudocode
            ----------
            read datafilepath_tmp (either pickle or h5ad file)
                create datapath_tmp if not exists
                if datapath != ''
                    read datafilepath
                    if error:
                        if datapath not existing
                            assume new version, copy from mounted path of old version to new version, then read
                        else stop without reading data
                    save data according to file extension of datafilepath_tmp (pickle or h5ad)
                else: no data because there was no data_root_path
            """
            # Create tmp variables
            datapath = self.paths["datapath"] 
            datapath_tmp = self.paths["datapath_tmp"]
            datafilepath_tmp = self.paths["datafilepath_tmp"]
            datafilepath = self.paths["datafilepath"]

            # Read data
            if(path.exists(datafilepath_tmp)):
                if not self.read_data(datafilepath_tmp):
                    print(f"Please make sure that datafilepath_tmp ({datafilepath_tmp}) either ends with '.pickle' or with '.h5ad'.")
            else:
                if not path.exists(datapath_tmp):
                    makedirs(datapath_tmp)
                if datapath != '':
                    try: 
                        self.data = sc.read(datafilepath, cache = True)
                    except OSError as e:
                        # probably a new version number
                        if not path.exists(datapath):
                            # get previous version number
                            import re
                            datapath = str(datapath)
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
                            return
                    if not self.save_data(datafilepath_tmp):
                        print(f"Please make sure that datafilepath_tmp ({datafilepath_tmp}) either ends with '.pickle' or with '.h5ad'.")
                    print("Data was read in from datapath and is now saved in datapath_tmp.")
                else: 
                    self.data = []
                    print("Please be aware that no data was read in as no data_root_path was provided.")
        init()
        self.clean_datasets()
        self.add_metadata()

    ### Processing functions ###       
    def clean_datasets(self) :
        """ correct datatypes """
        @loop(self.datasets, True)
        def clean (dataset) : 
            if len(dataset.data) >= 1:
                if "seurat_clusters" in dataset.data.obs.columns : dataset.data.obs["seurat_clusters"] = dataset.data.obs.seurat_clusters.astype("category")
                            # to be sure that the sources and targets from the model match with the names used in the data, convert to lowercase
                #dataset.data.var_names = [x.lower() for x in list(dataset.data.var_names)]
                #if dataset.data.raw != None:
                #    dataset.data.raw.var.index = [x.lower() for x in list(dataset.data.raw.var.index)]
        clean()

    def add_metadata(self):
        """read metadata and save it in data.obs
        """
        @loop(self.datasets, True)
        def add_meta (dataset) : 
            if len(dataset.data) >= 1:
                if len(dataset.data.obs.columns) <=2: # no metadata loaded, yet
                    if path.exists(dataset.paths['metadata_orig_filepath']) and len(dataset.data) >= 1:
                        metadata = pd.read_csv(dataset.paths['metadata_orig_filepath'], sep='\t')
                        # make sampleID same as sampleID in counts table
                        metadata['sampleID'] = metadata['sampleID'].astype('string') 

                        # prepare obs
                        obs = pd.DataFrame(dataset.data.obs, columns = ['sampleID']) # data.obs has only an index so far
                        isbulk = False
                        if 'sampleID' not in dataset.data.obs.columns:
                            obs['sampleID'] = obs.index
                            obs = obs.reset_index(drop=True)
                            isbulk = True
                        # add metadata
                        
                        if isbulk == True:
                            dataset.data.obs = obs.merge(metadata, on='sampleID', how='left')
                            dataset.data.obs.index = dataset.data.obs['sampleID']
                            # test
                            dataset.data.obs.shape == metadata.shape
                        else:
                            obs['cellID'] = obs.index
                            dataset.data.obs = obs.merge(metadata, on='sampleID', how='left')
                            dataset.data.obs.index = obs['cellID']
                            dataset.data.obs = dataset.data.obs.drop('cellID', axis = 1)

                    else: 
                        print(f"Either no data was read in or no metadata file exists here:{dataset.paths['metadata_orig_filepath']}. You can use the add_metadata() function of your analysis object and the read_data() function of your dataset(s) to fix this.")
        add_meta()

    def get(self, dataset_name):
        for i in range(len(self.datasets)):
            if self.datasets[i].name == dataset_name:
                return self.datasets[i]

    ### Save & Load ###
    # These functions are for saving an intermediate status as pickle file. 
    def save_datasets(self):
        with open(path.join(self.paths["datapath"], "datasets.pickle"), "wb") as dill_file:
            dill.dump(self.datasets, dill_file)

    def load_datasets(self):
        return dill.load(path.join(self.paths["datapath"], "datasets.pickle"))


                     



