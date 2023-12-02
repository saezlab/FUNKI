#from re import S
import sys, yaml, dill, re
#sys.path.insert(1, '/Users/hanna/Documents/projects/Workflows/Python/scUtilities/v4')
#from scUtilities import analysis_params
#import sc_analysis_loops as scl
from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import path, makedirs
import collections
import scanpy as sc, numpy as np, decoupler as dc, matplotlib.pyplot as plt, seaborn as sns, matplotlib as mpl, pandas as pd
from .sc_analysis_baseclass import AnalysisI
from .sc_analysis_baseclass import Baseanalysis
from . import sc_classes
from . import sc_analysis_loops as scl
import pandas as pd


class Liana(AnalysisI):
    """ This class is a wrapper around the tool Liana with some comfort functions. """
    def __init__(self):
        super().__init__()
        # add an object that holds the results data
        (self.paths).update({'liana_resdir': path.join(self.paths['resultpath'], 'log', 'liana')})
        self.paths.update({'liana_figdir': path.join(self.paths['figpath'], 'log', 'liana')})
        self.paths.update({'liana_aggdir': path.join(self.paths['liana_resdir'], 'aggregated')})

class Decoupler(AnalysisI):
    """ This class is a wrapper around the tool Decoupler with some comfort functions. """

    def __init__(self):
        """ A Decoupler object has a list of activity objects and additional paths. """
        super().__init__()
        self.acts = []
        self.paths.update({'dc': 'decoupler'})
        self.paths.update({'actsdir': path.join(self.paths['resultpath'], 'log', self.paths['dc'])})
        self.paths.update({'pseudobulkdir': path.join(self.paths['figpath'], 'counts', self.paths['dc'], 'pseudobulk')})
        # liana
        self.paths.update({'liana_resdir': path.join(self.paths['resultpath'], 'log', 'liana')})
        self.paths.update({'liana_figdir': path.join(self.paths['figpath'], 'log', 'liana')})
        self.paths.update({'liana_aggdir': path.join(self.paths['liana_resdir'], 'aggregated')})

    def get_all_acts(self, new = False):
        """ Create new Activity object for all parameter combinations. """
        pk = self.analysis_params['priorKnowledge']
        for modelcategory in pk: # pathways, tfs etc.
            for modeltype in pk[modelcategory]: # collectri, progeny, ...
                for modelparams in pk[modelcategory][modeltype].values():
                    for param in modelparams: 
                        for method in self.analysis_params['decoupler']['methods']:
                            model = self._getmodel(modeltype, param)
                            print(f'Activity calculation starts for modeltype **{modeltype}** with parameter **{param}** and method **{method}**')
                            self._get_acts(model, modeltype, param, method, deepcopy(self.paths), new)   

    # priorKnowledge
    def _getmodel(self, modeltype, param, filetype = 'csv'):
        """ Get prior knowledge model with the fitting get method of decoupler. """
        # pickle to  decoupler/priorKnowledge/progeny_topvalue
        # add to paths
        param_name = param
        if type(param_name) == tuple:
            param_name = str.lower(''.join(param_name))

        dirpath = path.join(self.paths['priorknowledge'], modeltype)
        filepath = path.join(dirpath, f"{param_name}.{filetype}")
        
        if(exists(filepath)):
            if filetype == 'pickle':
                with open(filepath, 'rb') as file:
                    model = dill.load(file)
                print('Prior Knowledge was read in from pickle files.')
            else:
                model = pd.read_csv(filepath)
                print('Prior Knowledge was read in from csv files.')
        else:
            organism = self.organism
            if str(param) == 'false':
                param = 'False'
            if str(param) == 'true':
                param = 'True'
            model = eval('dc.get_' + modeltype + '(organism,' + str(param) + ')')
            if not exists(dirpath):
                makedirs(dirpath)
            if filetype == 'pickle':
                with open(filepath, "wb") as dill_file:
                    dill.dump(model, dill_file)
            else:
                model.to_csv(filepath)
            print('Prior Knowledge was read in from Omnipath via Decoupler.')
        return model

    def get_msigdb(self, files:list, substr = 'immune')-> pd.DataFrame:
        """Creates a table (source, target) from multiple MSigDB gmt files for use with Decoupler.

        Args:
            files (list): gmt files from MSigDB
            substr (str, optional): Substring that must be contained in collection. Defaults to 'immune'.

        Returns:
            pd.DataFrame: source, target table for decoupler
        """
        dirpath = path.join(self.get_paths()['priorknowledge'], 'MSigDB')
        if not path.exists(dirpath):
            makedirs(dirpath)
        net = pd.DataFrame()
        for filename in files:
            genesets_dict = {}
            netpath = path.join(dirpath, filename)
            with open(netpath) as genesets:
                for line in genesets:
                    entries = line.strip().split("\t")
                    genesets_dict[entries[0]] = set(entries[2:]) 
            for key in genesets_dict.keys():
                if substr in key.lower():
                    df = pd.DataFrame({'target': list(genesets_dict[key])})
                    df ['source'] = key
                    net = pd.concat([net, df], axis=0)
        net.to_csv(path.join(dirpath, f'MSigDB_{substr}.csv'))
        return net    

    def _get_acts(self, model, modeltype, modelparams, methods, paths, new):
        """ Create new Activity object and add it to dataset. """
        # define activity class
        act = sc_classes.Analysis.new_dataset(Activity)

        # example: sn -> all_t -> raw -> decoupler -> Progeny -> 50_ulmmlm_estimate.pickle
        paths.update({'actsdir': path.join(paths['actsdir'], modeltype)})
        dirpath = self.paths['actsdir']
        if not exists(dirpath):
                makedirs(dirpath)

        # Define method names
        methods = list(methods) # Decoupler needs list as input, not tuple
        methodnames = deepcopy(methods)

        if(len(methods) >= 2):
            is_consensus = True
            # add name for consensus method; it is the concatenation of the methods
            cmethod = ''.join(methods)
            cmethod = cmethod + 'consensus'
            methodnames.insert(0, cmethod) # insert needs the names as a list, not tuple
            methods.insert(0, 'consensus')
        else: 
            is_consensus = False

        # Define model param names
        if type(modelparams) in [tuple, list]:
            paramname = str.lower(''.join(modelparams))
        else:
            paramname = str(modelparams)#[0])

        # Assumption: if consensus exists then other results exist as well. 
        # Assumption: Either one method or multiple with consensus
        # Caution: Consensus can get out of date when the results that it's based on get updated. 

        def create_actanndata(base, estimatekey, paths):
            data = dc.get_acts(base, obsm_key= estimatekey) # with estimates and drop estimates from obsm and save acts
            data.obsm.pop(estimatekey)
            # Write anndata
            data.write(paths['acts_data'])
            return data

        namedmethods = list(zip(methodnames, methods))
        for namedmethod in namedmethods: 
            modelname = modeltype + paramname
            pvalkey = modelname + '_' + namedmethod[0] + '_pvals'
            estimatekey = modelname + '_' + namedmethod[0] + '_estimate'

            actpaths = deepcopy(paths)
            actpath = path.join(paths['actsdir'], f'{paramname}_{namedmethod[0]}')    
            if not exists(actpath):
                makedirs(actpath) 
            paths.update({'acts_data': path.join(actpath, 'data.h5ad')})
            paths.update({'mean_acts': path.join(actpath, 'mean_acts')}) 
            actpaths.update({'acts_data': path.join(actpath, 'data.h5ad')})
            actpaths.update({'mean_acts': path.join(actpath, 'mean_acts')}) 

            actpaths.update({'acts_estimate': path.join(actpath, 'estimate.csv')})
            actpaths.update({'acts_pvals': path.join(actpath, 'pvals.csv')})
            
            actpaths.update({'actdir': actpath})

            if(exists(actpaths['acts_data']) and not new): 
                data = sc.read(actpaths['acts_data'], cache = True)
                print('Activitiy objects were read in from h5ad files.')
            elif(exists(actpaths['acts_estimate']) and not new): # Assumption: When there is an estimate file there is a pvals file, too.
                pvals = pd.read_csv(actpaths['acts_pvals'], index_col=0)
                estimates = pd.read_csv(actpaths['acts_estimate'], index_col=0)
                self.data.obsm[estimatekey] = estimates
                self.data.obsm[pvalkey] = pvals
                data = create_actanndata(self.data, estimatekey, actpaths)
            else:
                if(namedmethod[1] == 'consensus'): 
                    methods.pop(0) # delete the 'consensus' entry
                    new = False # When consensus is newly calculated 'new' is set to false so that the single methods that are part of the consensus calculation don't get recalculated. 

                # to be sure that the sources and targets from the model match with the names used in the data, convert to lowercase
                #model['target'] = [x.lower() for x in list(model['target'])]
                #model['source'] = [x.lower() for x in list(model['source'])]
                print(model)
                print(self.data)
                print(self.data.var)
                dc.decouple(mat = self.data, net = model, methods=methods, consensus = is_consensus, use_raw = True) 
                # rename new obsm properties (consensus and all the rest)
                self.data.obsm[pvalkey] = self.data.obsm.pop(namedmethod[1] + '_pvals')
                self.data.obsm[estimatekey] = self.data.obsm.pop(namedmethod[1] + '_estimate')

                # Create anndata obj with activity as .X
                data = create_actanndata(self.data, estimatekey, actpaths)

                # save estimates and pvals separately, delete results from dataset obj
                pvals = self.data.obsm.pop(pvalkey)
                estimates = self.data.obsm.pop(estimatekey)
                pvals.to_csv(actpaths['acts_pvals'])
                estimates.to_csv(actpaths['acts_estimate'])
                #estimate_cols = list(filter(re.compile(".*estimate").match, self.data.obsm_keys()))
                #pval_cols = list(filter(re.compile(".*pvals").match, self.data.obsm_keys()))
                #estimates = {k: dataset.data.obsm.pop(k) for k in estimate_cols}
                #pvals = {k: dataset.data.obsm.pop(k) for k in pval_cols}

            # Add activity                
            self.acts.append(act(data, modeltype, modelparams, namedmethod, deepcopy(actpaths)))
                


   


    
# an activity is a Dataset that is created from another dataset and adds paths to activitiy files
# would make sense to do decoupler.consensus on read in estimate data instead of recalculating everyting -> ulmmlm_consensus, viperulm_consensus
class Activity(): 
    """ An Activity is a Dataset that holds an Anndata object that stores activity estimations in its X property. These estimations are derived from the tool Decoupler. """
    def __init__(self, data, modeltype, modelparams, namedmethod, paths):
        """ Adds paths and reads in anndata object if it already exists or creates it newly via Decoupler. """
        super().__init__()
        self.data = data
        self.modeltype = modeltype
        self.modelparams = modelparams
        self.namedmethod = namedmethod
        self.paths = paths
        self.mean_acts = {}
        
    
    def get_mean_acts(self, minstd, groupby):
        meanactname = f'{groupby}_{minstd}'
        dirpath = self.paths['mean_acts']
        meanactpath = path.join(self.paths['mean_acts'], f'{meanactname}.csv')
        if(exists(meanactpath)):
            self.mean_acts.update({meanactname: pd.read_csv(meanactpath, index_col=0)})
        else:
            obsvalues = self.data.obs[groupby]
            counts = collections.Counter(obsvalues)
            if(len(counts) is not sum(counts.values())):
                if(len(counts) <= 20): 
                    if not exists(dirpath):
                        makedirs(dirpath)
                    dcresult = dc.summarize_acts(self.data, groupby=groupby, min_std=minstd)
                    if(not dcresult.empty):
                        self.mean_acts.update({meanactname: dcresult})
                        self.mean_acts[meanactname].to_csv(meanactpath)

    # The first column is called Unnamed: 0  because these are the row names that were saved beforehand
    # TODO: show == TRUE ??
    def plot_mean_acts(self):
        dirpath = self.paths['mean_acts']
        dirpath = str(dirpath).replace('results', 'figures')
        if not exists(dirpath):
            makedirs(dirpath)
        import matplotlib.pyplot as plt

        @scl.loop(self.mean_acts, True)
        def __plot_mean_acts(meanactname, self): 
            meanactpath = str(Path(dirpath,f'{meanactname}.pdf')) # str / str doesn't work, when using / one of the elements must be of type Path
            mean_acts = self.mean_acts[meanactname]

            # plot factor levels
            obsname = re.search('(.*)_.*', meanactname).group(1)
            obsvalues = self.data.obs[obsname]
            counts = collections.Counter(obsvalues)
            
            #if(len(counts) is not sum(counts.values())):
            #    plt.figure(figsize=(4,4))
            #    plt.xticks(rotation='vertical')
           #     plt.bar(counts.keys(), counts.values())
           #     plt.title(obsname) 
           #     plt.draw() 
           #     #plt.show() 
            if(len(counts) <= 20): 
                #print(counts)
                if(len(counts.keys()) >= 2): 
                    if(exists(meanactpath)):
                        import glob
                        import numpy
                        import matplotlib.pyplot

                        filenames = sorted(glob.glob(meanactpath))
                        s = ''
                        ##filenames = filenames[0:3]
                        for filename in filenames: 
                            s = s + str('\n###' + obsname + '  \n')
                            s = s + str('<p><img src="'+ filename+ '" style="width:500px;height:600px;"></p>')
                            ##data = numpy.loadtxt(fname=filename, delimiter=',')
                        #print(s)
                            ##from IPython.display import SVG, display
                        ##def show_svg():
                        ##    display(SVG(url='...'))

                            ##fig.tight_layout()
                            ##matplotlib.pyplot.show()
                    else:
                        ##from matplotlib import rcParams
                        ## USE FIGURE PATH
                        ##rcParams['figure.figsize'] = 15,20
                        ##fig.tick_params(labelsize=1)
                        sns.set(font_scale=0.5)
                        sns.set(rc={'figure.figsize':(15,20)})
                        if(len(mean_acts.columns) >= 2):
                            #.iloc[: , 1:]
                            fig = sns.clustermap(mean_acts, xticklabels=(mean_acts).columns, vmin=-2, vmax=2, cmap='coolwarm')
                            fig = fig.fig
                            if not len(mean_acts.columns) >= 50: 
                                fig.set_figwidth(len(mean_acts.columns))
                            else:
                                fig.set_figwidth(100)
                            fig.suptitle(obsname) #self.modeltype = modeltype
        #self.modelparams = modelparams
        #self.namedmethod = namedmethod
                            fig.savefig(meanactpath)
                        else: 
                            fig = sns.heatmap(mean_acts, xticklabels=(mean_acts).columns, vmin=-2, vmax=2, cmap='coolwarm')
                        
                        
                        ##fig.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 16)
                        # commented out because it doesn't work well together with a lot of datasets and multithreading
                        #if not len(mean_acts.columns) >= 50: 
                        #    plt.show()
        __plot_mean_acts(self)

