import dill, re, collections
from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import path, makedirs
import scanpy as sc, decoupler as dc, pandas as pd
from .analysis_baseclass import AnalysisI
from . import baseclasses
from . import analysis_loops as al
import anndata
from functools import partial
from IPython.display import display, Markdown 
import seaborn as sns

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
        #self.paths.update({'actdir': path.join('log', self.paths['dc'])})
        self.paths.update({'pseudobulkdir': path.join(self.paths['figpath'], 'counts', self.paths['dc'], 'pseudobulk')})
        
    def get_all_acts(self, new = True):
        """ Create new Activity object for all parameter combinations. """
        pk = self.analysis_params['priorKnowledge']
        for modelcategory in pk: # pathways, tfs etc.
            for modeltype in pk[modelcategory]: # collectri, progeny, ...
                all_params = pk[modelcategory][modeltype] # dict of params for msigdb, otherwise not relevant
                for modelparams in pk[modelcategory][modeltype].values():
                    for min_n in self.analysis_params['decoupler']['min_n']:
                        for method in self.analysis_params['decoupler']['methods']:
                            if modeltype == 'MSigDB':
                                    model = eval('self.get_' + modeltype.lower() +'('+ '**all_params' + ')')
                                    print(f'Activity calculation starts for modeltype **{modeltype}** with parameter **{"_".join(all_params["substr"])}** and method **{method}**')
                                    self._get_acts(model, modeltype, '_'.join(all_params['substr']), method, min_n, deepcopy(self.paths), new)   
                            else:    
                                for param in modelparams: 
                                    model = self._getmodel(modeltype, param)
                                    print(f'Activity calculation starts for modeltype **{modeltype}** with parameter **{param}** and method **{method}**')
                                    self._get_acts(model, modeltype, param, method, min_n, deepcopy(self.paths), new)   
        
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
            if modeltype == 'liana':
                model = eval('self.get_' + modeltype.lower() + '(' + str(param) + ')')
            else:
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
    
    def get_liana(self, resources)->pd.DataFrame:
        """Prepare ligand receptor network. Take data from Liana and format it fitting for Decoupler.

        Args:
            resources (str): Liana resource

        Returns:
            pd.DataFrame: network that can be used as input for Decoupler
        """
        import liana as ln
        dirpath = path.join(self.get_paths()['priorknowledge'], 'liana')
        resource_names = '_'.join(resources)
        filepath = path.join(dirpath, f'{resource_names}.csv')
        if not path.exists(dirpath):
            makedirs(dirpath)
        if path.exists(filepath):
            net = pd.read_csv(filepath)
        else:
            nets = []
            for resource in resources:
                nets.append(ln.resource.select_resource(resource))
            net = pd.concat(nets)
            net = ln.resource.explode_complexes(net)
            # Create two new DataFrames, each containing one of the pairs of columns to be concatenated
            df1 = net[['interaction', 'ligand']]
            df2 = net[['interaction', 'receptor']]
            # Rename 
            df1.columns = ['source', 'target']
            df2.columns = ['source', 'target']
            # Concatenate 
            net = pd.concat([df1, df2], axis=0)
            net['weight'] = 1
            # Find and remove duplicated rows
            duplicates = net.duplicated()
            net = net[~duplicates]
            net.to_csv(filepath)
        return net

    def get_msigdb(self, files:list, substr = ['immune', 'inflammatory'])-> pd.DataFrame:
        """Creates a table (source, target) from multiple MSigDB gmt files for use with Decoupler.
        ToDo: if filepath exists read file
        ToDo: for each substr in substrings:...

        Args:
            files (list): gmt files from MSigDB
            substr (str, optional): Substring that must be contained in collection. Defaults to 'immune'.

        Returns:
            pd.DataFrame: source, target table for decoupler
        """
        dirpath = path.join(self.get_paths()['priorknowledge'], 'MSigDB')
        subsname = '_'.join(substr)
        filepath = path.join(dirpath, f'{subsname}.csv')
        if not path.exists(dirpath):
            makedirs(dirpath)
        net = pd.DataFrame()
        if path.exists(filepath):
            net = pd.read_csv(filepath)
        else:
            for filename in files:
                genesets_dict = {}
                netpath = path.join(dirpath, filename)
                with open(netpath) as genesets:
                    for line in genesets:
                        entries = line.strip().split("\t")
                        genesets_dict[entries[0]] = set(entries[2:]) 
                for key in genesets_dict.keys():
                    for subs in substr:
                        if subs in key.lower():
                            df = pd.DataFrame({'target': list(genesets_dict[key])})
                            df ['source'] = key
                            net = pd.concat([net, df], axis=0)
            net.to_csv(filepath)
        return net    

    def _get_acts(self, model, modeltype, modelparams, methods, min_n, paths, new):
        """ Create new Activity object and add it to dataset. """
  
        def create_actanndata(base, estimatekey, paths):
            data = dc.get_acts(base, obsm_key= estimatekey) # with estimates and drop estimates from obsm and save acts
            data.obsm.pop(estimatekey)
            # Write anndata
            data.write(paths['acts_data'])
            return data
        
        def calc_dc_acts(dc_input, use_raw, methods, method_pval, method_estimate, pvalkey, estimatekey, namedmethod, actpaths, is_consensus, weight):
            """Runs decoupler

            Args:
                use_raw (_type_): _description_
                methods (_type_): _description_
                method_pval (_type_): _description_
                method_estimate (_type_): _description_
                pvalkey (_type_): _description_
                estimatekey (_type_): _description_
                namedmethod (_type_): _description_
                actpaths (_type_): _description_
                is_consensus (bool): _description_
                weight (_type_): _description_

            Returns:
                _type_: _description_
            """
            if type(dc_input) == pd.DataFrame:
                input = dc_input
            else:
                input = dc_input.data
            result = dc.decouple(mat = input, net = model, methods=methods, min_n=min_n, consensus = is_consensus, weight = weight, use_raw = use_raw) 
            act_class = baseclasses.Analysis.new_dataset(Activity)   
            
            paths.update({'acts_data': path.join(actpaths['actdir'], 'data.h5ad')})
            paths.update({'mean_acts': path.join(actpaths['actdir'], 'mean_acts')}) 
            actpaths.update({'acts_data': path.join(actpaths['actdir'], 'data.h5ad')})
            actpaths.update({'mean_acts': path.join(actpaths['actdir'], 'mean_acts')}) 
            actpaths.update({'acts_estimate': path.join(actpaths['actdir'], 'estimate.csv')})
            actpaths.update({'acts_pvals': path.join(actpaths['actdir'], 'pvals.csv')})
            if not path.exists(actpaths['actdir']):
                makedirs(actpaths['actdir'])
            if result == None:
                # rename new obsm properties (consensus and all the rest)
                dc_input.data.obsm[pvalkey] = dc_input.data.obsm.pop(method_pval)
                dc_input.data.obsm[estimatekey] = dc_input.data.obsm.pop(method_estimate)

                # Create anndata obj with activity as .X
                data = create_actanndata(dc_input.data, estimatekey, actpaths)

                # save estimates and pvals separately, delete results from dataset obj
                pvals = dc_input.data.obsm.pop(pvalkey)
                estimates = dc_input.data.obsm.pop(estimatekey)

                # Add activity              
                dc_input.acts.append(act_class(data, modeltype, modelparams, namedmethod, deepcopy(actpaths)))
                pvals.to_csv(actpaths['acts_pvals'])
                estimates.to_csv(actpaths['acts_estimate']) 
                for cond in dc_input.analysis_params['diffExpr']['conditions']:
                    ranked = dc.rank_sources_groups(data, groupby=cond, reference='rest', method='t-test_overestim_var')
                    ranked.to_csv(path.join(dc_input.acts[-1].paths['actdir'], f'{cond}.csv'))
            else: 
                result[method_pval].to_csv(actpaths['acts_pvals'])
                result[method_estimate].to_csv(actpaths['acts_estimate']) 
                return act_class(result, modeltype, modelparams, namedmethod, deepcopy(actpaths))
            return self
        
        def calc_dds_acts(self, methods, method_pval, method_estimate, pvalkey, estimatekey, namedmethod, actpaths, is_consensus, weight):
            self.contr_acts=[]
            for contrast_name in self.contrasts.keys():
                contrast_data = self.contrasts[contrast_name]['data']
                actpaths['actdir'] = path.join(actpaths['actdir_part1'], 'contrasts', contrast_name, actpaths['actdir_part2'])
                dc_input = contrast_data[['stat']].T.rename(index={'stat': 'treatment.vs.control'})
                use_raw = False
                self.contrasts[contrast_name]['acts'].append(calc_dc_acts(dc_input, use_raw, methods, method_pval, method_estimate, pvalkey, estimatekey, namedmethod, actpaths, is_consensus, weight))

        def prepare_calc_acts(self, filt_type, calc_func, methods, new):
            # example: sn -> all_t -> raw -> decoupler -> Progeny -> 50_ulmmlm_estimate.pickle
            paths['actdir_part1'] = path.join(self.paths['resultpath'])
            paths['actdir_part2'] = path.join(filt_type, self.paths['dc'], modeltype)
            paths['actdir'] = path.join(self.paths['resultpath'], filt_type, self.paths['dc'], modeltype)
            dirpath = paths['actdir']
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

            namedmethods = list(zip(methodnames, methods))
            for namedmethod in namedmethods: 
                modelname = modeltype + paramname
                pvalkey = modelname + '_' + namedmethod[0] + '_pvals'
                estimatekey = modelname + '_' + namedmethod[0] + '_estimate'

                actpaths = deepcopy(paths)
                actpath = path.join(paths['actdir'], f'{paramname}_{namedmethod[0]}')    
                if not exists(actpath):
                    makedirs(actpath) 
                actpaths.update({'actdir': actpath})
                actpaths.update({'acts_data': path.join(actpath, 'data.h5ad')})
                actpaths.update({'acts_estimate': path.join(actpath, 'estimate.csv')})
                actpaths.update({'acts_pvals': path.join(actpath, 'pvals.csv')})

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
                    if 'weight' not in model.columns:
                        weight = None
                    else: 
                        weight = 'weight'
                    method_pval = namedmethod[1] + '_pvals'
                    method_estimate = namedmethod[1] + '_estimate'
                    calc_func(self, methods=methods, method_pval=method_pval, method_estimate=method_estimate, pvalkey=pvalkey, estimatekey=estimatekey, namedmethod=namedmethod, actpaths=actpaths, is_consensus=is_consensus, weight=weight)

        if hasattr(self, 'ddss'):
            for key in self.ddss.keys():
                calc_func = calc_dds_acts
                prepare_calc_acts(self.ddss[key], '', calc_func, methods, new)
                
        else:
            calc_func = partial(calc_dc_acts, use_raw = True)
            prepare_calc_acts(self, self.analysis_params['rawType'], calc_func, methods, new)
            

    def plot_acts_perDds(self, topn_targets, topn_associations):
        for dds_name, dds in self.ddss.items():
            display(Markdown(f'# DeSeq2 Dataset: {dds_name}'))
            for cont_name, contrast in dds.contrasts.items():
                display(Markdown(f'## Contrast: {cont_name}'))
                contrast_data = contrast['data']
                for act in contrast['acts']:
                    display(Markdown(f'### Activity'))
                    display(Markdown(f'**modeltype:** {act.modeltype}'))
                    display(Markdown(f'**modelparams:** {act.modelparams}'))
                    display(Markdown(f'**method(s):** {act.namedmethod}'))
                    method = act.namedmethod[0]
                    estimate = act.data[f'{method}_estimate']
                    pvals = act.data[f'{method}_pvals']
                    dirpath = act.paths['actdir'].replace('results', 'figures')
                    print(dirpath)
                    if not path.exists(dirpath):
                        makedirs(dirpath)
                    dc.plot_barplot(estimate, estimate.index[0], top=topn_associations, vertical=True, return_fig=False, save=path.join(dirpath,'barplot.svg'))
                    #plt.title(modeltype + ' estimate', fontsize=16)
                    #fig.savefig(act.paths['actdir']+'barplot.svg')
                    if act.modeltype != 'MSigDB':
                        if act.modeltype != 'ksn':
                            # top targets
                            targets_to_inspect = estimate.T.sort_values(by=estimate.T.columns[0], axis=0, ascending=False, key = abs)[:topn_targets].index
                            if type(act.modelparams) in [tuple, list]:
                                modelparams_name = str.lower('_'.join(act.modelparams))
                            else:
                                modelparams_name = act.modelparams
                            net = pd.read_csv(path.join(act.paths['priorknowledge'],act.modeltype, f'{modelparams_name}.csv'))
                            for source_name in targets_to_inspect:  #[:9]: #iloc[:,0][:5]:
                                targetfigurepath = path.join(dirpath,f'{source_name}.svg')
                                dc.plot_targets(contrast_data, stat=contrast_data.columns[1], source_name=source_name, net=net, top=topn_targets, return_fig = False, save = targetfigurepath)





# an activity is a Dataset that is created from another dataset and adds paths to activitiy files
# would make sense to do decoupler.consensus on read in estimate data instead of recalculating everyting -> ulmmlm_consensus, viperulm_consensus
class Activity(): 
    """ An Activity is a Dataset that holds an Anndata object that stores activity estimations in its X property. These estimations are derived from the tool Decoupler. """
    def __init__(self, data, modeltype, modelparams, namedmethod, paths):
        """ Adds paths and reads in anndata object. """
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

        @al.loop(self.mean_acts, True)

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

