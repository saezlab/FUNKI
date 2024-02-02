import functools, inspect
import scanpy as sc
from os.path import exists
from os import makedirs, path
from IPython.display import display, Markdown # to display Markdown in code chunk
from pathlib import Path
from copy import deepcopy
from multiprocessing import Pool
from toolz import compose
from pathos.threading import ThreadPool as Pool
import matplotlib.pyplot as plt
import numpy as np

class AnalysisI():
    def __init__(self) -> None:
        self.datasets = []

analysis = AnalysisI()

def loop(over, is_parallel):
    """ A looping decorator. """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            s = over
            if(type(over) != str):
                # attaches quotes around parameter  
                s = lambda_to_expr_str(lambda: over)
            if is_parallel == 'miau': 
                number_of_elements = len(eval(s))
                pool = Pool(nodes =  number_of_elements if number_of_elements <= 10 else 10)
                t = compose(
                    lambda x: func(x, *args, **kwargs)
                )
                list(pool.map(t, eval(s)))
            else:
                [func(elem, *args, **kwargs) for elem in eval(s)]
        return wrapper
    return decorator


def lambda_to_expr_str(lambda_fn):
    """c.f. https://stackoverflow.com/a/52615415/134077"""
    if not lambda_fn.__name__ == "<lambda>":
        raise ValueError('Tried to convert non-lambda expression to string')
    else:
        lambda_str = inspect.getsource(lambda_fn).strip()
        expression_start = lambda_str.index(':') + 1
        expression_str = lambda_str[expression_start:].strip()
        if expression_str.endswith(')') and '(' not in expression_str:
            # i.e. l = lambda_to_expr_str(lambda x: x + 1) => x + 1)
            expression_str = expression_str[:-1]
        return expression_str


# Calculate
# Decoupler: get_acts
# @CAUTION: When this is repeated the acts are added additionally to the previous acts to the act list of the dataset
@loop('analysis.datasets', True)
def get_acts(data):
    """ Calculate Decoupler.get_acts() """
    print('**', data.name, '**  ')
    data.get_all_acts(new = False)
    print('\n')


@loop('analysis.datasets', True)
def get_mean_acts(dataset):
    """ Loops over datasets and parameters to calculate the mean activities per activity. """
    print('\n**', dataset.name, '**  ')
    @loop(dataset.analysis_params['decoupler']['meanacts']['minstd'], False)
    @loop(dataset.acts, True)
    def get_peract_perminstd(acts, minstd):
        print('- *', acts.modeltype, '* -', sep='')
        @loop(acts.data.obs, False)
        def get_pergroupby(groupby):
            print(groupby)
            acts.get_mean_acts(float(minstd), groupby)
        get_pergroupby()
    get_peract_perminstd()


@loop('analysis.datasets', True)
def plot_mean_acts(dataset):
    """ Generates a barplot per obs and a heatmap per act. """
    print('### ', dataset.name)
    @loop(dataset.acts, True)
    def plot_peract(acts):
        print('Mean acts are plotted for ', acts.modeltype, acts.modelparams, '.  \n', sep='')

        #print(json.dumps((acts.mean_acts).keys(), indent=4, sort_keys=True, default=str))
        #print(acts.mean_acts).keys(), sep='')
        acts.plot_mean_acts()
        print('\n')
    plot_peract() 


def split_umap(adata, split_by, ncol=2, nrow=None, **kwargs):
    categories = adata.obs[split_by].cat.categories
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5*ncol, 4*nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        sc.pl.umap(adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs)
    plt.tight_layout()

@loop('analysis.datasets', False)
def plot_umap(data, groupby = None):
    display(Markdown(f'## Estimated activities for {data.name}'))
    @loop(data.acts, False)
    def _plot_umap(act):
        display(Markdown(f'{act.modeltype} with parameter {act.modelparams} and method {act.namedmethod}.'))
        dirpath = path.join(act.paths['actdir'], 'umaps')
        dirpath = str(dirpath).replace('results', 'figures')
        if not exists(dirpath):
            makedirs(dirpath)
        sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='pdf')
        sc.settings.figdir = dirpath
        sc.settings.verbosity = 3    
        
        if(len(act.data.var_names) <= 16):
            if not exists(f'{dirpath}/umap.pdf'):
                try:
                    sc.pl.umap(act.data, color=act.data.var_names, vcenter=0, cmap='coolwarm', save = True, show = True)
                except Exception as e: 
                    import sys
                    def eprint(*args, **kwargs):
                        print(*args, file=sys.stderr, **kwargs)
                    eprint('File should be saved here: ', dirpath, 'as umap.pdf\n', e)
                    raise
        else: 
            @loop(act.data.var_names, False) # don't write: 'act.data.var_names' -> act wouldn't be found anymore 
            def _plot_pervar(varname):
                if not exists(f'{dirpath}/umap_{varname}.pdf'):
                    try: 
                        if groupby == None:
                            sc.pl.umap(act.data, color=varname, vcenter=0, cmap='coolwarm', save = f'_{varname}.pdf', show = False)
                        else:
                            split_umap(act.data, color = varname, split_by=groupby,legend_loc = "right margin")
                    except Exception as e: 
                        import sys
                        def eprint(*args, **kwargs):
                            print(*args, file=sys.stderr, **kwargs)
                        eprint('File should be saved here: ', dirpath, 'as', varname, '\n', e)
                        raise
            _plot_pervar()
        
    _plot_umap()


## TODO: REDUNDANT TO UMAP :`((
@loop('analysis.datasets', False)
def plot_violin(data, groupby):
    display(Markdown(f'## Estimated activities for {data.name}'))
    @loop(data.acts, False)
    def _plot_violin(act):
        display(Markdown(f'{act.modeltype} with parameter {act.modelparams} and method {act.namedmethod}.'))
        dirpath = path.join(act.paths['actdir'], 'violins')
        dirpath = str(dirpath).replace('results', 'figures')
        dirpath = path.join(dirpath, groupby)
        if not exists(dirpath):
            makedirs(dirpath)
        sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='pdf')
        sc.settings.figdir = dirpath
        sc.settings.verbosity = 3    
        
        if(len(act.data.var_names) <= 16):
            if not exists(f'{dirpath}/violin.pdf'):
                try:
                    sc.pl.violin(act.data, keys=act.data.var_names, groupby=groupby, save = True, show = True)
                except Exception as e: 
                    import sys
                    def eprint(*args, **kwargs):
                        print(*args, file=sys.stderr, **kwargs)
                    eprint('File should be saved here: ', dirpath, 'as violin.pdf\n', e)
                    raise
        else: 
            @loop(act.data.var_names, False) # don't write: 'act.data.var_names' -> act wouldn't be found anymore 
            def _plot_pervar(varname):
                if not exists(f'{dirpath}/violin_{varname}.pdf'):
                    try: 
                        sc.pl.violin(act.data, keys=varname, groupby=groupby, save = f'_{varname}.pdf', show = False)
                    except Exception as e: 
                        import sys
                        def eprint(*args, **kwargs):
                            print(*args, file=sys.stderr, **kwargs)
                        eprint('File should be saved here: ', dirpath, 'as', varname, '\n', e)
                        raise
            _plot_pervar()
        
    _plot_violin()


@loop ('analysis.datasets', True)
def plot_meta_barplots(data, col1, col2):
    meta = data.data.obs
    if col2 != None: 
        ax = meta.groupby([col1, col2])[col2].count().unstack(0).plot.bar(title="Frequencies", figsize=(18,6))
        _ = ax.set_xlabel(col1)
        _ = ax.set_ylabel('Frequency')
    else: 
        d = meta.groupby(col1)[col1].count()
        if(len(d)<=30):
            ax =d.plot.bar(title="Frequencies", figsize=(18,6), color = 'orange')
            _ = ax.set_xlabel(col1)
            _ = ax.set_ylabel('Frequency')
        else: 
            ax = None
    import matplotlib.pyplot as plt 
    plt.show()


@loop ('analysis.datasets', True)
def get_subsets(ds, dataset_class) -> tuple:
    """ Writes subsets to local subset folder if non existant. 

    Returns
    -------
    Tuples with information needed for Dataset init. 
    """
    for metaCol in ds.analysis_params['subs']: 
        print('*The following subsets are calculated for dataset ', ds.name, ':*\n', ds.analysis_params['subs'][metaCol], sep='')
        subsConds = deepcopy(ds.analysis_params['subs'][metaCol]) # subsetConditions, 'cluster': ['clust', (1,2)]
        subsNamePrefix = subsConds[0] # 'clust'
        del subsConds[0] # [(1,2)]

        def subset(cond): # cond = (1,2)
            """Does not overwrite existing subsets"""
            # The name must contain the name of the dataset and the meta col as well so that it can be handled in the same list as the dataset that it's coming from.
            subsName = ds.name + '_' + subsNamePrefix + ''.join([str(e) for e in cond]) # all_t_clust12
            subs = ds.data[ds.data.obs[metaCol].isin(cond)]
            path = Path(ds.paths['subsetspath'], subsName, 'data') # MBEN/v01/human/sn/clust12/ 
            if not exists(path):
                makedirs(path)
            path = str(path) + '/' + subsName + '.h5ad' # MBEN/v01/human/sn/all_t/subsets/clust/all_t_clust12.h5ad
            if not exists(path):
                subs.write(path)
            return (subsName, ds.seq_type, ds.organism, dataset_class) # (all_t_clust12, 'sn', 'human', <<dc_dataset>>)
        
        [subset(cond) for cond in subsConds]



@loop ('analysis.datasets', True)
def prepare_nfcore(dataset, pipeline, new = True):
    """Prepare run, create sample sheet, get reference genome. 
    Initialise 'update_content' with the pipeline params from the ap file as new attribute of the dataset

    Args:
        dspos (int): dataset position in the dataset list of analysis obj.
        pipeline (str): name of pipeline to prepare for
    """
    dataset.prepare_run(pipeline)
    dataset.create_sample_sheet(pipeline, new)
    dataset.get_reference()
    dataset.nfc_custom_params = {}
    dataset.nfc_custom_params['custom_default'] = dataset.analysis_params['nfcore']['pipeline'][pipeline]['params']


