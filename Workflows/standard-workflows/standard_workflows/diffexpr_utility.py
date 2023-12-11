from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import path, makedirs, cpu_count
import decoupler as dc, pandas as pd
from IPython.display import display, Markdown 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import dill, itertools
from .analysis_baseclass import AnalysisI
from standard_workflows import analysis_baseclass as baseclasses
from standard_workflows import decoupler_utility as dcu
from standard_workflows import analysis_loops as al

def init_contrasts(dds, design_factor):
    """Formats contrasts for Deseq2.

    Args:
        dds (Deseq2DataSet): _description_
        design_factor (str): _description_

    Returns:
        list: contrasts that are usable with DeseqStats function
    """
    contrasts = list(itertools.combinations(dds.obs[design_factor].unique(), 2))
    def compare_strings(s):
        """ To be used as *key* of the funciton *sorted* to sort a list of string by last character."""
        return s[-1]
    # sort contrasts by last character to order like this: treated vs. untreated, for example ['m_1', 'm_0']
    contrasts = [sorted(list(cont), key=compare_strings, reverse = True) for cont in contrasts]
    # add the design_factor to each element
    return [[design_factor] + list(cont) for cont in contrasts]

def subset_contrasts(dds, contrasts):
    """Filter contrasts by the ones that are actually available.

    Args:
        dds (Deseq2DataSet): _description_
        contrasts (str): _description_

    Returns:
        list: _description_
    """
    contrasts_str = [f'{c0}_{c1}_vs_{c2}' for c0,c1,c2 in contrasts] # make strings
    flags = [e in dds.varm["LFC"].columns for e in contrasts_str] # compare with availabe LFC data
    return [x for x, flag in zip(contrasts, flags) if flag]

class DiffExpr(AnalysisI):
    """ This class provides methods to run DeSeq2 with bulk data. 
        DiffExpr works on the raw data in the AnnData object and saves the results under ap['rawType']. 
    """

    def __init__(self):
        """set paths
        """
        super().__init__()
        self.paths['diffExpr'] = {}
        input_type = self.analysis_params['rawType']
        self.paths['diffExpr']['deseq2_fig'] = path.join(self.paths['figpath'], input_type, 'diffExpr', 'deseq2')
        self.paths['diffExpr']['deseq2_res'] = path.join(self.paths['resultpath'], input_type, 'diffExpr', 'deseq2')
        self.ddss = {} # dict of deseq datasets

    # Build DESeq2 object
    # A DeseqDataSet fits dispersion and log-fold change (LFC) parameters from the data, and stores them.
    # X stores the count data,
    # obs stores design factors,
    # obsm stores sample-level data, such as "design_matrix" and "size_factors",
    # varm stores gene-level data, such as "dispersions" and "LFC".
    def run_deseq(self, design_factor, filtering, new = True, refit_cooks = True, min_replicates = 7):
        """Full differential expression workflow with Deseq2
        Creates filepath to dds file and either reads in existing pickle file or creates new one.

        Args:
            data (AnalysisI): The raw data is taken. Counts should be rather raw. It is advised to filter like this: At least 10 counts per gene per group.
            design_factors (List): If you want to test differences based on sex (female/male) and treatment (true/false) then it is preferable to combine the columns into one (sex_treatment) instead of giving both columns as design factor (=interaction term). 
            ref_level (List): The current pyDeseq implementation might not calculate all interactions when a column has more than two levels. To get the missing interactions, rerun with one of the levels that were not considered as reference level. For example given the condition sex_treatment and the reference level f_1, you get m0_f0, m1_f0, f1_f0
            new (bool, optional): If false, diffexpr_path_results is used to read in an existing deseq object. Defaults to True.
            filtering (String): for example 'filtab10pre5'
            refit_cooks (bool, optional): Defaults to True.
            min_replicates (int, optional): If there are less replicates, no cooks filtering is done. Defaults to 7.

        Returns:
            _type_: _description_
        """
        data = self.data.raw.to_adata().copy()
        paths = self.get_paths()['diffExpr']

        dds_dirpath = path.join(paths['deseq2_res'], filtering, design_factor[0], ''.join(['refLev', design_factor[1], '_cooks', str(int(refit_cooks)), '_minRep', str(min_replicates)]))
        dds_filepath = path.join(dds_dirpath, 'dds.pickle')
        if(path.exists(dds_filepath) and not new):
            with open(dds_filepath, "rb") as f:
                dds = dill.load(f)
        else:
            # remove samples that have NA for the condition
            #for design_factor in design_factors:
            #    is_na = data.obs[design_factor].isna()
            #    data.obs = data.obs.loc[~is_na]
            #    display(Markdown(f'The following samples have *NA* for the condition "{design_factor}"\n: {data.obs.loc[is_na]}'))

            # create Deseq dataset
            dds = DeseqDataSet(
                adata=data, 
                design_factors= design_factor[0],
                ref_level = design_factor,
                refit_cooks=refit_cooks,
                min_replicates = min_replicates, 
                n_cpus=cpu_count(),
            )

            # fit dispersions and LFCs
            dds.deseq2()
            display(Markdown('**LFC**  '))
            display(dds.varm["LFC"])

            # save
            if not path.exists(dds_dirpath):
                makedirs(dds_dirpath)
            with open(dds_filepath, "wb") as f:
                dill.dump(dds, f)

        dds_class = baseclasses.Analysis.new_dataset(Dds)
        dds_paths = deepcopy(self.paths)
        dds_paths['resultpath'] = dds_dirpath
        self.ddss[f'{design_factor[0]}_{design_factor[1]}'] = dds_class(dds, design_factor[0], design_factor[1], deepcopy(self.analysis_params), dds_paths)
        print(self)


    def get_contrasts(self, key, gene_symbols = 'gene_name', new = True):
        """Calculate contrasts

        Args:
            design_factor (_type_): _description_
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
        """
        dds = self.ddss[key]
        design_factor = dds.design_factor
        display(Markdown(f'Working on design factor {design_factor}'))
        contrasts = subset_contrasts(dds.data, init_contrasts(dds.data, design_factor))

        @al.loop (contrasts, True)
        def calc_contrasts(contrast, dds, new):
            # Save contrast
            contrast_str = f'{contrast[0]}_{contrast[1]}_vs_{contrast[2]}'
            dirpath = path.join(dds.paths['resultpath'], 'contrasts', contrast_str)
            file = 'contrast.csv'
            filepath = path.join(dirpath, file)
            if exists(filepath) and not new:
                results_df = pd.read_csv(filepath)
            else:
                display(Markdown(f'Contrast: {contrast}'))

                # Extract contrast 
                stat_res = DeseqStats(dds.data, contrast = contrast, n_cpus=8)#, contrast=[design_factor, 'm_0', 'f_1'], n_cpus=8)
                display(Markdown('**stat_res.summary**  '))
                # Compute Wald test to get pvalues
                # runs the whole statistical analysis, cooks filtering and multiple testing adjustement included
                # result of this is in stat_res.results_df
                stat_res.summary()
                d = pd.DataFrame(stat_res.results_df)
                d.insert(0, d.index.name, d.index)

                def save_tocsv(data, dirpath, file):
                    filepath = path.join(dirpath, file)
                    if not path.exists(dirpath):
                        makedirs(dirpath)
                    data.to_csv(filepath, index = False, header=True)
                save_tocsv(d, dirpath, file)

                # Shrink LFCs
                # Running lfc_shrink() will overwrite a DeseqStats’ log fold changes (and standard errors) with shrunk values. This can be checked using the shrunk_LFCs flag.
                # print(stat_res.shrunk_LFCs)  will be true when lfc_shrink was run
                stat_res.lfc_shrink(coeff=contrast_str)
                # Extract results
                results_df = stat_res.results_df
                # make index to column and remove index
                results_df[gene_symbols] = results_df.index
                results_df = results_df.reset_index(drop=True)
                #gene_transl = gene_transl.reset_index(drop=True)
                #save_tocsv(gene_transl, dirpath, 'gene_translation.csv')
                # merge with gene names/symbols
                #results_df = results_df.merge(gene_transl , on='gene_id', how='left')
                #col = results_df.pop('gene_id')
                #results_df.insert(0, 'gene_id', col)

                col = results_df.pop(gene_symbols)
                results_df.insert(0, gene_symbols, col)
                save_tocsv(results_df, dirpath, 'contrast_shrunk.csv')

                # for FUNKI
                genelist = results_df.loc[:,[gene_symbols, 'log2FoldChange']]
                save_tocsv(genelist, dirpath, 'genelist.csv')

            #display(Markdown('**Results**'))
            #display(results_df)
            dds.contrasts[contrast_str] = {'data': results_df, 'acts':[]}
            display(Markdown(f'**{contrast}**'))
            dds.plot_diffExpr(contrast_str, new = new)

        calc_contrasts(dds = dds, new = new)


class Dds(): 
    """ Dds is a Dataset that holds an Anndata object that stores differential expression results."""
    def __init__(self, data, design_factor, ref_level, ap, paths):
        """ Adds paths and reads in anndata object if it already exists or creates it newly via Decoupler. """
        super().__init__()
        self.data = data
        self.design_factor = design_factor
        self.ref_level = ref_level
        self.analysis_params = ap
        self.paths = paths
        #self.paths['resultpath'] = obj_path
        self.paths['figpath'] = self.paths['resultpath'].replace("results", "figures")
        self.contrasts = {}

    def plot_diffExpr(self, contrast_str, x_vals='log2FoldChange', y_vals=['padj', 'pvalue'], top=30, sign_thr=0.05, gene_symbols = 'gene_name', new = True):
        """ Volcano plots for a specific contrast

        Args:
            contrast (_type_): _description_
            x_vals (str, optional): _description_. Defaults to 'log2FoldChange'.
            y_vals (list, optional): _description_. Defaults to ['padj', 'pvalue'].
            top (int, optional): _description_. Defaults to 30.
            sign_thr (float, optional): _description_. Defaults to 0.05.
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
        """
        contrast = self.contrasts[contrast_str]['data']
        figpath = path.join(self.paths['figpath'], 'contrasts', contrast_str, 'volcano')
        contrast.index = contrast[gene_symbols]
        if not path.exists(figpath):
            makedirs(figpath)
        figpath = path.join(figpath, f'top{top}_sig{int(sign_thr*100)}')
        for y in y_vals:
            filepath = figpath + '_' + y + '.pdf'
            if exists(filepath) and not new:
                ...#read figures
            else:
                dc.plot_volcano_df(contrast, x=x_vals, y=y, top=top, sign_thr=sign_thr, save=filepath, return_fig=False)
       
