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
import .sc_analysis_baseclass as sc_classes
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import dill
from IPython.display import display, Markdown 
import decoupler as dc

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
                n_cpus=None,
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

        dds_class = sc_classes.Analysis.new_dataset(Dds)
        self.ddss[design_factor[0]] = dds_class(dds, design_factor[0], design_factor[1], deepcopy(self.paths), dds_dirpath)
        print(self)


class Dds(): 
    """ Dds is a Dataset that holds an Anndata object that stores differential expression results."""
    def __init__(self, data, design_factor, ref_level, paths, obj_path):
        """ Adds paths and reads in anndata object if it already exists or creates it newly via Decoupler. """
        super().__init__()
        self.data = data
        self.design_factor = design_factor
        self.ref_level = ref_level
        self.paths = paths
        self.paths['resultpath'] = obj_path
        self.paths['figpath'] = self.paths['resultpath'].replace("results", "figures")
        self.results = pd.DataFrame()

    def plot_diffExpr(self, x_vals='log2FoldChange', y_vals=['padj', 'pvalue'], top=30, sign_thr=0.05, gene_symbols = 'gene_name'):
            figpath = path.join(self.paths['figpath'], 'volcano')
            self.results.index = self.results[gene_symbols]
            if not path.exists(figpath):
                    makedirs(figpath)
            figpath = path.join(figpath, f'top{top}_sig{int(sign_thr*100)}')
            for y in y_vals:
                filepath = figpath + '_' + y + '.pdf'
                dc.plot_volcano_df(self.results, x='log2FoldChange', y=y, top=top, sign_thr=sign_thr, save=filepath )
       
