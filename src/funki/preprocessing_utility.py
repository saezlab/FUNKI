import scanpy as sc, numpy as np, decoupler as dc, matplotlib.pyplot as plt, pandas as pd
from .analysis_baseclass import AnalysisI
from . import utility_functions
import pathlib # create dirs
import os # cpu_count
import matplotlib

class Preprocessing(AnalysisI):
    """_summary_

    Args:
        AnalysisI (_type_): _description_
    """
    def __init__(self):
        super().__init__()

        #input_type = self.analysis_params['xType']
        #use_hvg = True
        #self.paths['figure_pca_path'] = pathlib.Path(self.paths['figpath'], input_type, f'pca_hvg{use_hvg}')
        #self.paths['result_pca_path'] = pathlib.Path(self.paths['resultpath'], input_type, f'pca_hvg{use_hvg}')

    def plot_obs(self):
        for obs in self.data.obs.columns:
            fig = plt.figure()
            plt.hist(self.data.obs[obs], bins=self.data.n_obs, color='orange')
            plt.xlabel(obs)
            plt.ylabel('N obs')
            fig

    def get_highest_expr_genes(self, topn = 30, input_type = 'counts', gene_symbols = 'gene_name'):
        dirpath = pathlib.Path(self.paths['figpath'], input_type, 'highest_expr_genes')
        dirpath.mkdir(parents=True, exist_ok=True)        
        sc.settings.figdir = dirpath
        sc.pl.highest_expr_genes(self.data, n_top=topn, gene_symbols = gene_symbols, save= ".pdf", show=False)

    def plot_violins(self, use_raw=True, input_type = 'counts'):
        """ Plots violins 
        based on ap['preprocessing']['qc_cols_obs'] 
        grouped by ap['diffExpr']['conditions']
        

        Args:
            use_raw (str): Defines if .X or .raw is used. Results are saved under ap[xType] or ap[rawType] accordingly.
            input_type(str): Only used if not use_raw 
        """
        if use_raw: 
            input_type = self.analysis_params['rawType']
        fig_violin_path = pathlib.Path(self.paths['figpath'], input_type, 'violin')
        fig_violin_path.mkdir(parents=True, exist_ok=True)        
        sc.settings.figdir = str(fig_violin_path)
        
        conditions = self.analysis_params['diffExpr']['conditions']

        def plot_violin(cond):
            self.data.obs[cond] = self.data.obs[cond].astype('category')
            fig_width= 2.5 + len(self.analysis_params['preprocessing']['qc_cols_obs'])
            with matplotlib.pyplot.rc_context({'figure.figsize': (fig_width, 3)}):
                sc.pl.violin(self.data, self.analysis_params['preprocessing']['qc_cols_obs'], jitter=1, multi_panel=True, groupby = cond, stripplot = True, save=f'_{cond}.pdf', rotation=45, use_raw=use_raw, show=False)
        # joblib needs more time 
        #import time
        #start_time = time.time()
        #joblib.Parallel(n_jobs=os.cpu_count())(joblib.delayed(plot_violin)(cond) for cond in conditions)
        #print("--- %s seconds ---" % (time.time() - start_time))
        #start_time = time.time()
        [plot_violin(cond) for cond in conditions]
        #print("--- %s seconds ---" % (time.time() - start_time))
        #print(f'violins are saved under {input_type}/violin')


    def plot_filter_expr(self, raw = True, input_type='counts'):
        """Code simplified from dc.plot_filter_by_expr()

        Args:
            raw (bool, optional): _description_. Defaults to True.
        """
        import seaborn as sns
        fig = None
        fig, ax = plt.subplots(1, 1)
        
        if raw == True:
            data = self.data.raw.X
            input_type = self.analysis_params['rawType']
            filtered = False
        else:
            data = self.data.X
            filtered = True
        dirpath = pathlib.Path(self.paths['figpath'], input_type)
        dirpath.mkdir(parents=True, exist_ok=True)        

        params = self.analysis_params['preprocessing']['basicFilt']
        min_count = params['min_count'] 
        min_total_count = params['min_total_count'] 
        min_sample_size = params['min_sample_size']
        min_sample_size -= 1e-14
        min_total_count -= 1e-14  
        lib_size    = np.sum(data, axis=1)
        median_lib_size = np.median(lib_size)
        cpm         = (data * 1e6)/ lib_size.reshape(-1, 1) 
        cpm_cutoff  = (min_count * 1e6) / median_lib_size 
        total_counts = np.sum(data, axis=0)
        sample_size = np.sum(cpm >= cpm_cutoff, axis=0)

        if not filtered:
            total_counts[total_counts < 1.] = np.nan  # handle 0
        
        # plot
        sns.set_palette('Paired')
        sns.histplot(x=np.log10(total_counts), y=sample_size, cbar=True,
                        cbar_kws=dict(shrink=.75, label='Number of features'),
                        discrete=(False, True), ax=ax, palette='inferno')
        ax.axhline(y=min_sample_size - 0.5, c='gray', ls='--')
        ax.set_xlabel('Log10 total sum of counts')
        ax.set_ylabel('Number of samples')
        if not filtered:
            ax.axvline(x=np.log10(min_total_count), c='gray', ls='--')

        fig.savefig(os.path.join(dirpath, 'gene_expr.pdf'), bbox_inches='tight')


    def get_hvgs(self, batch_key, input_type):
        """ Using a batch key results a simple batch correction where the hvgs are chosen per batch instead of overall.
        """
        sc.settings.figdir = pathlib.Path(self.paths['figpath'], input_type, 'highly_variable_genes')
        sc.settings.figdir.mkdir(parents=True, exist_ok=True)        
        conds = self.analysis_params['diffExpr']['conditions']
        for cond in conds: # compute, plot and save different versions/batch key results
            if len(set(self.data.obs[cond]))<=10:
                sc.pp.highly_variable_genes(self.data, batch_key=cond)
                sc.pl.highly_variable_genes(self.data, save=f'_{cond}.pdf')
        sc.pp.highly_variable_genes(self.data, batch_key=batch_key) # (re)compute with the chosen batch_key to have this as data in the AnnData object


    def plot_pcas(self, top_pca_genes, based_on_hvg, input_type, gene_symbols='gene_name', size = 1500, dimensions=[], new = True):
        """Plots pca for all genes and cols in ['preprocessing']['qc_cols_obs'] and ['diffExpr']['conditions'].

        Args:
            top_pca_genes (_type_): _description_
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
            size (int, optional): _description_. Defaults to 1500.
        """
        if new:
            # set/create paths
            figpath = pathlib.Path(self.paths['figpath'], input_type, f'pca_hvg{based_on_hvg}')
            figpath.mkdir(parents=True, exist_ok=True)

            figpath_genes = os.path.join(figpath, 'genes')
            figpath_metadata = os.path.join(figpath, 'metadata')

            import joblib
            # loop over genes
            for pc, genes in top_pca_genes.items():
                def plot_per_gene(gene):
                    sc.settings.figdir = figpath_genes
                    sc.pl.pca(self.data, color=gene, dimensions = [(pc-1,pc)], gene_symbols = gene_symbols, size = size, annotate_var_explained=True, save = f'_pc{pc}_{gene}.pdf', show = False)
                joblib.Parallel(n_jobs=os.cpu_count())(joblib.delayed(plot_per_gene)(gene) for gene in genes)
            # loop over metadata
            dimensions = [(0,1), (0,2), (1,2), (2,3), (1,3), (3,4)] + dimensions
            sc.settings.figdir = figpath_metadata
            for o in self.analysis_params['preprocessing']['qc_cols_obs'] + self.analysis_params['diffExpr']['conditions']:
                sc.pl.pca(self.data, color=o, dimensions=dimensions, save=f"_{o}.pdf", size = size, show = False)
                sc.pl.pca(self.data, color=o, dimensions=(1,2,3), projection='3d', save=f"_{o}3d.pdf", size = size, show = False)
        else:
            print("Step 'plot_pcas' was skipped.")


    def pc_loadings(self, based_on_hvg, input_type, top = 20, gene_symbols = 'gene_name', pc_max = 5, pca_key = 'X_pca')->dict:
        """Based on: self.data.varm['PCs']. 
        Plots heatmap for top genes per pc ordered by first condition in ['diffExpr']['conditions]
        Saves top genes per pc.

        Args:
            top (int, optional): _description_. Defaults to 20.
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
            use_raw (bool, optional): _description_. Defaults to False.
        """
        # set paths
        sc.settings.figdir = pathlib.Path(self.paths['figpath'], input_type, f'pca_hvg{based_on_hvg}')
        sc.settings.figdir.mkdir(parents=True, exist_ok=True)  
        sc.settings.writedir = pathlib.Path(self.paths['resultpath'], input_type, f'pca_hvg{based_on_hvg}')
        sc.settings.writedir.mkdir(parents=True, exist_ok=True)  

        groupby = self.analysis_params['diffExpr']['conditions'][0]
        if gene_symbols != None:
            genes = self.data.var[gene_symbols]
        else:
            genes = self.data.var_names
        top_pca_genes = {}
        for pc in range(1,pc_max):
            # sort PC data (the pcs containing the loadings.)
            pc_sorted = np.argsort(self.data.varm['PCs'][:,pc-1]) 
            top_genes = genes[
                    np.concatenate((pc_sorted[:top],pc_sorted[-top:])).tolist()
                ].tolist()
            # order by position on that pc
            tempdata = self.data[np.argsort(
                        self.data.obsm[pca_key][:,pc-1] # X_pca: pca representation of data
                ),]
            sc.pl.heatmap(tempdata, var_names = top_genes, groupby=groupby, save=f'_pc{pc}_{groupby}.pdf', swap_axes = True, use_raw=False, gene_symbols=gene_symbols, show=False)
            top_pca_genes[pc] = top_genes
            pd.DataFrame(top_genes).to_csv(os.path.join(sc.settings.writedir,f'pca_topGenes_pc{pc}_{groupby}.csv'))
        return top_pca_genes
    
    def get_pcas(self, input_type, gene_symbols = 'gene_name', dimensions = [], pc_max = 5, newplots = True, pca_key = 'X_pca'):
        """Calculates pca without using hvg. Saves loadings. Plots pca

        Args:
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
        """
        for use_hvg in [False, True]:
            sc.settings.figdir = pathlib.Path(self.paths['figpath'], input_type, f'pca_hvg{use_hvg}')
            sc.settings.figdir.mkdir(parents=True, exist_ok=True)  
            sc.settings.writedir = pathlib.Path(self.paths['resultpath'], input_type, f'pca_hvg{use_hvg}')
            sc.settings.writedir.mkdir(parents=True, exist_ok=True)  
            sc.tl.pca(self.data, svd_solver='arpack', random_state=0, use_highly_variable=use_hvg)
            # Get loadings for each gene for each PC
            if gene_symbols != None:
                index = self.data.var[gene_symbols]
            else:
                index = self.data.var_names
            df_loadings = pd.DataFrame(self.data.varm['PCs'], index=index)
            df_loadings.to_csv(os.path.join(sc.settings.writedir, 'loadings.csv'))
            top_pca_genes = self.pc_loadings(input_type=input_type, based_on_hvg = use_hvg, gene_symbols=gene_symbols, pc_max=pc_max, pca_key = pca_key)
            size = self.analysis_params['preprocessing']['basicFilt']['pca_dot_size']
            self.plot_pcas(top_pca_genes, based_on_hvg = use_hvg, input_type = input_type, gene_symbols=gene_symbols, size = size, dimensions = dimensions, new = newplots)
        self.paths['figure_pca_path'] = sc.settings.figdir
        self.paths['result_pca_path'] = sc.settings.writedir
        
    
    def get_pca_meta_associations(self, pca_key = 'X_pca', new = True):
        """Does not recalculate if plot is existing."""
        if new:
            keys = self.analysis_params['diffExpr']['conditions']
            keys = [key for key in keys if len(set(self.data.obs[key])) <= 10]
            dc.get_metadata_associations(
                self.data,
                obs_keys = keys, #metadata columns to associate to PCs
                obsm_key=pca_key,  # where the PCs are stored
                uns_key=f'{pca_key}_anova',  # where the results are stored
                inplace=True
            )
            self.data.uns[f'{pca_key}_anova'].to_csv(os.path.join(self.paths['result_pca_path'], '{pca_key}_anova_metadata_associations.csv'))

            plt.figure(figsize=(7,10))
            ax, legend_axes = dc.plot_associations(
                self.data,
                uns_key=f'{pca_key}_anova',  # summary statistics from the anova tests
                obsm_key=pca_key,  # where the PCs are stored
                stat_col='p_adj',  # which summary statistic to plot
                obs_annotation_cols = keys, # which sample annotations to plot
                titles=['Adjusted p-values from ANOVA', 'Principle component scores']
            )
            
            plt.savefig(os.path.join(self.paths['figure_pca_path'], f'{pca_key}_anova_metadata_associations.pdf'), bbox_inches='tight')
            plt.show()
        else:
            print("Step 'meta_associations' was skipped.")

        # def plot_associations(self, keys):
        #     """Adjusted from dc.plot_associations. 
        #           Problems: PCs on x axis start with PC0 instead of PC1. 
        #           Colours go from hot to cold instead of cold to hot.
        #           PCA plot is in the wrong order.

        #     Args:
        #         keys (_type_): _description_
        #     """
        #     pvals_df = self.data.uns['pca_anova'].pivot(index='variable', columns='factor', values='p_adj').dropna()
        #     neglog_p_df = pd.DataFrame(data=-np.log10(pvals_df.values),
        #                             index=list(pvals_df.index),
        #                             columns=list(pvals_df.columns))
        #     pcs_df = pd.DataFrame(data=self.data.obsm['X_pca'],
        #                   columns=[f'PC{ind}' for ind in range(self.data.obsm['X_pca'].shape[1])],
        #                   index=self.data.obs.index)#.reset_index().drop('index', axis=1)

        #     import marsilea as ma
        #     import marsilea.plotter as mp
        #     from marsilea.plotter import MarkerMesh

        #     h1 = ma.Heatmap(neglog_p_df, cmap="Reds", width=4, height=1, name="h1")
        #     h1.add_title(top="Adjusted p-values from ANOVA", align="center")
        #     h1.add_legends()
        #     h1.add_dendrogram("left")
        #     h1.add_right(mp.Labels(neglog_p_df.index, fontsize=10))

        #     h2 = ma.Heatmap(pcs_df.values, cmap="RdBu", width=0.4, height=4, name="h2")
        #     h2.add_legends()
        #     h2.add_dendrogram("left")
        #     h2.add_title(top="Principle components scores", align="center")

        #     for key in keys:
        #         key_cat = list(self.data.obs[key].values)
        #         key_colors = mp.Colors(key_cat, cmap="tab20", label=key)
        #         h2.add_right(key_colors, pad=0.1, size=0.2)
        #     h2.add_bottom(mp.Labels(pcs_df.columns, fontsize=10))

        #     c = (h1 / .2 / h2)
        #     c.add_legends(side="right",
        #                 stack_by='row', stack_size=2)
        #     c.render()
        #     c.save(os.path.join(self.paths['figure_pca_path'], 'anova_metadata_associations.pdf'))
        # plot_associations(self, keys)

    def preprocess(self, gene_symbols = 'gene_name', input_type='counts'):
        """
        1. Plots highest_expr_genes
        2. Filters genes and 'cells', adds filtering attributes to ap
        Adds mt to mark mitochondiral genes starting with mt-/MT- and 'noname' to mark mouse genes starting with Gm or ending with Rik.
        Args:
            gene_symbols (str, optional): _description_. Defaults to 'gene_name'.
            top (int, optional): _description_. Defaults to 30.
        """
        params = self.analysis_params['preprocessing']
        if len(set(self.data.var[gene_symbols])) != self.data.n_vars: # if not unique
            self.get_highest_expr_genes(gene_symbols=None)
        else:
            self.get_highest_expr_genes(gene_symbols=gene_symbols)

        # var: add 'mt'
        self.data.var['mt'] = self.data.var[gene_symbols].str.startswith('mt-') | self.data.var[gene_symbols].str.startswith('MT-')
        params['qc_cols_var'] = ['mt']
        sc.pp.calculate_qc_metrics(self.data, qc_vars=params['qc_cols_var'], percent_top=None, log1p=False, inplace=True, parallel=True, use_raw=False)

        # obs: add 'isHighMT' (binary column for MT depending on otsu threshold)
        otsu_threshold = utility_functions.get_otsu_threshold(np.array(self.data.obs['pct_counts_mt']*100, dtype='int'))
        self.data.obs['isHighMT'] = np.where(self.data.obs['total_counts_mt']>= otsu_threshold, 1, 0)
        params['qc_cols_obs'] += ['isHighMT']

        # var: add 'noname' (genes without a real gene_symbol/name)
        if self.organism == 'mouse': 
            self.data.var['noname'] = self.data.var[gene_symbols].str.startswith('Gm') | self.data.var[gene_symbols].str.endswith('Rik')
            params['qc_cols_var'] += ['noname']
        
        # var: add qc cols
        sc.pp.calculate_qc_metrics(self.data, qc_vars=params['qc_cols_var'], percent_top=None, log1p=False, inplace=True, parallel=True, use_raw=False)
        qc_cols = params['qc_cols_var'].copy()
        for var in qc_cols:
            params['qc_cols_obs'] += [f'total_counts_{var}'] + [f'pct_counts_{var}'] 

        # violins
        self.plot_violins(use_raw=False, input_type=input_type)

        sc.pl.scatter(self.data, x='total_counts', y='pct_counts_mt')
        sc.pl.scatter(self.data, x='total_counts', y='n_genes_by_counts')
        sc.pl.scatter(self.data, x='pct_counts_mt', y='n_genes_by_counts')

    def filter(self, prev, gene_symbols = 'gene_name', pca_dims = [], pc_max=5, newpcaplots = True, skipviolins = False):
        """Execute filtering + hvg + pca

        Args:
            prev (int): prevalence (large_n)
            gene_symbols (str, optional): Defaults to 'gene_name'.
            pca_dims (list, optional): PCA dimensions. Defaults to [].
            pc_max (int, optional): Maximal PCA dimension for plotting heatmaps. Defaults to 5.
        """
        min_prop = 0
        params = self.analysis_params['preprocessing']
        input_type = f"{self.analysis_params['xType']}_prev{prev}"

        if self.seq_type.lower() not in ['sc', 'scrna', 'sn', 'snrna', 'scrnaseq', 'snrnaseq']:
            # prevalence filtering
            dc.plot_filter_by_expr(self.data, min_count = params['basicFilt']['min_count'], min_total_count=params['basicFilt']['min_total_count'], large_n=prev, min_prop = min_prop)
            genes_to_keep = dc.filter_by_expr(self.data, min_count = params['basicFilt']['min_count'], min_total_count=params['basicFilt']['min_total_count'], large_n=prev, min_prop=min_prop)
            print(f'Number of genes after prevalence filtering: {len(genes_to_keep)}')
            self.data = self.data[:, genes_to_keep]
        else:
            print("No prevalence filtering is applied for single cell data. You can remove the parameters 'min_total_count' and 'min_count'.")

        # filter: cells, genes, mito
        # params['n_vars_diff'] shows number of filtered out features
        n_vars_pre = self.data.n_vars

        #if params['basicFilt']['max_genes'] == '':
        #        max_genes = max(self.data.obs['n_genes_by_counts'])
        #else: 
        #    max_genes = params['basicFilt']['max_genes']
        sc.pp.filter_cells(self.data, min_genes = params['basicFilt']['min_genes'])
        #sc.pp.filter_cells(self.data, max_genes = max_genes)
        sc.pp.filter_genes(self.data, min_cells = params['basicFilt']['min_cells'])
        self.data = self.data[self.data.obs.pct_counts_mt < 5, :]
        n_vars_diff = n_vars_pre - self.data.n_vars
        params['n_vars_diff'] = n_vars_diff
        print(f'Filtered out features after filter_cells, filter_genes and pct_counts_mt < 5:{n_vars_diff}')

        # prepare downsample
        if params['basicFilt']['downsample_to'] == '':
            subs_to = min(self.data.obs['total_counts'])
        else: 
            subs_to = params['basicFilt']['downsample_to']
        # norm (+downsample), log  
        sc.pp.normalize_total(self.data, target_sum = subs_to)
        sc.pp.log1p(self.data)
        print(f'After downsampling: {self.data}')

        # highly variable genes
        self.get_hvgs(self.analysis_params['diffExpr']['conditions'][0], input_type=input_type)

        self.data.layers['log'] = self.data.X

        # regress out and scale each gene to unit variance. Clip values exceeding standard deviation 10
        sc.pp.regress_out(self.data, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(self.data, max_value=10)
        if not skipviolins:
            self.plot_violins(use_raw=False, input_type = input_type)
        if len(set(self.data.var[gene_symbols])) != self.data.n_vars: # if not unique
            gene_symbols = None
        self.get_pcas(input_type=input_type, gene_symbols=gene_symbols, dimensions=pca_dims, pc_max=pc_max, newplots=newpcaplots)

        self.get_pca_meta_associations(new=newpcaplots)
        # save
        self.save_data(os.path.join(self.paths["datapath"], f'{input_type}.pickle'))
        