import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

from .analysis import sc_trans_qc_metrics


def plot_pca(
        data,
        color=None,
        use_highly_variable=True,
        recalculate=False,
        **kwargs,
):
    '''
    Plots the dimensionality reduction PCA results of a data set.

    :param data: The data set from which to compute the PCA
    :type data: :class:`funki.input.DataSet`
    :param color: Variable to color from, defaults to ``None``
    :type color: str, optional
    :param use_highly_variable: Whether to use highly variable genes only or all
        genes available, defaults to ``True``
    :type use_highly_variable: bool, optional
    :param recalculate: Whether to recalculate the dimensionality reduction,
        defaults to ``False``
    :type recalculate: bool, optional
    :param \*\*kwargs: Other keyword arguments that can be passed to
        `scanpy.pp.pca()`_
    :type \*\*kwargs: optional

    :returns: The figure contataining the scatter plot showing the PCA embedding
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    .. _scanpy.pp.pca(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.pp.pca.html
    '''

    if recalculate:

        data._del_meta({'obsm': 'X_pca'})

    if use_highly_variable:

        sc.pp.highly_variable_genes(data, inplace=True)

    if 'X_pca' not in data.obsm:

        sc.pp.pca(data, use_highly_variable=use_highly_variable, **kwargs)

    color_vals = data.obs[color].values if color in data.obs_keys() else None
    
    if color_vals is None:

        colors = 'C0'

    # Categorical variable
    elif color_vals.dtype is np.dtype(object):

        cmap = {v: f'C{i % 10}' for i, v in enumerate(sorted(set(color_vals)))}
        colors = [cmap[c] for c in color_vals]

    # Numerical variable
    else:

        colors = color_vals

    fig, ax = plt.subplots()

    df = data.obsm.to_df()[['X_pca1', 'X_pca2']]

    im = ax.scatter(
        x=df.X_pca1.values,
        y=df.X_pca2.values,
        c=colors,
    )

    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')

    if color_vals.dtype is np.dtype(object):

        ax.legend(loc=0, handles=[
            Line2D([0], [0], label=k, marker='.', ms=10, mfc=v, mec=v, ls='')
            for k, v in cmap.items()
        ])

    else:

        fig.colorbar(im)

    return fig


def plot_tsne(
    data,
    color=None,
    perplexity=30,
    recalculate=False
):
    '''
    Plots the dimensionality reduction t-SNE results of a data set.

    :param data: The data set from which to compute the t-SNE
    :type data: :class:`funki.input.DataSet`
    :param color: Variable to color from, defaults to ``None``
    :type color: str, optional
    :param perplexity: Perplexity hyperparmaeter for the t-SNE representation.
        Relates to the number of nearest neighbours, defaults to ``30``
    :type perplexity: int, optional
    :param recalculate: Whether to recalculate the dimensionality reduction,
        defaults to ``False``
    :type recalculate: bool, optional

    :returns: The figure contataining the scatter plot showing the tSNE
        embedding
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if recalculate:

        data._del_meta({'obsm': ['X_pca', 'X_tsne'], 'uns': 'tsne'})

    if 'X_tsne' not in data.obsm:

        sc.tl.tsne(data, perplexity=perplexity)

    # TODO: Could probably move this to an external function?
    color_vals = data.obs[color].values if color in data.obs_keys() else None
    
    if color_vals is None:

        colors = 'C0'

    # Categorical variable
    elif color_vals.dtype is np.dtype(object):

        cmap = {v: f'C{i % 10}' for i, v in enumerate(sorted(set(color_vals)))}
        colors = [cmap[c] for c in color_vals]

    # Numerical variable
    else:

        colors = color_vals

    fig, ax = plt.subplots()

    df = data.obsm.to_df()[['X_tsne1', 'X_tsne2']]

    im = ax.scatter(
        x=df.X_tsne1.values,
        y=df.X_tsne2.values,
        c=colors,
    )

    ax.set_xlabel('tSNE 1')
    ax.set_ylabel('tSNE 2')

    if color_vals.dtype is np.dtype(object):

        ax.legend(loc=0, handles=[
            Line2D([0], [0], label=k, marker='.', ms=10, mfc=v, mec=v, ls='')
            for k, v in cmap.items()
        ])

    else:

        fig.colorbar(im)

    return fig


def plot_umap(
    data,
    color=None,
    min_dist=0.5,
    spread=1.0,
    alpha=1.0,
    gamma=1.0,
    recalculate=False,
    **kwargs
):
    '''
    Plots the dimensionality reduction UMAP results of a data set.

    :param data: The data set from which to compute the UMAP
    :type data: :class:`funki.input.DataSet`
    :param color: Variable to color from, defaults to ``None``
    :type color: str, optional
    :param min_dist: Effective minimum distance between the embedded points
    :type min_dist: float, optional
    :param spread: Effective scale of embedded points
    :type spread: float, optional
    :param alpha: Initial learning rate for the optimization
    :type alpha: float, optional
    :param gamma: Weighting applied to negative samples for the optimization
    :type gamma: float, optional
    :param recalculate: Whether to recalculate the dimensionality reduction,
        defaults to ``False``
    :type recalculate: bool, optional
    :param \*\*kwargs: Other keyword arguments that can be passed to
        `scanpy.tl.umap()`_
    :type \*\*kwargs: optional

    :returns: The figure contataining the scatter plot showing the UMAP
        embedding
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    .. _scanpy.tl.umap(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.tl.umap.html
    '''

    if recalculate:

        data._del_meta({
            'obsm': ['X_pca', 'X_umap'],
            'obsp': ['distances', 'connectivities'],
            'uns': ['umap', 'neighbors']
        })
    
    if not 'neighbors' in data.uns:

        sc.pp.neighbors(data)

    if 'X_umap' not in data.obsm:

        sc.tl.umap(data, min_dist=min_dist, spread=spread, alpha=alpha,
                   gamma=gamma, **kwargs)

    color_vals = data.obs[color].values if color in data.obs_keys() else None
    
    if color_vals is None:

        colors = 'C0'

    # Categorical variable
    elif color_vals.dtype is np.dtype(object):

        cmap = {v: f'C{i % 10}' for i, v in enumerate(sorted(set(color_vals)))}
        colors = [cmap[c] for c in color_vals]

    # Numerical variable
    else:

        colors = color_vals

    fig, ax = plt.subplots()

    df = data.obsm.to_df()[['X_uamp1', 'X_umap2']]

    im = ax.scatter(
        x=df.X_umap1.values,
        y=df.X_umap2.values,
        c=colors,
    )

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    if color_vals.dtype is np.dtype(object):

        ax.legend(loc=0, handles=[
            Line2D([0], [0], label=k, marker='.', ms=10, mfc=v, mec=v, ls='')
            for k, v in cmap.items()
        ])

    else:

        fig.colorbar(im)

    return fig


def plot_highest_expr(data, top=10):
    '''
    Generates a box plot of the top expressed genes (based on mean expression).

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param top: Number of top genes to represent, defaults to ``10``
    :type top: int, optional

    :returns: The figure contataining the resulting box plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    mean = data.X.mean(axis=0)
    inds = np.argpartition(mean, -top)[-top:]
    inds = inds[np.argsort(mean[inds])]
    
    usegenes = data.var_names[inds].values[::-1]

    fig = px.box(data.to_df(), y=usegenes, title=f'Top {top} expressed genes')
    fig.update_layout(xaxis_title='Genes', yaxis_title='Expression')
    
    return fig

def plot_n_genes(data):
    '''
    Generates a violin plot displaying the number of genes by counts. This is,
    number of genes per cell that have non-zero counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting violin plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if 'n_genes_by_counts' not in data.obs.keys():
        data = sc_trans_qc_metrics(data)

    fig = px.violin(
        data.obs['n_genes_by_counts'],
        y='n_genes_by_counts',
        title='Number of genes',
    )
    fig.update_layout(yaxis_title='Genes')

    return fig

def plot_total_counts(data):
    '''
    Generates a violin plot displaying the total gene counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting violin plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if 'total_counts' not in data.var.keys():
        data = sc_trans_qc_metrics(data)

    fig = px.violin(
        data.var['total_counts'],
        y='total_counts',
        title='Total counts'
    )
    fig.update_layout(yaxis_title='Counts')

    return fig

def plot_pct_counts_mito(data):
    '''
    Generates a violin plot displaying the percentage of mitochondrial genes.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting violin plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if 'pct_counts_mito' not in data.obs.keys():
        data = sc_trans_qc_metrics(data)

    fig = px.violin(
        data.obs['pct_counts_mito'],
        y='pct_counts_mito',
        title='Pct. of mitochondrial genes',
    )
    fig.update_layout(yaxis_title='Pct. mito. genes')

    return fig

def plot_counts_vs_pct_mito(data):
    '''
    Generates a scatter plot displaying the percentage of mitochondrial genes
    versus total gene counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting scatter plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if (
        'pct_counts_mito' not in data.obs.keys()
        or 'total_counts' not in data.obs.keys()
    ):
        data = sc_trans_qc_metrics(data)

    df = pd.DataFrame([data.obs['pct_counts_mito'], data.obs['total_counts']])

    fig = px.scatter(
        df.T,
        x='total_counts',
        y='pct_counts_mito',
        title='Total counts vs. pct. mitochondrial genes',
    )
    fig.update_layout(xaxis_title='Counts', yaxis_title='Pct. mito. genes')

    return fig

def plot_counts_vs_n_genes(data):
    '''
    Generates a scatter plot displaying the number of genes by counts versus
    total gene counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting scatter plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if (
        'n_genes_by_counts' not in data.obs.keys()
        or 'total_counts' not in data.obs.keys()
    ):
        data = sc_trans_qc_metrics(data)

    df = pd.DataFrame([data.obs['n_genes_by_counts'], data.obs['total_counts']])

    fig = px.scatter(
        df.T,
        x='total_counts',
        y='n_genes_by_counts',
        title='Total counts vs. number of genes',
    )
    fig.update_layout(xaxis_title='Counts', yaxis_title='Genes')

    return fig

def plot_dex(data, logfc_thr=1.0, fdr_thr=0.05):
    '''
    Plots the results of the differential expression analisis as a volcano plot.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param logfc_thr: Threshold for signifacnce based on the log2(FC) value,
        defaults to ``1.0``
    :type logfc_thr: float, optional
    :param fdr_thr: Threshold for signifacnce based on the FDR value, defaults
        to ``0.05``
    :type fdr_thr: float, optional

    :returns: The figure contataining the resulting scatter plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if any(x not in data.var_keys() for x in ['log2FoldChange', 'padj']):
        raise KeyError(
            'Results of differential expression not found in the DataSet '
            'provided, please run funki.analysis.diff_exp() first'
        )
    
    df = pd.DataFrame(
        [data.var['log2FoldChange'], -np.log10(data.var['padj'])]
    ).T
    df.columns = ['log2(FC)', '-log10(FDR)']
    df['sig'] = [
        ('UP' if r['log2FoldChange'] >= logfc_thr else 'DW')
        if r['padj'] <= fdr_thr and abs(r['log2FoldChange']) >= logfc_thr
        else 'NS'
        for i, r in data.var.iterrows()
    ]

    fig = px.scatter(
        df,
        x='log2(FC)',
        y='-log10(FDR)',
        color='sig',
        color_discrete_map={
            'NS': 'gray',
            'UP': 'red',
            'DW': 'blue'
        },
        title='Differential expression',
    )

    return fig

# TODO: implement for multiple methods (and significance?)
def plot_enrich(
    data,
    top=10,
    #methods=None
):
    '''
    Generates a horizontal barplot displaying the top results of an enrichment
    analysis based on the consensus score across methods.

    :param data: The data set from which to generate the figure (it is assumed
        that ``funki.analysis.enrich()`` as been performed beforehand).
    :type data: :class:`funki.input.DataSet`
    :param top: Number of top enriched gene sets to display based on their
        consensus score. If a negative number is provided, the bottom ones will
        be displayed instead.
    :type top: int

    :returns: The figure contataining the resulting bar plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    try:
        avail = data.uns['funki']['enrich']['methods']

    except KeyError as ke:
        raise ke(
            'Enrichment results not found in DataSet, please run ,'
            '`funki.analysis.enrich()` beforehand.'
        )

    # Ensuring list
#    methods = methods if type(methods) is list else [methods]
    # Ensuring methods are available
#    methods = [m for m in methods if m in avail]
    # If none available/provided default to all available ones
#    methods = methods or avail

    res = data.obsm['consensus_estimate'].mean(axis=0)
    res.sort_values(ascending=top < 0, inplace=True)
    res = res.head(abs(top))[::-1] if len(res) > abs(top) else res[::-1]

    fig = px.bar(
        res,
        orientation='h',
    )

    fig.update_layout(
        xaxis={'title': {'text': 'Consensus score'}},
        yaxis={'title': {'text': 'Gene set'}},
        showlegend=False,
    )

    return fig