import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px

from .analysis import sc_trans_qc_metrics


def plot_pca(data, color=None, use_highly_variable=True, recalculate=False,
             **kwargs):
    '''
    Plots the dimensionality reduction PCA results of a data set.

    :param data: The data set from which to compute the PCA
    :type data: :class:`funki.input.DataSet`
    :param color: Variables or observations to color from, defaults to ``None``
    :type color: str | list[str], optional
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
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    .. _scanpy.pp.pca(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.pp.pca.html
    '''

    if recalculate:
        data._del_meta({'obsm': 'X_pca'})

    if use_highly_variable:
        sc.pp.highly_variable_genes(data, inplace=True)

    if 'X_pca' not in data.obsm:
        sc.pp.pca(data, use_highly_variable=use_highly_variable, **kwargs)

    colors = data.obs[color].values if color in data.obs_keys() else None
    
    fig = px.scatter(
        data.obsm.to_df(),
        x='X_pca1',
        y='X_pca2',
        color=colors,
        labels={
            'X_pca1': 'PC 1',
            'X_pca2': 'PC 2',
            'color': color
        },
        width=800,
        height=500
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)

    return fig

def plot_tsne(data, color=None, perplexity=30, recalculate=False):
    '''
    Plots the dimensionality reduction t-SNE results of a data set.

    :param data: The data set from which to compute the t-SNE
    :type data: :class:`funki.input.DataSet`
    :param color: Variables or observations to color from, defaults to ``None``
    :type color: str | list[str], optional
    :param perplexity: Perplexity hyperparmaeter for the t-SNE representation.
        Relates to the number of nearest neighbours, defaults to ``30``
    :type perplexity: int, optional
    :param recalculate: Whether to recalculate the dimensionality reduction,
        defaults to ``False``
    :type recalculate: bool, optional

    :returns: The figure contataining the scatter plot showing the tSNE
        embedding
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    if recalculate:
        data._del_meta({'obsm': ['X_pca', 'X_tsne'], 'uns': 'tsne'})

    if 'X_tsne' not in data.obsm:
        sc.tl.tsne(data, perplexity=perplexity)

    colors = data.obs[color].values if color in data.obs_keys() else None
    
    fig = px.scatter(
        data.obsm.to_df(),
        x='X_tsne1',
        y='X_tsne2',
        color=colors,
        labels={
            'X_tsne1': 'tSNE 1',
            'X_tsne2': 'tSNE 2',
            'color': color
        },
        width=800,
        height=500
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)

    return fig


def plot_umap(data, color=None, min_dist=0.5, spread=1.0, alpha=1.0, gamma=1.0,
              recalculate=False, **kwargs):
    '''
    Plots the dimensionality reduction UMAP results of a data set.

    :param data: The data set from which to compute the UMAP
    :type data: :class:`funki.input.DataSet`
    :param color: Variables or observations to color from, defaults to ``None``
    :type color: str | list[str], optional
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
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    .. _scanpy.tl.umap(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.tl.umap.html
    '''

    if recalculate:
        data._del_meta({'obsm': ['X_pca', 'X_umap'],
                        'obsp': ['distances', 'connectivities'],
                        'uns': ['umap', 'neighbors']})
    
    if not 'neighbors' in data.uns:
        sc.pp.neighbors(data)

    if 'X_umap' not in data.obsm:
        sc.tl.umap(data, min_dist=min_dist, spread=spread, alpha=alpha,
                   gamma=gamma, **kwargs)

    colors = data.obs[color].values if color in data.obs_keys() else None
    
    fig = px.scatter(
        data.obsm.to_df(),
        x='X_umap1',
        y='X_umap2',
        color=colors,
        labels={
            'X_umap1': 'UMAP 1',
            'X_umap2': 'UMAP 2',
            'color': color
        },
        width=800,
        height=500
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)

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
        points='all',
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
        points='all',
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
        points='all',
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
