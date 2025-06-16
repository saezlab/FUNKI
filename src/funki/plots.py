import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import decoupler as dc

from .analysis import sc_trans_qc_metrics
from .preprocessing import sc_trans_filter
from .common import _colors, is_numeric


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

        colors = _colors['blue']

    # Categorical variable
    elif not is_numeric(color_vals):

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

    fig.tight_layout()

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

        colors = _colors['blue']

    # Categorical variable
    elif not is_numeric(color_vals):

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

    fig.tight_layout()

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

        colors = _colors['blue']

    # Categorical variable
    elif not is_numeric(color_vals):

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

    fig.tight_layout()

    return fig


def plot_highest_expr(data, top=10, ax=None):
    '''
    Generates a box plot of the top expressed genes (based on mean expression).

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param top: Number of top genes to represent, defaults to ``10``
    :type top: int, optional
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting box plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    mean = data.X.mean(axis=0)
    inds = np.argpartition(mean, -top)[-top:]
    inds = inds[np.argsort(mean[inds])]
    
    usegenes = data.var_names[inds].values[::-1]

    df = data.to_df().loc[:, usegenes]

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    ax.boxplot(
        df,
        patch_artist=True,
        boxprops={'color': _colors['blue'], 'facecolor': _colors['white']},
        capprops={'color': _colors['blue']},
        whiskerprops={'color': _colors['blue']},
        flierprops={'markeredgecolor': _colors['blue']},
        medianprops={'color': _colors['red']},
    )

    ax.set_title(f'Top {top} expressed genes')
    ax.set_xlabel('Genes')
    ax.set_ylabel('Expression')

    ax.set_xticks(range(1, top + 1))
    ax.set_xticklabels(usegenes, rotation=90)

    if return_fig:

        fig.tight_layout()

        return fig


def plot_n_genes(data, ax=None):
    '''
    Generates a violin plot displaying the number of genes by counts. This is,
    number of genes per cell that have non-zero counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting violin plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if 'n_genes_by_counts' not in data.obs.keys():

        data = sc_trans_qc_metrics(sc_trans_filter(data, mito_pct=100))

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    parts = ax.violinplot(
        data.obs['n_genes_by_counts'].values,
        widths=1,
        showmedians=True,
    )

    for p in parts['bodies']:

        p.set_facecolor(_colors['blue'])
        p.set_edgecolor(_colors['blue'])

    ax.set_title('Number of genes per cell')
    ax.set_ylabel('No. of genes')

    ax.xaxis.set_visible(False)

    if return_fig:

        fig.tight_layout()

        return fig


def plot_total_counts(data, ax=None):
    '''
    Generates a violin plot displaying the total counts per gene.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting violin plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if 'total_counts' not in data.var.keys():

        data = sc_trans_qc_metrics(sc_trans_filter(data, mito_pct=100))

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    parts = ax.violinplot(
        data.var['total_counts'].values,
        widths=1,
        showmedians=True,
    )

    for p in parts['bodies']:

        p.set_facecolor(_colors['blue'])
        p.set_edgecolor(_colors['blue'])

    ax.set_title('Counts per gene')
    ax.set_ylabel('Counts')

    ax.xaxis.set_visible(False)

    if return_fig:

        fig.tight_layout()

        return fig


def plot_pct_counts_mito(data, ax=None):
    '''
    Generates a violin plot displaying the percentage of mitochondrial genes.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting violin plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if 'pct_counts_mito' not in data.obs.keys():

        data = sc_trans_qc_metrics(sc_trans_filter(data, mito_pct=100))

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    parts = ax.violinplot(
        data.obs['pct_counts_mito'].values,
        widths=1,
        showmedians=True,
    )

    for p in parts['bodies']:

        p.set_facecolor(_colors['blue'])
        p.set_edgecolor(_colors['blue'])

    ax.set_title('Mitochondrial genes per cell')
    ax.set_ylabel('Pct. of mito. genes')

    ax.xaxis.set_visible(False)

    if return_fig:

        fig.tight_layout()

        return fig


def plot_counts_vs_pct_mito(data, ax=None):
    '''
    Generates a scatter plot displaying the percentage of mitochondrial genes
    versus total gene counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting scatter plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if (
        'pct_counts_mito' not in data.obs.keys()
        or 'total_counts' not in data.obs.keys()
    ):

        data = sc_trans_qc_metrics(sc_trans_filter(data, mito_pct=100))

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    ax.scatter(
        x=data.obs['total_counts'].values,
        y=data.obs['pct_counts_mito'].values,
        c=_colors['blue'],
    )
    ax.set_title('Total counts vs. pct. mitochondrial genes')
    ax.set_ylabel('Pct. of mito. genes')
    ax.set_xlabel('Counts')

    if return_fig:

        fig.tight_layout()

        return fig


def plot_counts_vs_n_genes(data, ax=None):
    '''
    Generates a scatter plot displaying the number of genes by counts versus
    total gene counts.

    :param data: The data set from which to generate the figure
    :type data: :class:`funki.input.DataSet`
    :param ax: Matplotlib Axes instance where to draw the plot. Defaults to
        ``None``, meaning a new figure and axes will be generated.
    :type ax: `matplotlib.axes.Axes`_

    :returns: The figure contataining the resulting scatter plot. If an axes is
        passed, nothing is returned.
    :rtype: `matplotlib.figure.Figure`_ | None

    .. _matplotlib.axes.Axes: https://matplotlib.org/stable/api/_as_gen/matplot\
        lib.axes.Axes.html#matplotlib.axes.Axes
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if (
        'n_genes_by_counts' not in data.obs.keys()
        or 'total_counts' not in data.obs.keys()
    ):

        data = sc_trans_qc_metrics(sc_trans_filter(data, mito_pct=100))

    if ax is None:

        fig, ax = plt.subplots()
        return_fig = True

    else:

        return_fig = False

    ax.scatter(
        x=data.obs['total_counts'].values,
        y=data.obs['n_genes_by_counts'].values,
        c=_colors['blue'],
    )
    ax.set_title('Total counts vs. number of genes')
    ax.set_ylabel('No. of genes')
    ax.set_xlabel('Counts')

    if return_fig:

        fig.tight_layout()

        return fig


def plot_dex(data, logfc_thr=1.0, fdr_thr=0.05, top=15):
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
    :param top: Number of top genes for which to display their gene name,
        defaults to ``15``.
    :type top: int, optional

    :returns: The figure contataining the resulting scatter plot
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    if any(x not in data.var_keys() for x in ['log2FoldChange', 'padj']):

        raise KeyError(
            'Results of differential expression not found in the DataSet '
            'provided, please run funki.analysis.diff_exp() first'
        )
    
    fig = dc.pl.volcano(
        data=data.var,
        x='log2FoldChange',
        y='padj',
        thr_stat=logfc_thr,
        thr_sign=fdr_thr,
        top=top,
        return_fig=True,
        color_neg=_colors['red'],
        color_null=_colors['gray'],
        color_pos=_colors['blue'],
    )

    return fig


def plot_enrich(
    data,
    top=10,
    method=None,
):
    '''
    Generates a dotplot displaying the top results of an enrichment analysis
    based on the provided methods.

    :param data: The data set from which to generate the figure (it is assumed
        that ``funki.analysis.enrich()`` as been performed beforehand).
    :type data: :class:`funki.input.DataSet`
    :param top: Number of top enriched gene sets to display based on their
        score, defaults to ``10``.
    :type top: int, optional
    :param method: Name of the enrichment method to plot the results from.

    :returns: The figure contataining the resulting bar plot
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/_as_gen/mat\
        plotlib.figure.Figure.html#matplotlib.figure.Figure
    '''

    try:

        # Raise error if method not present or no enrichment has ben run at all
        if not method in data.uns['funki']['enrich']['methods']:

            raise KeyError

    except KeyError as ke:

        raise ke(
            'Enrichment results for method not found in DataSet, please run ,'
            '`funki.analysis.enrich()` beforehand.'
        )

    score = data.obsm[f'score_{method}']
    pval = data.obsm[f'padj_{method}']

    res = (
        score
        .melt(value_name='score')
        .merge(
            pval
            .melt(value_name='pvalue')
            .assign(logpval=lambda x: x['pvalue'].clip(2.22e-4, 1))
            .assign(logpval=lambda x: np.log10(x['logpval']))
        )
    )

    fig, ax = plt.subplots()

    dc.pl.dotplot(
        df=res,
        x='score',
        y='variable',
        s='logpval',
        c='score',
        scale=1,
        top=top,
        ax=ax
    )
    ax.set_title(f'Enrichment results from method {method}')
    
    fig.tight_layout()

    return fig
