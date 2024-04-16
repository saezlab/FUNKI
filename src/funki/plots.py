import scanpy as sc

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
    :returns: The figure contataining the resulting PCA
    :rtype: `matplotlib.figure.Figure`_

    .. _scanpy.pp.pca(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.pp.pca.html
    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/figure_api.\
        html#matplotlib.figure.Figure
    '''

    if recalculate:
        data._del_meta({'obsm': 'X_pca'})

    if use_highly_variable:
        sc.pp.highly_variable_genes(data, inplace=True)

    if 'X_pca' not in data.obsm:
        sc.pp.pca(data, use_highly_variable=use_highly_variable, **kwargs)

    return sc.pl.pca(data, color=color)

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
    :returns: The figure contataining the resulting t-SNE
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/figure_api.\
        html#matplotlib.figure.Figure
    '''

    if recalculate:
        data._del_meta({'obsm': ['X_pca', 'X_tsne'], 'uns': 'tsne'})

    if 'X_tsne' not in data.obsm:
        sc.tl.tsne(data, perplexity=perplexity)

    return sc.pl.tsne(data, color=color)

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
        `scanpy.pp.umap()`_
    :type \*\*kwargs: optional
    :returns: The figure contataining the resulting UMAP
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/figure_api.\
        html#matplotlib.figure.Figure
    .. _scanpy.pp.umap(): https://scanpy.readthedocs.io/en/latest/generated/scan\
        py.pp.umap.html
    '''

    if recalculate:
        data._del_meta({'obsm': ['X_pca', 'X_umap'],
                        'obsp': ['distances', 'connectivities'],
                        'uns': ['umap', 'neighbors']})

    if not ('distances' in data.obsp and 'connectivities' in data.obsp):
        sc.pp.neighbors(data)

    if 'X_umap' not in data.obsm:
        sc.tl.umap(data, min_dist=min_dist, spread=spread, alpha=alpha,
                   gamma=gamma)

    return sc.pl.umap(data, color=color, **kwargs)
