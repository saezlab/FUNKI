import scanpy as sc

def plot_pca(data, color=None, use_highly_variable=True):
    '''
    Plots the PCA results of a data set.

    :param data: The data set from which to compute the PCA
    :type data: :class:`funki.input.DataSet`
    :param color: Variables or observations to color from, defaults to ``None``
    :type color: str | list[str], optional
    :param use_highly_variable: Whether to use highly variable genes only or all
        genes available, defaults to ``True``
    :type use_highly_variable: bool, optional
    :returns: The figure contataining the resulting PCA
    :rtype: `matplotlib.figure.Figure`_

    .. _matplotlib.figure.Figure: https://matplotlib.org/stable/api/figure_api.\
        html#matplotlib.figure.Figure
    '''

    if 'X_pca' not in data.obsm:
        sc.tl.pca(data)

    if use_highly_variable:
        sc.pp.highly_variable_genes(data, inplace=True)

    return sc.pl.pca(data, color=color)
