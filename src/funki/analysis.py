import scanpy as sc

def sc_trans_qc_metrics(data, var_name='mito'):
    '''
    Takes a single-cell transcriptomics data set and computes several quality
    control (QC) metrics according to the provided variable annotation.

    :param data: A single-cell transcriptomic data set instance
    :type data: :class:`funki.input.DataSet`
    :param var_name: The variable over which to compute the QC metrics. This
        name should exist in ``DataSet.var_names``, defaults to ``'mito'``
    :type var_name: str, optional

    :returns: The resulting filtered data set after applying the specified
        thresholds
    :rtype: :class:`funki.input.DataSet`
    '''

    aux = data.copy()
    sc.pp.calculate_qc_metrics(aux, qc_vars=[var_name], percent_top=None,
                               log1p=False, inplace=True)

    return aux

def sc_clustering(data, alg='leiden', resolution=1.0):
    '''
    Computes the clustering of the cells according to the selected algorithm and
    resolution parameter. You can plot the resulting clustering with your
    dimensionality reduction approach of choice (:func:`funki.plots.plot_pca`,
    :func:`funki.plots.plot_tsne` or :func:`funki.plots.plot_umap`) by passing
    the clustering algorithm name to the keyword argument ``color`` (e.g.:
    ``plot_tsne(mydata, color='leiden')``)

    :param data: A single-cell transcriptomic data set instance
    :type data: :class:`funki.input.DataSet`
    :param alg: The algorithm to use for computing the clusters, can be either
        ``'leiden'`` or ``'louvain'``, defaults to ``'leiden'``
    :type alg: str, optional
    :param resolution: The resolution of the clustering, a higher number
        generates more and smaller clusters, defaults to ``1.0``
    :type resolution: float, optional
    '''

    if not ('distances' in data.obsp and 'connectivities' in data.obsp):
        sc.pp.neighbors(data)

    if alg == 'leiden':
        sc.tl.leiden(data, resolution=resolution)

    elif alg == 'louvain':
        sc.tl.louvain(data, resolution=resolution)

    else:
        print('Algorithm not recognized, please use "leiden" or "louvain".')
