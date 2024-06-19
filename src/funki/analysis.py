import scanpy as sc
import anndata as ad
import decoupler as dc

from .input import DataSet

def enrich(data, net, methods=None, weight=None, **kwargs):
    '''
    Performs enrichment analysis using `Decoupler`_ based on a given network
    (e.g. gene set collection) and statistical method(s).

    :param data: The data set from which to perform the enrichment
    :type data: :class:`funki.input.DataSet`
    :param net: The network linking the features of the data to the attributes
        (e.g. pathways, gene sets, transcription factors, etc.)
    :type net: `pandas.DataFrame`_
    :param methods: Which statistical method(s) to use in order to compute the
        enrichment, defaults to ``None``. If none is provided, uses ``'mlm'``,
        ``'ulm'`` and ``'wsum'``. The option ``'all'`` performs all methods. To
        see all the available methods, you can run `decoupler.show_methods()`_
        function
    :type methods: NoneType | str | list[str]
    :param weight: Defines the column in the network containing the weights to
        use in the enrichment, defaults to ``None``.
    :type weight: NoneType | str
    :param \*\*kwargs: Other keyword arguments that passed to
        `decoupler.decouple()`_ function
    :type \*\*kwargs: optional

    :returns: ``None``, results are stored inplace of the passed ``data``
        object, which is a :class:`funki.input.DataSet` instance. Estimates,
        p-values and consensus scores (in case of multiple methods) are stored
        as part of the ``obsm`` attribute of the object.

    :rtype: NoneType

    .. _Decoupler: https://decoupler-py.readthedocs.io/en/latest/index.html
    .. _pandas.DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.D\
        ataFrame.html
    .. _decoupler.show_methods(): https://decoupler-py.readthedocs.io/en/latest\
        /generated/decoupler.show_methods.html#decoupler.show_methods
    .. _decoupler.decouple(): https://decoupler-py.readthedocs.io/en/latest/gen\
        erated/decoupler.decouple.html#decoupler.decouple
    '''

    # Making copy as AnnData to bypass Decoupler type checks
    aux = ad.AnnData(data.copy())

    dc.decouple(aux, net, methods=methods, use_raw=False, weight=weight,
                **kwargs)
    
    # Updating back the results to the original DataSet object
    data.obsm.update(aux.obsm)
 
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

    return DataSet(aux)

def sc_clustering(data, alg='leiden', resolution=1.0, neigh_kwargs={},
                  alg_kwargs={}):
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
    :param \*\*neigh_kwargs: Other keyword arguments that are passed to
        `scanpy.pp.neighbors()`_ function, defaults to ``{}``
    :type \*\*neigh_kwargs: dict, optional
    :param \*\*alg_kwargs: Other keyword arguments that are passed to the
        clustering algorithm `scanpy.tl.louvain()`_ or `scanpy.tl.leiden()`_
        functions, defaults to ``{}``
    :type \*\*alg_kwargs: dict, optional
    
    :returns: ``None``, results are stored inplace of the passed ``data`` object
    :rtype: NoneType

    .. _scanpy.pp.neighbors(): https://scanpy.readthedocs.io/en/latest/api/gene\
        rated/scanpy.pp.neighbors.html
    .. _scanpy.tl.louvain(): https://scanpy.readthedocs.io/en/latest/api/genera\
        ted/scanpy.tl.louvain.html
    .. _scanpy.tl.leiden(): https://scanpy.readthedocs.io/en/latest/api/generat\
        ed/scanpy.tl.leiden.html
    '''

    if not 'neighbors' in data.uns:
        sc.pp.neighbors(data, **neigh_kwargs)

    if alg == 'leiden':
        sc.tl.leiden(data, resolution=resolution, **alg_kwargs)

    elif alg == 'louvain':
        sc.tl.louvain(data, resolution=resolution, **alg_kwargs)

    else:
        print('Algorithm not recognized, please use "leiden" or "louvain".')
