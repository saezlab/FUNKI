import scanpy as sc
import anndata as ad
import pandas as pd
import decoupler as dc
from decoupler.mt import _methods
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from .input import DataSet


def enrich(
    data,
    net,
    contrast=None,
    method=None,
    source='source',
    target='target',
    weight=None,
    **kwargs
):
    '''
    Performs enrichment analysis using `Decoupler`_ based on a given network
    (e.g. gene set collection) and statistical method(s).

    :param data: The data set from which to perform the enrichment
    :type data: :class:`funki.input.DataSet`
    :param net: The network linking the features of the data to the attributes
        (e.g. pathways, gene sets, transcription factors, etc.)
    :type net: `pandas.DataFrame`_
    :param contrast: Which result of the differential expression to use for the
        enrichment. Must be present in ``data.varm_keys`` named with the format
        ``'{contrast_var}_vs_{ref_var}'``. Defaults to ``None``.
    :type contrast: str
    :param method: Which statistical method to use in order to compute the
        enrichment, defaults to ``None``. If none is provided, uses ``'ulm'``.
        To see all the available methods, you can run `decoupler.mt.show()`_
        function.
    :type method: NoneType | str
    :param source: Column name from the provided ``net`` containing the gene
        sets to enrich for. Defaults to ``'source'``.
    :type source: str
    :param target: Column name from the provided ``net`` containing the gene set
        components (e.g. gene/protein names) that can be mapped back to the data
        set variable names. Defaults to ``'target'``.
    :type target: str
    :param weight: Defines the column in the network containing the weights to
        use in the enrichment, defaults to ``None``.
    :type weight: NoneType | str
    :param \*\*kwargs: Other keyword arguments that passed to
        `decoupler.mt.decouple()`_ function
    :type \*\*kwargs: optional

    :returns: ``None``, results are stored inplace of the passed ``data``
        object, which is a :class:`funki.input.DataSet` instance. Estimates,
        p-values and consensus scores (in case of multiple methods) are stored
        as part of the ``obsm`` attribute of the object.

    :rtype: NoneType

    .. _Decoupler: https://decoupler-py.readthedocs.io/en/latest/index.html
    .. _pandas.DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.D\
        ataFrame.html
    .. _decoupler.mt.show(): https://decoupler-py.readthedocs.io/en/latest/api/\
        mt.html
    .. _decoupler.mt.decouple(): https://decoupler-py.readthedocs.io/en/latest/\
        api/generated/decoupler.mt.decouple.html#decoupler.mt.decouple
    '''

    # Checking data
    if contrast not in data.varm_keys():

        raise KeyError(
            f'Results of the contrast {contrast} not found in the DataSet '
            'provided, please run funki.analysis.diff_exp() first or make sure '
            'it is properly written.'
        )

    # Checking method
    method = method if method in {mt.name for mt in _methods} else 'ulm'

    if not hasattr(dc.mt, method.lower()):

        raise KeyError(
            'Method not implemented in Decoupler, see `decoupler.mt.show()` '
            'for a list of available methods.'
        )

    # Storing parameters
    if 'enrich' not in data.uns['funki']:

        data.uns['funki']['enrich'] = dict()

    data.uns['funki']['enrich'][contrast] = {
        'method': method,
        'weight': weight,
        'source': source,
        'target': target,
        **kwargs,
    }

    # Preparing network matrix
    net = net.rename(columns={
        source: 'source',
        target: 'target',
        weight: 'weight'
    })

    # Making copy as AnnData to bypass Decoupler type checks
    stat = data.varm[contrast][['stat']].T.rename(index={'stat': contrast})

    score, padj = getattr(dc.mt, method.lower())(
        stat,
        net,
        raw=False,
        **kwargs
    )

    # Updating back the results to the original DataSet object
    res = pd.merge(
        score.T,
        padj.T,
        how='outer',
        left_index=True,
        right_index=True
    )
    res.columns = ['score', 'pval']

    if 'enrich' not in data.uns:

        data.uns['enrich'] = dict()

    data.uns['enrich'][contrast] = res.copy()



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

    # Storing parameters
    data.uns['funki']['sc_trans_qc_metrics'] = {
        'var_name': var_name,
    }

    aux = data.copy()
    sc.pp.calculate_qc_metrics(
        aux,
        qc_vars=[var_name],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    return DataSet(aux)


def clustering(
    data,
    alg='leiden',
    resolution=1.0,
    neigh_kwargs={},
    alg_kwargs={}
):
    '''
    Computes the clustering of the cells/samples according to the selected
    algorithm and resolution parameter. You can plot the resulting clustering
    with your dimensionality reduction approach of choice
    (:func:`funki.plots.plot_pca`, :func:`funki.plots.plot_tsne` or
    :func:`funki.plots.plot_umap`) by passing the clustering algorithm name to
    the keyword argument ``color`` (e.g.: ``plot_tsne(mydata, color='leiden')``)

    :param data: The data set instance from which to compute the clustering.
    :type data: :class:`funki.input.DataSet`
    :param alg: The algorithm to use for computing the clusters, can be either
        ``'leiden'`` or ``'louvain'``, defaults to ``'leiden'``
    :type alg: str, optional
    :param resolution: The resolution of the clustering, a higher number
        generates more and smaller clusters, defaults to ``1.0``
    :type resolution: float, optional
    :param \*\*neigh_kwargs: Other keyword arguments that are passed to
        `scanpy.pp.neighbors()`_ function, defaults to ``{}``
    :type \*\*neigh_kwargs: dict[str, any], optional
    :param \*\*alg_kwargs: Other keyword arguments that are passed to the
        clustering algorithm `scanpy.tl.louvain()`_ or `scanpy.tl.leiden()`_
        functions, defaults to ``{}``
    :type \*\*alg_kwargs: dict[str, any], optional

    :returns: ``None``, results are stored inplace of the passed ``data`` object
    :rtype: NoneType

    .. _scanpy.pp.neighbors(): https://scanpy.readthedocs.io/en/latest/api/gene\
        rated/scanpy.pp.neighbors.html
    .. _scanpy.tl.louvain(): https://scanpy.readthedocs.io/en/latest/api/genera\
        ted/scanpy.tl.louvain.html
    .. _scanpy.tl.leiden(): https://scanpy.readthedocs.io/en/latest/api/generat\
        ed/scanpy.tl.leiden.html
    '''

    # Storing parameters
    data.uns['funki']['clustering'] = {
        'alg': alg,
        'resolution': resolution,
        'neigh_kwargs': neigh_kwargs,
        'alg_kwargs': alg_kwargs,
    }

    if not 'neighbors' in data.uns:

        sc.pp.neighbors(data, **neigh_kwargs)

    if alg == 'leiden':

        sc.tl.leiden(data, resolution=resolution, **alg_kwargs)

    elif alg == 'louvain':

        sc.tl.louvain(data, resolution=resolution, **alg_kwargs)

    else:

        print('Algorithm not recognized, please use "leiden" or "louvain".')


def label_transfer(data, ref_data, transfer_label, **kwargs):
    '''
    Performs label transfer between the provided reference and target data sets.

    :param data:
        The target data set for the labels to be transferred to.
    :type data: :class:`funki.input.DataSet`
    :param ref_data:
        The reference data set to transfer the labels from. The labels to
        transfer must be present in the ``obs`` table.
    :type ref_data: :class:`funki.input.DataSet` | `anndata.AnnData`
    :param transfer_label:
        The column from the ``obs`` table which contains the labels that are to
        be transferred.
    :type transfer_label: str
    :param \*\*kwargs: Other keyword arguments that are passed to the
        label transfer function `scanpy.tl.ingest()`_.
    :type \*\*kwargs: optional

    .. _scanpy.tl.ingest(): https://scanpy.readthedocs.io/en/stable/generated/s\
        canpy.tl.ingest.html
    '''

    # Storing parameters
    data.uns['funki']['label_transfer'] = {
        'transfer_label': transfer_label,
        **kwargs
    }

    if transfer_label not in ref_data.obs:

        raise KeyError(f'{transfer_label} could not be found in reference data')

    if 'neighbors' not in ref_data.uns:

        sc.pp.neighbors(ref_data)

    # Subsetting to common genes in both reference and target data sets
    common = sorted(set(data.var_names).intersection(ref_data.var_names))
    aux = data[:, common].copy()
    ref_data = ref_data[:, common]

    sc.tl.ingest(aux, ref_data, obs=transfer_label, **kwargs)

    # Updating back the results to the original DataSet object
    data.obs[transfer_label] = aux.obs[transfer_label]


def diff_exp(
    data,
    design_factor,
    contrast_var,
    ref_var,
    method='pydeseq2',
):
    '''
    Computes differential expression analysis on the provided data based on a
    given design factor and both the contrast and reference variables (e.g.
    treatment and control).

    :param data: The data from which to compute the differential expression
    :type data: :class:`funki.input.DataSet`
    :param design_factor: Name of the column containing the variables which the
        contrasting samples are assigned. The column must be present in the
        ``data.obs`` table
    :type design_factor: str
    :param contrast_var: The variable value(s) that defines the samples that are
        to be contrasted against the reference (e.g. ``'treatment'``). The value
        must be present in the specified ``design_factor`` column
    :type contrast_var: any | list[any]
    :param ref_var: The variable value(s) that defines the refence samples (e.g.
        ``'control'``). The value must be present in the specified
        ``design_factor`` column
    :type ref_var: any | list[any]
    :param method: Which method to use for computing the differential
        expression. Available methods are ``'pydeseq2'`` or ``'limma'``,
        defaults to ``'pydeseq2'``.
    :type method: str

    :returns: ``None``, results are stored inplace of the passed ``data`` object
    :rtype: NoneType
    '''

    contrast_name = f'{contrast_var}_vs_{ref_var}'

    # Storing parameters
    if not 'diff_exp' in data.uns['funki']:

        data.uns['funki']['diff_exp'] = dict()

    data.uns['funki']['diff_exp'][contrast_name] = {
        'method': method,
        'design_factor': design_factor,
        'contrast_var': contrast_var,
        'ref_var': ref_var,
    }

    if design_factor not in data.obs_keys():

        msg = f'Design factor {design_factor} not found in provided DataSet'

        raise KeyError(msg)

    if method == 'pydeseq2':

        result = _dex_pydeseq2(data, design_factor, contrast_var, ref_var)

    elif method == 'limma':

        result = _dex_limma(data, design_factor, contrast_var, ref_var)

    # Adding results to DataSet.varm table
    data.varm[contrast_name] = result.loc[data.var_names, :]


def _dex_pydeseq2(data, design_factor, contrast_var, ref_var, n_cpus=8):

    # Converting to list if not already
    ref = ref_var if type(ref_var) is list else [ref_var]
    contrast = contrast_var if type(contrast_var) is list else [contrast_var]

    if not all(x in data.obs[design_factor].values for x in ref + contrast):

        msg = 'Contrast and/or reference value(s) not found in design factor'
        raise ValueError(msg)

    # Using PyDESeq2 to calculate differential expression
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        adata=data.copy(),
        design_factors=design_factor,
        ref_level=[design_factor] + ref,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    result = DeseqStats(
        dds,
        contrast=[design_factor] + contrast + ref,
        inference=inference
    )
    result.summary()

    res = result.results_df[[
        'baseMean',
        'log2FoldChange',
        'stat',
        'pvalue',
        'padj'
    ]]

    return res


def _dex_limma(data, design_factor, contrast_var, ref_var):

    def convert(df):

        with (ro.default_converter + pandas2ri.converter).context():

            conv = getattr(
                ro.conversion.get_conversion(),
                'rpy2py' if isinstance(df, ro.vectors.DataFrame) else 'py2rpy'
            )

            return conv(df)

    # Creating design matrix
    groups = [contrast_var, ref_var]#sorted(set(data.obs[design_factor]))

    design = pd.DataFrame(0, index=data.obs_names, columns=groups)

    for g in groups:

        design.loc[data.obs[design_factor] == g, g] = 1

    # Creating contrast matrix
    contrast = pd.DataFrame(
        [1, -1],
        index=[contrast_var, ref_var],
        columns=[f'{contrast_var}_vs_{ref_var}']
    )

    # Running limma
    base = importr('base')
    limma = importr('limma')

    rdata = convert(data.to_df().T)
    rdesign = convert(design)
    rcontrast = base.sapply(convert(contrast), base.as_numeric)

    fit = limma.lmFit(rdata, rdesign)
    fit = limma.contrasts_fit(fit, rcontrast)
    fit = limma.eBayes(fit)

    result = limma.topTable(fit, coef=1, adjust='BH', number=data.n_vars)
    res = convert(result)[['AveExpr', 'logFC', 't', 'P.Value', 'adj.P.Val']]
    # Making columns consistent with PyDESeq2
    res.columns = ['baseMean', 'log2FoldChange', 'stat', 'pvalue', 'padj']

    return res
