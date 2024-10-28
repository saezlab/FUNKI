import scanpy as sc
import harmonypy as hm

from .analysis import sc_trans_qc_metrics
from .input import DataSet

def sc_trans_filter(data, min_genes=None, max_genes=None, mito_pct=None):
    '''
    Applies quality control filters to a given single-cell transcriptomic data
    set. Can filter out cells based on a minimum and maximum number of genes as
    well as based on the percentage of mitochondrial genes.

    :param data: A single-cell transcriptomic data set containing raw counts
    :type data: :class:`funki.input.DataSet`
    :param min_genes: Minimum number of different genes for a cell. If the
        number is below, that cell will be filtered out, defaults to ``None``
    :type min_genes: int, optional
    :param max_genes: Maximum number of different genes for a cell. If the
        number is above, that cell will be filtered out, defaults to ``None``
    :type max_genes: int, optional
    :param mito_pct: Percentage of mitochondrial genes. If a cell has a
        percentage above the threshold, it will be filtered out, defaults to
        ``None``
    :type mito_pct: int, optional
    :returns: The resulting filtered data set after applying the specified
        thresholds
    :rtype: :class:`funki.input.DataSet`
    '''

    # Storing parameters
    data.uns['funki']['sc_trans_filter'] = {
        'min_genes': min_genes,
        'max_genes': max_genes,
        'mito_pct': mito_pct,
    }

    aux = data.copy()

    if min_genes:
        sc.pp.filter_cells(aux, min_genes=min_genes, inplace=True)

    if max_genes:
        sc.pp.filter_cells(aux, max_genes=max_genes, inplace=True)

    if mito_pct:
        aux.var['mito'] = aux.var_names.str.upper().str.startswith('MT-')
        aux = sc_trans_qc_metrics(aux, var_name='mito')
        aux = aux[aux.obs['pct_counts_mito'] < mito_pct, :]

    return DataSet(aux)

def sc_trans_normalize_total(data, target_sum=1e6, log_transform=False):
    '''
    Normalizes the total counts per cell in a single-cell data set. The
    normaliztion scales the counts so that the sum of all genes in a cell add up
    to the specified ``target_sum`` (1e6 by default, equivalent to CPM
    normalization). If ``log_transform=True``, it also applies a
    :math:`\log(X+1)` transformation to the resulting normailzed data.

    :param data: A single-cell transcriptomic data set containing raw counts
    :type data: :class:`funki.input.DataSet`
    :param target_sum: The targeted total counts per cell to normalize for,
        defaults to ``None``, which is equivalent to CPM normalization
    :type target_sum: int | float, optional
    :param log_transform: Whether to apply log-transformation after normalizing
        the data, defaults to ``False``
    :type log_transform: bool, optional
    :returns: The resulting normalized (and log-transformed, if applicable)
        data set
    :rtype: :class:`funki.input.DataSet`
    '''

    # Storing parameters
    data.uns['funki']['sc_trans_normalize_total'] = {
        'target_sum': target_sum,
        'log_transform': log_transform,
    }

    aux = data.copy()
    
    if target_sum:
        sc.pp.normalize_total(aux, target_sum=target_sum, inplace=True)

    if log_transform:
        sc.pp.log1p(aux)

    return DataSet(aux)

def harmonize(data, vars_use, use_highly_variable=True, recalculate=False,
              **kwargs):
    '''
    Executes `Harmony`_ batch correction on the data set PCA embedding. NOTE:
    this method will overwrite the :attr:`DataSet.obsm['X_pca']` matrix.

    :param data: The data set that which `Harmony`_ will be executed on
    :type data: :class:`funki.input.DataSet`
    :param vars_use: Variables over which to correct for (i.e. batches). Must
        correspond to column(s) defined in :attr:`DataSet.obs`
    :type vars_use: list[str]
    :param use_highly_variable: Whether to use highly variable genes only or all
        genes available. Only used if PCA has not been computed previously or if
        ``recalculate=True``, defaults to ``True``
    :type use_highly_variable: bool, optional
    :param recalculate: Whether to recalculate the PCA dimensionality reduction,
        defaults to ``False``
    :type recalculate: bool, optional
    :param \*\*kwargs: Other keyword arguments that can be passed to
        ``harmonypy.run_harmony()``
    :type \*\*kwargs: optional

    .. _Harmony: https://portals.broadinstitute.org/harmony/
    '''

    # Storing parameters
    data.uns['funki']['harmonize'] = {
        'vars_use': vars_use,
        'use_highly_variable': use_highly_variable,
        'recalculate': recalculate,
        **kwargs
    }

    if recalculate:
        data._del_meta({'obsm': 'X_pca'})

    if use_highly_variable:
        sc.pp.highly_variable_genes(data, inplace=True)

    if 'X_pca' not in data.obsm:
        sc.pp.pca(data, use_highly_variable=use_highly_variable)

    ho = hm.run_harmony(data.obsm['X_pca'], data.obs, vars_use, **kwargs)
    data.obsm['X_pca'] = ho.Z_corr.T
