import scanpy as sc

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
    '''

    aux = data.copy()

    if min_genes:
        sc.pp.filter_cells(aux, min_genes=min_genes, inplace=True)

    if max_genes:
        sc.pp.filter_cells(aux, max_genes=max_genes, inplace=True)

    if mito_pct:
        aux.var['mt'] = aux.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(aux, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)
        aux =  aux[aux.obs.pct_counts_mt < mito_pct, :]

    return aux
