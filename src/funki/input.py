import os

import anndata

# TODO: Implement concat function returning DataSet
# TODO: Re-implement __repr__

_read_ext = {
    '.csv': anndata.read_csv,
    '.h5ad': anndata.read_h5ad,
    '.xlsx': anndata.read_excel,
    '.h5': anndata.read_hdf,
    '.loom': anndata.read_loom,
    '.mtx': anndata.read_mtx,
    '.gz': anndata.read_umi_tools,
    '.txt': anndata.read_text, # Wildcard
}


class DataSet(anndata.AnnData):
    '''
    Data set class, inherits from `anndata.AnnData`_ class.

    :param X: An instance of `anndata.AnnData`_ or any other class accepted by
        `anndata.AnnData`_, defaults to ``None``.
    :type X: `anndata.AnnData`_, optional
    :param \*\*kwargs: Other keyword arguments that can be passed to
        `anndata.AnnData`_ class.
    :type \*\*kwargs: optional

    .. _anndata.AnnData: https://anndata.readthedocs.io/en/latest/generated/ann\
        data.AnnData.html
    '''

    def __init__(self, X=None, **kwargs):
        super().__init__(X, **kwargs)

    def __getitem__(self, index):
        '''
        Returns a sliced view of the DataSet.
        '''

        oidx, vidx = self._normalize_indices(index)

        X = self.X[oidx, vidx]
        obs = self.obs.iloc[oidx]
        var = self.var.iloc[vidx]
        uns = self.uns
        obsm = {k: v[odsx] for k, v in self.obsm.items()}
        varm = {k: v[vidx] for k, v in self.varm.items()}
        raw = self.raw
        layers = self.layers
        obsp = {k: v[odsx] for k, v in self.obsp.items()}
        varp = {k: v[vidx] for k, v in self.varp.items()}

        return DataSet(
            X=X,
            obs=obs,
            var=var,
            uns=uns,
            obsm=obsm,
            varm=varm,
            raw=raw,
            layers=layers,
            obsp=obsp,
            varp=varp,
        )

    def _del_meta(self, attrs):
        '''
        Deletes metadata table(s) from a DataSet object attributes. If the
        attribute keys are not found, captures the ``KeyError`` and moves on.

        :param attrs: The attributes and corresponding metadata table(s) to
            delete. Keys should correspond to the attribute name (``'obsm'``,
            ``'obsp'``, ``'varm'``, ``'varp'`` or ``'uns'``) and values should
            correspond to the key(s) in that attribute
        :type attrs: dict[str, str | list[str]]
        '''

        for attr, keys in attrs.items():
            if isinstance(keys, str):
                keys = [keys]

            for key in keys:
                try:
                    del(self.__dict__['_%s' % attr][key])

                except KeyError:
                    continue

    def copy(self):
        '''
        Creates a copy of the current :class:`DataSet` instance.

        :returns: The copy of the current object
        :rtype: :class:`DataSet`
        '''

        return DataSet(self._mutated_copy())

    def sum_duplicated_gene_counts(self):
        '''
        Takes duplicated genes (as defined in :attr:`DataSet.var_names`) and
        sums their counts. NOTE: This should be applied to raw counts, prior to
        any normalization or preprocessing as is it assumed that e.g. different
        splicing variants of the same gene can be added  towards the counts of
        the same gene (e.g. gene symbol, UniProt ID, etc.).
        '''

        aux = DataSet(self.to_df().groupby(self.var_names, axis=1).sum())
        self.__dict__.update(aux.__dict__)


def read(path, *args, **kwargs):
    '''
    Function takes a path to a file and uses the appropiate `anndata reading`_
    function according to the file extension provided in ``path``. You can find
    further information on the available formats in the
    `anndata reading`_ documentation. If file extension not recognized, it
    defaults to `anndata.read_text`_.

    :param path: Path to the data file.
    :type path: str
    :param \*args: Other positional arguments that can be required by the
        reading function, according to the data file format.
    :type \*args: optional
    :param \*\*kwargs: Other keyword arguments that can be required by the
        reading function, according to the data file format.
    :type \*\*kwargs: optional
    :returns: The loaded data set object
    :rtype: :class:`funki.input.DataSet`

    .. _anndata reading: https://anndata.readthedocs.io/en/latest/api.html/#rea\
        ding
    .. _anndata.read_text: https://anndata.readthedocs.io/en/latest/generated/a\
        nndata.read_text.html
    '''

    fname, ext = os.path.splitext(path)

    if ext in _read_ext.keys():
        return DataSet(_read_ext[ext](path, *args, **kwargs))

    else:
        return DataSet(_read_ext['.txt'](path, *args, **kwargs))
