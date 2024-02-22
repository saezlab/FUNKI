import os

import anndata


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
    :type X: anndata.AnnData, optional

    .. _anndata.AnnData: https://anndata.readthedocs.io/en/latest/generated/ann\
        data.AnnData.html
    '''

    def __init__(self, X=None):
        super(self.__class__, self).__init__(X)


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

    .. _anndata reading: https://anndata.readthedocs.io/en/latest/api.htm\
        l\#reading
    .. _anndata.read_text: https://anndata.readthedocs.io/en/latest/generated/a\
        nndata.read_text.html
    '''

    fname, ext = os.path.splitext(path)

    if ext in _read_ext.keys():
        return DataSet(_read_ext[ext](path, *args, **kwargs))

    else:
        return DataSet(_read_ext['.txt'](path, *args, **kwargs))
