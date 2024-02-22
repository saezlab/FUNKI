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
    def __init__(self, X=None):
        super(self.__class__, self).__init__(X)


def read(path, *args, **kwargs):
    '''
    Function takes a path to a file and uses the appropiate `anndata` reading
    function according to the file extension.
    '''

    fname, ext = os.path.splitext(path)

    if ext in _read_ext.keys():
        return DataSet(_read_ext[ext](path, *args, **kwargs))

    else:
        return DataSet(_read_ext['.txt'](path, *args, **kwargs))
