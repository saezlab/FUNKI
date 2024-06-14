import os
import base64
import io

import numpy as np
import pandas as pd

from funki.input import DataSet


def parse_contents(content, filename):
    ext = os.path.splitext(filename)[-1]

    if ext not in ('.csv', '.txt', '.xlsx'):
        raise NotImplementedError(
            f"File format {ext} not supported,"
            + " please use '.csv', '.txt' or '.xlsx'"
        )

    content_type, content_string = content.split(',')

    decoded = base64.b64decode(content_string)

    if ext in ('.csv', '.txt'):
        f = io.StringIO(decoded.decode('utf-8'))
        df = pd.read_csv(f, index_col=0)
    
    elif ext == '.xlsx':
        f = io.BytesIO(decoded)
        df = pd.read_excel(f, index_col=0)

    return df

def dataframe_to_serial(df):
    return {'index': df.index, 'records': df.to_dict('records')}

def serial_to_dataframe(data):
    df = pd.DataFrame(data['records'])
    df.index = data['index']

    return df

def serial_to_dataset(data):
    df = serial_to_dataframe(data)

    kwargs = {
        k: serial_to_dataframe(data[k])
        for k in ('obs', 'var')
        if k in data.keys()
    }
    kwargs.update({
        k: {mk: np.array(mv) for mk, mv in data[k].items()}
        for k in ('obsm', 'varm', 'obsp', 'varp')
        if k in data.keys()
    })

    return DataSet(df, **kwargs)

def dataset_to_serial(dset):
    data = dataframe_to_serial(dset.to_df())
    data.update({
        k: dataframe_to_serial(getattr(dset, k))
        for k in ('obs', 'var')
    })
    data.update({
        k: {mk: mv.tolist() for mk, mv in getattr(dset, k).items()}
        for k in ('obsm', 'varm', 'obsp', 'varp')
    })

    return data
