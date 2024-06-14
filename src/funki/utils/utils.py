import os
import base64
import io

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

    del data['index'], data['records']

    kwargs = {k: serial_to_dataframe(v) for k, v in data.items()}

    return DataSet(df, **kwargs)

def dataset_to_serial(dset):
    data = dataframe_to_serial(dset.to_df())
    data.update({
        k: dataframe_to_serial(getattr(dset, k))
        for k in ('obs', 'var')
    })

    return data
