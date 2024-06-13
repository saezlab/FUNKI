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

def serial_to_dataframe(data):
    df = pd.DataFrame(data['records'])
    df.index = data['index']

    return df

def serial_to_dataset(data, annot):
    df = serial_to_dataframe(data)
    ann = serial_to_dataframe(annot)

    return DataSet(df, obs=ann)