import os
import base64
import io

import numpy as np
import pandas as pd
from dash import html

from funki.input import DataSet
from .info_msg import msg


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
    return {
        'obs_names': df.index,
        'var_names': df.columns,
        'X': df.values.tolist()
    }

def serial_to_dataframe(data):
    df = pd.DataFrame(
        np.array(data['X']),
        index=data['obs_names'] if 'obs_names' in data.keys() else None,
        columns=data['var_names'] if 'var_names' in data.keys() else None,
    )

    return parse_types(df)

def serial_to_dataset(data):
    df = serial_to_dataframe(data)

    kwargs = {
        k: serial_to_dataframe(data[k])
        for k in ('obs', 'var')
        if k in data.keys()
    }
    kwargs.update({
        k: {
            mk: np.array(mv)
            if type(mv) is list
            else serial_to_dataframe(mv)
            for mk, mv in data[k].items()
        }
        for k in ('obsm', 'varm', 'obsp', 'varp')
        if k in data.keys()
    })

    return DataSet(df, **kwargs)

def dataset_to_serial(dset):
    data = dataframe_to_serial(dset.to_df())

    data.update({
        k: dataframe_to_serial(getattr(dset, k))
        for k in ('obs', 'var')
        if not getattr(dset, k).empty
    })

    attrs = {}

    for k in ('obsm', 'varm', 'obsp', 'varp'):
        attrs[k] = dict()

        for mk, mv in getattr(dset, k).items():
            if type(mv) is np.ndarray:
                res = mv.tolist()
            
            elif type(mv) is pd.DataFrame:
                res = dataframe_to_serial(mv)
            
            else:
                res = mv.toarray().tolist()
            
            attrs[k][mk] = res
            
    data.update(attrs)

    return data

def md_to_str(path):
    with open(path) as f:
        md = '\n \t'
        for line in f.read():
            if '\n' in line:
                md += '\n \t'

            else:
                md += line

    return md

def parse_types(df):
    result = {}
    
    for c in df.columns:
        try:
            result[c] = pd.to_numeric(df[c])

        except ValueError:
            result[c] = df[c].astype(object)

    return pd.DataFrame(result)

def info(key):
    return html.Abbr(
        '\u2753',
        title=msg[key] if key in msg.keys() else None,
        style={
            'text-decoration': 'none',
            'font-size': 12,
            'vertical-align': 'super'
        }
    )
