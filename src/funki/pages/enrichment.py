import numpy as np
import pandas as pd
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
from dash.dash_table import DataTable
import decoupler as dc

from utils import serial_to_dataset
from utils import dataset_to_serial
from utils import serial_to_dataframe
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
from funki import _colors
import funki.analysis as fan
import funki.plots as fpl


# ================================== LAYOUT ================================== #

tab_enrichment = dcc.Tab(
    label='Enrichment',
    value='tab-enrichment',
    children=html.Div(
        children=[
            html.H1('Clustering', style=header_style),
            html.Br(),
            html.Div(
                children=[
                    html.Div(
                        children=[
                            'Choose a gene set collection: ',
                            html.Br(),
                            dcc.Dropdown(
                                id='gset-collection',
                                options=[
                                    {'label': i, 'value': i}
                                    for i in dc.show_resources()
                                ],
                                searchable=True,
                                clearable=True,
                            ),
                            html.Br(),
                            dcc.Loading(
                                DataTable(
                                    id='table-gset',
                                    fixed_rows={'headers': True, 'data': 0},
                                    fixed_columns={'headers': True, 'data': 1},
                                    style_table={
                                        'maxHeight': 500,
                                        'minWidth': '100%',
                                        'overflowY': 'auto',
                                        'overflowX': 'auto'
                                    },
                                    style_cell={
                                        'width': 100,
                                        'whiteSpace': 'normal'
                                    }
                                )
                            ),
                            html.Br(),
                            html.Div(
                                id='gset-excl-from-col',
                                hidden=True,
                                children=[
                                    '- Select variable to filter by:',
                                    dcc.Dropdown(
                                        id='gset-excl-from-col-select',
                                        searchable=True,
                                        clearable=True,
                                    ),
                                    html.Br(),
                                    html.Div(
                                        id='gset-excl-elems-num',
                                        hidden=True,
                                        children=[
                                            '- Select range of values to keep:',
                                            dcc.RangeSlider(
                                                id='gset-excl-elems-num-select',
                                                min=0,
                                                max=1,
                                                tooltip={
                                                    'always_visible': True,
                                                    'placement': 'top'
                                                },
                                            ),

                                        ]
                                    ),
                                    html.Div(
                                        id='gset-excl-elems-cat',
                                        hidden=True,
                                        children=[
                                            '- Select variable(s) to exclude:',
                                            dcc.Dropdown(
                                                id='gset-excl-elems-cat-select',
                                                searchable=True,
                                                clearable=True,
                                                multi=True,
                                            ),
                                        ]
                                    ),
                                    html.Br(),
                                    html.Button(
                                        'Apply filter',
                                        id='apply-gset-filter',
                                        disabled=True
                                    ),
                                ],
                            ),
                        ],
                        style={
                            'width': '49%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
                        }
                    ),
                    html.Div(
                        children=[

                        ],
                        style={
                            'width': '49%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
                        }                    
                    )
                ]
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('table-gset', 'columns', allow_duplicate=True),
    Output('table-gset', 'data', allow_duplicate=True),
    Output('gset-excl-from-col-select', 'options'),
    Output('gset-excl-from-col', 'hidden'),
    Input('gset-collection', 'value'),
    prevent_initial_call=True
)
def load_gset_table(gset):
    if gset is None:
        return None, None, [], True

    df = dc.get_resource(gset)

    if len(df) == 0:
        return None, [{0: 'Error downloading the data'}], [], True

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    options = [
        {'label': c, 'value': c}
        for c in df.columns
        if c != 'genesymbol'
    ]

    return table_columns, table_data, options, False

@callback(
    Output('gset-excl-elems-num-select', 'min'),
    Output('gset-excl-elems-num-select', 'max'),
    Output('gset-excl-elems-num', 'hidden'),
    Output('gset-excl-elems-cat-select', 'options'),
    Output('gset-excl-elems-cat', 'hidden'),
    Output('apply-gset-filter', 'disabled'),
    Input('gset-excl-from-col-select', 'value'),
    State('table-gset', 'data'),
    prevent_initial_call=True
)
def update_filter(col, data):    
    df = pd.DataFrame(data)
    
    # Fallback defaults
    min, max = 0, 1
    options = []
    hid_num, hid_cat = True, True
    dis_button = True

    if not col:
        pass

    elif pd.api.types.is_numeric_dtype(df[col]) and df[col].dtype!= bool:
        max, min = df[col].max(), df[col].min()
        hid_num = False
        hid_cat = True
        dis_button = False

    else:
        options = [
            {'label': str(i), 'value': str(i)}
            for i in df[col].unique()
        ]
        hid_cat = False
        hid_num = True
        dis_button = False

    return min, max, hid_num, options, hid_cat, dis_button
