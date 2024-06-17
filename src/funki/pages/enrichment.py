import numpy as np
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
                            '- Select variable to exclude gene sets from: ',
                            dcc.Dropdown(
                                id='gset-exclude-from-col',
                                searchable=True,
                                clearable=True,
                                disabled=True,
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
                                    style_cell={'width': '50px'}
                                )
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
    Output('table-gset', 'columns'),
    Output('table-gset', 'data'),
    Output('gset-exclude-from-col', 'options'),
    Output('gset-exclude-from-col', 'disabled'),
    Input('gset-collection', 'value')
)
def load_gset_table(gset):
    if gset is None:
        return None, None, [], True

    df = dc.get_resource(gset)

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    options = [
        {'label': c, 'value': c}
        for c in df.columns
        if c != 'genesymbol'
    ]

    return table_columns, table_data, options, False