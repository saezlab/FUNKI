import pandas as pd
import scanpy as sc
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from utils import serial_to_dataset
from utils import dataframe_to_serial
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
import funki.preprocessing as fpp
import funki.analysis as fan
import funki.plots as fpl
import funki.pipelines as fppl


tab_norm = tab_home = dcc.Tab(
    label='Filter & normalization',
    value='tab-norm',
    children=html.Div(
        children=[
            html.H1('Filtering and normalization', style=header_style),
            html.Br(),
            html.H3('Choose filters:'),
            html.Br(),
            '- Max. genes per cell: ',
            dcc.Input(
                id='max-genes',
                type='number',
                placeholder='e.g. 5000',
                value=None,
                min=0,
                style={'width': 100}
            ),
            html.Br(),
            '- Min. genes per cell: ',
            dcc.Input(
                id='min-genes',
                type='number',
                placeholder='e.g. 500',
                value=None,
                min=0,
                style={'width': 100}
            ),
            html.Br(),
            '- Max. % mitochondrial genes per cell: ',
            dcc.Input(
                id='mito-pct',
                type='number',
                placeholder='e.g. 5',
                value=None,
                min=0,
                max=100,
                style={'width': 100}
            ),
            html.Br(),
            html.Button(
                'Apply filters',
                id='apply-filter'
            ),
            dcc.Graph(id='plot-filter', style={'height': 1500})

        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

@callback(
    Output('plot-filter', 'figure'),
    Output('plot-filter', 'style'),
    Input('proc-data', 'data'),
    prevent_initial_call=True
)
def plot_filter(data):
    if data is None:
        raise PreventUpdate
    
    dset = serial_to_dataset(data)

    fig = fppl.sc_quality_control(dset)
    
    height = 800
    fig.update_layout(height=height)

    return fig, {'height': height}

@callback(
    Output('proc-data', 'data', allow_duplicate=True),
    Input('apply-filter', 'n_clicks'),
    State('raw-data', 'data'),
    State('max-genes', 'value'),
    State('min-genes', 'value'),
    State('mito-pct', 'value'),
    prevent_initial_call=True
)
def apply_filter(n_clicks, data, max_genes, min_genes, mito_pct):
    if data is None:
        raise PreventUpdate

    dset = serial_to_dataset(data)
    dset = fpp.sc_trans_filter(
        dset,
        min_genes=min_genes,
        max_genes=max_genes,
        mito_pct=mito_pct
    )

    serial = dataframe_to_serial(dset.to_df())

    return serial