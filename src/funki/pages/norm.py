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
from plotly.tools import mpl_to_plotly

from utils import serial_to_dataset
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
import funki.preprocessing as fpp


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
                id='pct-mito',
                type='number',
                placeholder='e.g. 5',
                value=None,
                min=0,
                max=100,
                style={'width': 100}
            ),
            html.Br(),
            html.Button(
                id='apply-filter'
            ),
            dcc.Graph(id='plot-filter')

        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

@callback(
    Output('plot-filter', 'figure'),
    Input('proc-data', 'data'),
    Input('ann-data', 'data')
)
def plot_filter(data, annot):
    if None in (data, annot):
        raise PreventUpdate
    
    dset = serial_to_dataset(data, annot)
    print(dset) # TODO