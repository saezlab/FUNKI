from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import pandas as pd

from utils import serial_to_dataset
from utils import dataset_to_serial
from utils import serial_to_dataframe
from utils import info
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
from funki import _colors


# ================================== LAYOUT ================================== #

tab_difexp = dcc.Tab(
    label='Differential expression',
    value='tab-difexp',
    children=html.Div(
        children=[
            html.H1('Differential expression analysis', style=header_style),
            html.Br(),
            'Choose group variable for the contrast:',
            dcc.Dropdown(
                id='obs-selector',
                clearable=True,
                multi=False,
                searchable=True,
                style={'width': '80%'}
            ),
            html.Div(
                children=[
                    'Choose group(s) of samples to contrast from '
                    '(e.g.: treatment):',
                    dcc.Dropdown(
                        id='group-selector-a',
                        clearable=True,
                        multi=True,
                        searchable=True,
                        style={'width': '80%'}
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                    'pad': 10
                }
            ),
            html.Div(
                children=[
                    'Choose group(s) of samples to contrast against '
                    '(e.g. control):',
                    dcc.Dropdown(
                        id='group-selector-b',
                        clearable=True,
                        multi=True,
                        searchable=True,
                        style={'width': '80%'}
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                    'pad': 10
                }
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('obs-selector', 'options'),
    Input('data', 'data')
)
def update_obs_selector(data):
    if data is None:
        raise PreventUpdate

    elif 'obs' not in data.keys():
        raise PreventUpdate

    df = serial_to_dataframe(data['obs'])
    options = [
        c for c in df.columns
        if not pd.api.types.is_numeric_dtype(df[c])
        or c in ('leiden', 'louvain')
    ]

    return options

@callback(
    Output('group-selector-a', 'options'),
    Output('group-selector-b', 'options'),
    Input('data', 'data'),
    Input('obs-selector', 'value'),
    Input('group-selector-a', 'value'),
    Input('group-selector-b', 'value'),
    prevent_initial_call=True,
)
def update_group_selector(data, obs_var, va, vb):
    if data is None:
        raise PreventUpdate

    elif 'obs' not in data.keys():
        raise PreventUpdate
    
    elif not obs_var:
        return [], []

    df = serial_to_dataframe(data['obs'])
    options = list(set(df[obs_var]))

    options_a = [o for o in options if o not in vb] if vb else options
    options_b = [o for o in options if o not in va] if va else options

    return options_a, options_b

