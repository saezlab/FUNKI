from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go 

from utils import serial_to_dataset
from utils import dataset_to_serial
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
            html.Div(
                children=[
                    'Choose samples to contrast from (e.g. treatment):',
                    dcc.Dropdown(
                        id='samples-a',
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
                    'Choose samples to contrast against (e.g. control):',
                    dcc.Dropdown(
                        id='samples-b',
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
    Output('samples-a', 'options'),
    Output('samples-b', 'options'),
    Input('data', 'data')
)
def update_dropdown_samples(data):
    try:
        options = list(data['obs_names'])

    except (KeyError, TypeError):
        options = []

    return options, options
