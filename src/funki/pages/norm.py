from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate

from utils import serial_to_dataset
from utils import dataset_to_serial
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
import funki.preprocessing as fpp
import funki.pipelines as fppl


# ================================== LAYOUT ================================== #

tab_norm = tab_home = dcc.Tab(
    label='Filter & normalization',
    value='tab-norm',
    children=html.Div(
        children=[
            html.H1('Filtering and normalization', style=header_style),
            html.Br(),
            html.Div(
                children=[
                    html.H3('Choose filters:'),
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
                    '- Max. % mito. genes per cell: ',
                    dcc.Input( # TODO: Switch to slider?
                        id='mito-pct',
                        type='number',
                        placeholder='e.g. 5',
                        value=None,
                        min=0,
                        max=100,
                        style={'width': 100}
                    ),
                    html.Br(),
                    html.Br(),
                    html.Button(
                        'Apply filters',
                        id='apply-filter'
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
                    html.H3('Choose normalization:'),
                    '- Size factor: ',
                    dcc.Input(
                        id='size-factor',
                        type='number',
                        placeholder='e.g. 1000000',
                        value=None,
                        min=0,
                        style={'width': 100}
                    ),
                    html.Br(),
                    html.Div(
                        '- Log-transform: ',
                        style={'display': 'inline-block'}
                    ),
                    html.Div(
                        dcc.Checklist(
                            id='log-transform',
                            options=[{'label': '', 'value': True}]
                        ),
                        style={'display': 'inline-block'}
                    ),
                    html.Br(),
                    html.Br(),
                    html.Button(
                        'Apply normalization',
                        id='apply-norm'
                    )
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                    'pad': 10
                }
            ),
            dcc.Graph(id='plot-filter', style={'height': 1500}),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('plot-filter', 'figure'),
    Output('plot-filter', 'style'),
    Input('data', 'data'),
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
    Output('data', 'data', allow_duplicate=True),
    Input('apply-filter', 'n_clicks'),
    State('data', 'data'),
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

    serial = dataset_to_serial(dset)

    return serial

@callback(
    Output('data', 'data', allow_duplicate=True),
    Input('apply-norm', 'n_clicks'),
    State('data', 'data'),
    State('size-factor', 'value'),
    State('log-transform', 'value'),
    prevent_initial_call=True
)
def apply_norm(n_clicks, data, size_factor, log_transform):
    if data is None:
        raise PreventUpdate

    dset = serial_to_dataset(data)
    dset = fpp.sc_trans_normalize_total(
        dset,
        target_sum=size_factor,
        log_transform=log_transform
    )

    serial = dataset_to_serial(dset)

    return serial