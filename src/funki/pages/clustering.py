import numpy as np
import pandas as pd
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
from funki import _colors
import funki.analysis as fan
import funki.plots as fpl


# ================================== LAYOUT ================================== #

tab_cluster = dcc.Tab(
    label='Clustering',
    value='tab-cluster',
    children=html.Div(
        children=[
            html.H1('Clustering', style=header_style),
            html.Br(),
            html.Div(
                children=[
                    '- Choose a clustering algorithm: ',
                    dcc.RadioItems(
                        id='cluster-algorithm',
                        options=[
                            {'label': 'Leiden', 'value': 'leiden'},
                            {'label': 'Louvain', 'value': 'louvain'}
                        ],
                        value='leiden'
                    ),
                    html.Br(),
                    '- Choose resolution: ',
                    html.Div(
                        dcc.Slider(
                            id='cluster-resolution',
                            min=0.0,
                            max=2.0,
                            step=0.01,
                            marks={
                                i if i % 1 else int(i): '%.1f' % i
                                for i in np.arange(0, 2.5, 0.5)
                            },
                            tooltip={
                                'always_visible': True,
                                'placement': 'top'
                            },
                            value=1.0,
                        ),
                        style={'width': 350, 'padding-top': 20}
                    ),
                    html.Button(
                        'Calculate',
                        id='apply-cluster'
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
                    html.Div(
                        '- Apply Harmony: ',
                        style={'display': 'inline-block'}
                    ),
                    html.Div(
                        dcc.Checklist( # TODO
                            id='harmony',
                            options=[{'label': '', 'value': True}]
                        ),
                        style={'display': 'inline-block'}
                    ),
                    html.Br(),
                    html.Br(),
                    '- Choose embedding: ',
                    dcc.RadioItems(
                        id='embedding',
                        options=[
                            {'label': 'PCA', 'value': 'pca'},
                            {'label': 'tSNE', 'value': 'tsne'},
                            {'label': 'UMAP', 'value': 'umap'}
                        ],
                        value='pca'
                    ),
                    html.Br(),
                    html.Div(
                        id='param-panel',
                        style={
                            'border': 'solid',
                            'border-color': _colors['teal'],
                            'background-color': _colors['aqua'],
                            'padding': 10,
                            'width': '80%'
                        }
                    ),
                    html.Br(),
                    '- Select variable to color the embedding by: ',
                    dcc.Dropdown(
                        id='color-embedding',
                        searchable=False,
                        clearable=True,
                    ),
                    html.Br(),
                    html.Button(
                        'Visualize',
                        id='apply-embedding'
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                }
            ),
            dcc.Loading(dcc.Graph(id='plot-embedding'))
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('param-panel', 'children'),
    Output('param-panel', 'hidden'),
    Input('embedding', 'value'),
    Input('data', 'data')
)
def update_param_panel(embedding, data):
    children = []

    if embedding != 'pca':
        children.extend(['Select embedding parameters: ', html.Br(), html.Br()])

    if embedding == 'tsne':
        max_per = len(data['index']) - 1 if data else 50
        min_per = 1 if max_per < 10 else 5
        step = 1 if max_per < 10 else 5

        children.extend([
            '- Perplexity',
            html.Br(),
            dcc.Slider(
                id='perplexity',
                min=min_per,
                max=max_per if max_per > 0 else 1,
                step=1,
                marks={
                    int(i): f'{i}'
                    for i in np.arange(min_per, max_per + step, step)},
                tooltip={
                    'always_visible': True,
                    'placement': 'top'
                },
                value=30 if max_per >= 30 else max_per,
            ),
        ])

    elif embedding == 'umap':
        children.extend([
            '- Minimum distance: ',
            dcc.Input(
                id='min-dist',
                type='number',
                placeholder='e.g. 0.5',
                value=0.5,
                min=0,
                style={'width': 50}
            ),
            html.Br(),
            '- Spread: ',
            dcc.Input(
                id='spread',
                type='number',
                placeholder='e.g. 1.0',
                value=1.0,
                min=0,
                style={'width': 50}
            ),
        ])

    return children, False if children else True

@callback(
    Output('data', 'data', allow_duplicate=True),
    Input('apply-cluster', 'n_clicks'),
    State('data', 'data'),
    State('cluster-algorithm', 'value'),
    State('cluster-resolution', 'value'),
    prevent_initial_call=True
)
def apply_clustering(n_clicks, data, algorithm, resolution):
    if data is None:
        raise PreventUpdate
    
    dset = serial_to_dataset(data)
    fan.sc_clustering(dset, alg=algorithm, resolution=resolution)
    
    return dataset_to_serial(dset)

@callback(
    Output('color-embedding', 'options'),
    Input('data', 'data')
)
def update_dropdown(data):
    try:
        options = list(data['obs']['records'][0].keys())

    except (KeyError, TypeError):
        options = []

    return options

@callback(
    Output('plot-embedding', 'figure'),
    Input('apply-embedding', 'n_clicks'),
    State('data', 'data'),
    State('embedding', 'value'),
    State('param-panel', 'children'),
    State('color-embedding', 'value'),
    prevent_initial_call=True,
)
def plot_embedding(n_clicks, data, embedding, param_panel, color):
    if data is None:
        raise PreventUpdate
    
    dset = serial_to_dataset(data)

    if embedding == 'pca':
        fig = fpl.plot_pca(dset, color=color)

    elif embedding == 'tsne':
        # TODO: There is probably a more elegant way to do this
        perplexity = param_panel[-1]['props']['value']
        print(perplexity)
        fig = fpl.plot_tsne(dset, perplexity=perplexity, color=color)

    elif embedding == 'umap':
        min_dist = param_panel[-4]['props']['value']
        spread = param_panel[-1]['props']['value']

        fig = fpl.plot_umap(dset, min_dist=min_dist, spread=spread, color=color)

    return fig
