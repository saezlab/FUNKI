import numpy as np
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
import funki.analysis as fan
import funki.preprocessing as fpp
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
                    html.H3(
                        children=[
                            'Embedding:',
                            info('embedding')
                        ]
                    ),
                    html.Div(
                        '- Apply Harmony: ',
                        style={'display': 'inline-block'}
                    ),
                    html.Div(
                        dcc.Checklist(
                            id='harmony',
                            options=[{'label': '', 'value': True}]
                        ),
                        style={'display': 'inline-block'}
                    ),
                    html.Br(),
                    html.Div(
                        children=[
                            '- Select variable(s) to correct for:',
                            html.Br(),
                            dcc.Dropdown(
                                id='harmonize-var',
                                multi=True,
                                style={'width': '80%'}
                            ),
                        ],
                        id='harmonize-var-panel',
                        hidden=True,
                    ),
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
                    # tSNE params
                    html.Div(
                        id='param-panel-tsne',
                        children=[
                            'Select embedding parameters: ',
                            html.Br(),
                            html.Br(),
                            '- Perplexity',
                            html.Br(),
                            dcc.Slider(
                                id='perplexity',
                                min=0,
                                max=1,
                                tooltip={
                                    'always_visible': True,
                                    'placement': 'top'
                                },
                            ),
                        ],
                        style={
                            'border': 'solid',
                            'border-color': _colors['teal'],
                            'background-color': _colors['aqua'],
                            'padding': 10,
                            'width': '80%'
                        },
                        hidden=True
                    ),
                    # UMAP params
                    html.Div(
                        id='param-panel-umap',
                        children=[
                            'Select embedding parameters: ',
                            html.Br(),
                            html.Br(),
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
                        ],
                        style={
                            'border': 'solid',
                            'border-color': _colors['teal'],
                            'background-color': _colors['aqua'],
                            'padding': 10,
                            'width': '80%'
                        },
                        hidden=True
                    ),
                    html.Br(),
                    '- Select variable to color the embedding by: ',
                    dcc.Dropdown(
                        id='color-embedding',
                        searchable=False,
                        clearable=True,
                        style={'width': '80%'}
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
            html.Div(
                children=[
                    html.H3(
                        children=[
                            'Clustering:',
                            info('clustering')
                        ]
                    ),
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
                    html.Br(),
                    html.Br(),
                    dcc.Loading(
                        html.Div(
                            id='cluster-complete',
                            children='Clustering computed successfully!',
                            hidden=True,
                        ),
                        color=_colors['teal']
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                }
            ),
            dcc.Loading(
                dcc.Graph(id='plot-embedding'),
                color=_colors['teal']
            ),
            html.Button(
                'Open plot in new tab',
                id='nw-plot-embedding',
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('param-panel-umap', 'hidden'),
    Output('param-panel-tsne', 'hidden'),
    Output('perplexity', 'min'),
    Output('perplexity', 'max'),
    Output('perplexity', 'step'),
    Output('perplexity', 'marks'),
    Output('perplexity', 'value'),
    Input('embedding', 'value'),
    Input('data', 'data')
)
def update_param_panel(embedding, data):
    if not data:
        raise PreventUpdate
    
    # Fallback defaults
    umap_hidden = True
    tsne_hidden = True
    per_min = 5
    per_max = 50
    per_step = 5
    per_marks = None
    per_value = 30

    if embedding == 'tsne':
        tsne_hidden = False
        per_max = len(data['obs_names']) - 1 if len(data['X']) < 50 else 50
        per_max = per_max if per_max > 0 else 1
        per_min = 1 if per_max < 10 else 5
        per_step = 1 if per_max < 10 else 5
        per_marks = {
            int(i): f'{i}'
            for i in np.arange(per_min, per_max + per_step, per_step)
        }
        per_value = 30 if per_max >= 30 else per_max
    

    elif embedding == 'umap':
        umap_hidden = False

    return (
        umap_hidden,
        tsne_hidden,
        per_min,
        per_max,
        per_step,
        per_marks,
        per_value
    )

@callback(
    Output('data', 'data', allow_duplicate=True),
    Output('cluster-complete', 'hidden'),
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
    
    return dataset_to_serial(dset), False

@callback(
    Output('color-embedding', 'options'),
    Input('data', 'data')
)
def update_dropdown_embedding(data):
    try:
        options = list(data['obs']['var_names'])

    except (KeyError, TypeError):
        options = []

    return options

@callback(
    Output('plot-embedding', 'figure'),
    Input('apply-embedding', 'n_clicks'),
    State('data', 'data'),
    State('embedding', 'value'),
    State('perplexity', 'value'),
    State('min-dist', 'value'),
    State('spread', 'value'),
    State('color-embedding', 'value'),
    State('harmony', 'value'),
    State('harmonize-var', 'value'),
    prevent_initial_call=True,
)
def plot_embedding(
    n_clicks,
    data,
    embed,
    perplexity,
    min_dist,
    spread,
    color,
    harmony,
    hvar
):
    if data is None:
        raise PreventUpdate
    
    dset = serial_to_dataset(data)
    print(dset)
    
    # Run Harmony?
    if harmony and hvar:
        fpp.harmonize(
            dset,
            hvar if type(hvar) is list else [hvar],
            recalculate=True
        )
    
    # Bypassing clusters being considered numerical
    if color in ('louvain', 'leiden'):
        dset.obs[color] = dset.obs[color].astype('category')

    # Embeddings
    if embed == 'pca':
        fig = fpl.plot_pca(dset, color=color)

    elif embed == 'tsne':
        fig = fpl.plot_tsne(dset, perplexity=perplexity, color=color)

    elif embed == 'umap':
        fig = fpl.plot_umap(dset, min_dist=min_dist, spread=spread, color=color)

    return fig

@callback(
    Output('harmonize-var-panel', 'hidden'),
    Input('harmony', 'value')
)
def update_harmony_var_panel(harmony):
    return not harmony

@callback(
    Output('harmonize-var', 'options'),
    Input('data', 'data')
)
def update_dropdown_embedding(data):
    try:
        options = list(data['obs']['var_names'])

    except (KeyError, TypeError):
        options = []

    return options

@callback(
    Input('nw-plot-embedding', 'n_clicks'),
    State('plot-embedding', 'figure')
)
def plot_embedding_new_tab(n_clicks, fig):
    if fig:
        return go.Figure(fig).show()