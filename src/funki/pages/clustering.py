from dash import html
from dash import dcc
import numpy as np

from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style


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
                        dcc.Checklist(
                            id='harmony',
                            options=[{'label': '', 'value': True}]
                        ),
                        style={'display': 'inline-block'}
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                }
            ),
            html.Button(
                'Calculate',
                id='apply-cluster'
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)