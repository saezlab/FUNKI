import dash
from dash import html
from dash import dcc

from utils.style import tab_style
from utils.style import tab_selected_style

tab_home = dcc.Tab(
    label='Home',
    value='tab-home',
    children=[
        html.H1('Welcome to FUNKI'),
        html.P(
            '''Welcome to FUNKI, the omics FUNctional analysis worKflows
            Interface tool. This Python package is intended to integrate
            different omic data analysis workflows including a graphical user
            interface (GUI), but also as a standalone Python package that users
            can integrate into their existing pipelines.''',
            style={
                'text-align': 'justify',
                'padding-right': 20
            }
        )
    ],
    style=tab_style,
    selected_style=tab_selected_style,
)