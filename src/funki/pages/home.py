from dash import html
from dash import dcc

from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style


# ================================== LAYOUT ================================== #

tab_home = dcc.Tab(
    label='Home',
    value='tab-home',
    children=html.Div(
        children=[
            html.H1('Welcome to FUNKI', style=header_style),
            html.Br(),
            '''Welcome to FUNKI, the omics FUNctional analysis worKflows
            Interface tool. This Python package is intended to integrate
            different omic data analysis workflows including a graphical user
            interface (GUI), but also as a standalone Python package that users
            can integrate into their existing pipelines.''',
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)