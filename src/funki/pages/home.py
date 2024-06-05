import dash
from dash import html
from dash import dcc

tab_home = dcc.Tab(
    label='Home',
    value='tab-home',
    children=[
        html.H1('Home'),
        html.Div('Lorem ipsum')
    ]
)