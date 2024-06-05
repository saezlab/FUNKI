import dash
from dash import html
from dash import dcc

tab_data_load = dcc.Tab(
    label='Data load',
    value='tab-data-load',
    children=[
        html.H1('Data loading'),
        html.Div('Lorem ipsum')
    ]
)