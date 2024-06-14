from dash import html
from dash import dcc

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
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)