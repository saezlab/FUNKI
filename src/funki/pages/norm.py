import pandas as pd

from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import callback
import plotly.express as px
import plotly.graph_objects as go

from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style


tab_norm = tab_home = dcc.Tab(
    label='Filter & normalization',
    value='tab-norm',
    children=html.Div(
        children=[
            html.H1('Filtering and normalization', style=header_style),
            html.Br(),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)