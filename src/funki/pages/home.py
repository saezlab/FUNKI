from dash import html
from dash import dcc

from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
from utils import md_to_str

# ================================== LAYOUT ================================== #

tab_home = dcc.Tab(
    label='Home',
    value='tab-home',
    children=html.Div(
        children=[
            html.H1('Welcome to FUNKI', style=header_style),
            dcc.Markdown(
                md_to_str('src/funki/assets/home.md'),
                style={'width': 'auto'}
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)
